"""
build_catalog.py — genera il catalogo metadati delle serie temporali (passo 1).

Input : file Excel grezzi (uno per anno/strumento) in INPUT_DIR
Output: catalog.json (stazioni, strumenti, serie, file sorgente)
        series.csv   (vista piatta delle serie, editabile in Excel)

RIFERIMENTO TEMPORALE: UTC+1 (ora solare italiana, CET, senza ora legale)
per tutte le misure. I timestamp sono scritti in ISO 8601 con offset "+01:00".

Le informazioni "curate" (stazioni, strumenti, regole di riconoscimento
delle colonne, note di qualità) sono definite qui sotto; periodo, passo,
lacune e statistiche sono calcolati automaticamente dai dati.
"""
import json, re, glob, os, math
from datetime import datetime, timezone, timedelta
import pandas as pd

INPUT_DIR = os.environ.get("INPUT_DIR", ".")
OUT_DIR = os.environ.get("OUT_DIR", ".")
TZ = timezone(timedelta(hours=1))           # UTC+1 fisso
GAP_REPORT_H = 24                           # lacune elencate nel catalogo se > 24 h

def iso(ts):
    return ts.to_pydatetime().replace(tzinfo=TZ).isoformat(timespec="minutes")

# ----------------------------------------------------------------------------
# 1. STAZIONI (siti di monitoraggio)
# ----------------------------------------------------------------------------
STATIONS = {}

# Stazioni Salerno Meter (Università di Salerno): logger METER con sensori TEROS 10 / TEROS 21
#   codice: (nome breve, ID logger, lat, lon, quota m s.l.m.)
METER_STATIONS = {
    "NV": ("Meter NV", "Z6-28531", 40.6930886, 14.7605388, 247),
    "V":  ("Meter V",  "Z6-28525", 40.6931392, 14.7605476, 255),
}
for _k, (_n, _id, _la, _lo, _el) in METER_STATIONS.items():
    STATIONS[f"SAL_{_k}"] = {
        "name": f"Salerno – {_n}", "site": "Salerno", "municipality": "Salerno (SA)",
        "lat": _la, "lon": _lo, "elevation_m": _el,
        "owner": "Università di Salerno", "external_code": _id,
        "notes": f"Datalogger METER (ID {_id}) con sensori TEROS 10 (contenuto d'acqua) e TEROS 21 "
                 "(potenziale matriciale e temperatura del terreno) a tre profondità.",
    }

# Pluviometri della rete del Centro Funzionale Multirischi (Regione Campania).
# Coordinate fornite in UTM WGS84 fuso 33N (EPSG:32633), convertite in WGS84 geografiche (pyproj).
def _utm_to_wgs84(e, n):
    from pyproj import Transformer
    lon, lat = Transformer.from_crs(32633, 4326, always_xy=True).transform(e, n)
    return round(lat, 6), round(lon, 6)

def _cf(code, name, municipality, lat, lon, elev, utm):
    return {"name": f"{name} (CF {code})", "site": municipality.split(" (")[0], "municipality": municipality,
            "lat": lat, "lon": lon, "elevation_m": elev,
            "owner": "Centro Funzionale Multirischi – Regione Campania", "external_code": code,
            "notes": "Stazione della rete regionale; dati forniti con sola colonna 'Time UTC+1'. "
                     f"Coordinate originali UTM 33N: E {utm[0]}, N {utm[1]}."}

CF_STATIONS = {   # codice: (nome, comune, (E, N) UTM33N, quota m s.l.m., pattern file)
    "21521": ("Cologna",         "Cologna (Pellezzano, SA)", (481135, 4507808), 125, r"^COLOGNA_PLUVIOMETRO_CF_\d{4}\.xlsx$"),
    "18925": ("Pellezzano",      "Pellezzano (SA)",          (479507, 4508839), 356, r"^PELLEZZANO_PLUVIOMETRO_CF_\d{4}\.xlsx$"),
    "18957": ("Salerno Genio Civile", "Salerno (SA)",        (479043, 4503222), 28,  r"^SALERNO_GENIO_CIVILE_PLUVIOMETRO_CF_\d{4}\.xlsx$"),
}
for _code, (_n, _m, _utm, _el, _pat) in CF_STATIONS.items():
    _la, _lo = _utm_to_wgs84(*_utm)
    STATIONS[f"CF_{_code}"] = _cf(_code, _n, _m, _la, _lo, _el, _utm)


# ----------------------------------------------------------------------------
# 2. STRUMENTI (logger / sensori fisici) e riconoscimento file
# ----------------------------------------------------------------------------
INSTRUMENTS = {}
for _k, (_n, _id, *_r) in METER_STATIONS.items():
    INSTRUMENTS[f"SAL_{_k}_METER"] = {
        "station_id": f"SAL_{_k}", "manufacturer": "METER Group",
        "type": "Datalogger con sensori TEROS 10 (contenuto d'acqua) e TEROS 21 (potenziale, temperatura)",
        "file_pattern": rf"^Salerno_Meter_{_k}_\d{{4}}\.xlsx$"}
for _code, (_n, _m, _utm, _el, _pat) in CF_STATIONS.items():
    INSTRUMENTS[f"CF_{_code}_RAIN"] = {"station_id": f"CF_{_code}", "manufacturer": "Centro Funzionale Regione Campania",
                                      "type": "Pluviometro", "file_pattern": _pat}


# Stazioni/strumenti "griglia" (pixel ERA5-Land, SMAP L4) creati da ingest_gee.py: non derivano da file Excel,
# si conservano dal catalogo precedente (REBUILD=1 li elimina; poi rieseguire ingest_gee.py).
_prev_cat = os.path.join(OUT_DIR, "catalog.json")
if os.path.exists(_prev_cat) and os.environ.get("REBUILD") != "1":
    _pc = json.load(open(_prev_cat, encoding="utf-8"))
    for _s in _pc.get("stations", []):
        if _s.get("kind") == "grid":
            STATIONS[_s["station_id"]] = {k: v for k, v in _s.items() if k not in ("station_id", "distance_km")}
    for _i in _pc.get("instruments", []):
        if _i["station_id"] in STATIONS and STATIONS[_i["station_id"]].get("kind") == "grid":
            INSTRUMENTS[_i["instrument_id"]] = {k: v for k, v in _i.items() if k != "instrument_id"}

# ----------------------------------------------------------------------------
# 3. GRANDEZZE (vocabolario controllato)
# ----------------------------------------------------------------------------
VARIABLES = {
    "rain":   {"label_it": "Pioggia", "label_en": "Rainfall", "unit": "mm",
               "aggregation": "sum", "plot": "bar",
               "notes": "Altezza di pioggia caduta nell'intervallo che termina al timestamp."},
    "vwc":    {"label_it": "Contenuto d'acqua volumetrico", "label_en": "Volumetric water content",
               "unit": "m3/m3", "aggregation": "mean", "plot": "line"},
    "psi":    {"label_it": "Potenziale matriciale", "label_en": "Matric potential",
               "unit": "kPa", "aggregation": "mean", "plot": "line",
               "notes": "Valori negativi = suzione. TEROS 21: accuratezza dichiarata tra -9 e -100 kPa; "
                        "valori > -9 kPa (prossimi alla saturazione) sono indicativi."},
    "t_soil": {"label_it": "Temperatura del terreno", "label_en": "Soil temperature",
               "unit": "°C", "aggregation": "mean", "plot": "line"},
    "t_air":  {"label_it": "Temperatura dell'aria (2 m)", "label_en": "Air temperature (2 m)",
               "unit": "°C", "aggregation": "mean", "plot": "line"},
}

# ----------------------------------------------------------------------------
# 4. REGOLE: (strumento, intestazione colonna) -> serie
#    series_id = <STAZIONE>_<GRANDEZZA>[_<PROFONDITÀ cm, 3 cifre>]
# ----------------------------------------------------------------------------
def classify(instrument_id, col):
    station = INSTRUMENTS[instrument_id]["station_id"]
    none = {"depth_m": None, "sensor_model": None, "logger_port": None}
    if instrument_id.endswith("_METER"):
        m = re.match(r"Port(\d)-(T10|T21)-(\d+)cm(?:\s*\((kPa|°C)\))?", col)
        if not m:
            return None
        port, sensor, depth, unit = int(m[1]), m[2], int(m[3]), m[4]
        if sensor == "T10":
            var, model = "vwc", "TEROS 10"
        else:
            var, model = ("psi" if unit == "kPa" else "t_soil"), "TEROS 21"
        code = {"vwc": "VWC", "psi": "PSI", "t_soil": "TSOIL"}[var]
        return f"{station}_{code}_{depth:03d}", {"variable": var, "depth_m": depth / 100,
                                                "sensor_model": model, "logger_port": port}
    if instrument_id.endswith("_RAIN") and col.startswith("Rain"):
        return f"{station}_RAIN", {"variable": "rain", **none}
    return None

# Note di qualità (interpretazione) emerse dall'analisi dei dati — orari in UTC+1
QC_NOTES = {
    "CF_21521_RAIN": ["Lacune: 2025-02-12 ~08:10–09:40 (1,5 h) e 2025-11-19 ~10:50–23:40 (12,8 h).",
                      "Alcuni timestamp con millisecondi spuri (es. 23:49:59.999): arrotondati al minuto."],
    "CF_18925_RAIN": ["Alcuni timestamp con millisecondi spuri (es. 00:10:00.001): arrotondati al minuto."],
    "CF_18957_RAIN": ["Alcuni timestamp con millisecondi spuri (es. 00:10:00.001): arrotondati al minuto."],
}
for _code in CF_STATIONS:
    QC_NOTES[f"CF_{_code}_RAIN"].insert(0, "Anno 2025 completo (passo 10 min). Totale annuo 2025: "
        + {"21521": "1167 mm (125 m)", "18925": "1227 mm (356 m)", "18957": "987 mm (28 m)"}[_code]
        + "; correlazione giornaliera Cologna–Pellezzano 0.95, Cologna–Salerno G.C. 0.83, Pellezzano–Salerno G.C. 0.75.")
_METER_NOTE_PSI = ("Registrazioni dal 2025-01-31 (installazione). In estate il terreno superficiale raggiunge valori "
                   "molto negativi (fino a circa -3000 kPa): oltre i -100 kPa il TEROS 21 è fuori dall'intervallo "
                   "di accuratezza dichiarato, i valori sono indicativi. Le rapide risalite (es. 2025-09-11, "
                   "2025-10-23) sono fronti di bagnamento reali dopo la pioggia, non spike.")
for _k in METER_STATIONS:
    _d = {"NV": ("015", "030", "060"), "V": ("015", "035", "060")}[_k]
    for _x in _d:
        QC_NOTES[f"SAL_{_k}_PSI_{_x}"] = [_METER_NOTE_PSI]
        QC_NOTES[f"SAL_{_k}_VWC_{_x}"] = ["Registrazioni dal 2025-01-31 (installazione)."]
        QC_NOTES[f"SAL_{_k}_TSOIL_{_x}"] = ["Registrazioni dal 2025-01-31 (installazione)."]
QC_NOTES["SAL_NV_PSI_015"].append("Gap 2025-12-18 ~02:40 (10,7 h) condiviso da tutte le serie della stazione.")
QC_NOTES["SAL_V_PSI_015"].append("Gap 2025-12-17 16:20 – 2025-12-18 17:25 (25 h) condiviso da tutte le serie della stazione.")


# ----------------------------------------------------------------------------
# 5. SCANSIONE FILE
# ----------------------------------------------------------------------------
def match_instrument(fname):
    for iid, inst in INSTRUMENTS.items():
        if inst.get("file_pattern") and re.search(inst["file_pattern"], fname, re.I):
            return iid
    return None

series, files, gaps = {}, [], {}
for path in sorted(glob.glob(os.path.join(INPUT_DIR, "*.xlsx"))):
    fname = os.path.basename(path)
    clean = re.sub(r"^[0-9a-f]{8}-", "", fname)          # rimuove prefisso di upload
    iid = match_instrument(clean)
    if iid is None:
        print("!! file non riconosciuto:", fname); continue
    year = int(re.search(r"(\d{4})\.xlsx$", clean)[1])
    for sheet, df in pd.read_excel(path, sheet_name=None, engine="openpyxl").items():
        # Colonna B = 'Time UTC+1' -> riferimento temporale unico (arrotondato al minuto)
        t = pd.to_datetime(df.iloc[:, 1], errors="coerce").dt.round("min")
        has_local = bool(df.iloc[:, 0].notna().any())
        for col in df.columns[2:]:
            if not isinstance(col, str) or col.startswith("Unnamed"):
                continue                                  # colonne di servizio (calcolo offset)
            cl = classify(iid, col)
            if cl is None:
                continue
            sid, attrs = cl
            v = pd.to_numeric(df[col], errors="coerce")
            ok = v.notna() & t.notna()
            tt, vv = t[ok].reset_index(drop=True), v[ok]
            step = tt.diff().dt.total_seconds().div(60)
            for i in step.index[step > GAP_REPORT_H * 60]:
                gaps.setdefault(sid, []).append({"from": iso(tt[i - 1]), "to": iso(tt[i]),
                                                 "hours": round(step[i] / 60, 1)})
            files.append({"series_id": sid, "file": clean, "year": year, "sheet": sheet,
                          "column": col, "time_column": "Time UTC+1",
                          "local_time_column": "Time" if has_local else None,
                          "start": iso(tt.min()), "end": iso(tt.max()),
                          "n_values": int(ok.sum()), "n_rows": int(len(df)),
                          "step_min": float(step.median()),
                          "gaps_over_1h": int((step > 60).sum()),
                          "max_gap_h": round(float(step.max()) / 60, 1),
                          "min": float(vv.min()), "max": float(vv.max()),
                          "mean": round(float(vv.mean()), 4)})
            series.setdefault(sid, {"series_id": sid, "instrument_id": iid,
                                    "station_id": INSTRUMENTS[iid]["station_id"], **attrs})

# Modalità incrementale: i file Excel già elaborati possono non essere presenti in INPUT_DIR.
# Si parte dal catalogo esistente e si aggiornano solo i file/le serie scansionati (REBUILD=1: da zero).
# Per aggiungere un anno a una serie esistente servono comunque tutti i suoi file Excel (li rielabora standardize.py).
prev = {}
_prev_path = os.path.join(OUT_DIR, "catalog.json")
if os.path.exists(_prev_path) and os.environ.get("REBUILD") != "1":
    prev = json.load(open(_prev_path, encoding="utf-8"))
_scanned = {(f["file"], f["sheet"], f["column"]) for f in files}
files += [f for f in prev.get("source_files", []) if (f["file"], f["sheet"], f["column"]) not in _scanned]

# ----------------------------------------------------------------------------
# 6. ASSEMBLAGGIO
# ----------------------------------------------------------------------------
def label(s):
    v = VARIABLES[s["variable"]]
    d = f" {int(s['depth_m']*100)} cm" if s.get("depth_m") is not None else ""
    return f"{v['label_it']}{d} – {STATIONS[s['station_id']]['name']}"

out_series = []
for sid, s in sorted(series.items()):
    fs = sorted((f for f in files if f["series_id"] == sid), key=lambda f: f["start"])
    v = VARIABLES[s["variable"]]
    out_series.append({
        **s,
        "label": label(s),
        "unit": v["unit"],
        "aggregation": v["aggregation"],
        "step_min": min(f["step_min"] for f in fs),
        "start": fs[0]["start"],
        "end": max(f["end"] for f in fs),
        "n_values": sum(f["n_values"] for f in fs),
        "min": min(f["min"] for f in fs),
        "max": max(f["max"] for f in fs),
        "status": "active",
        "gaps_over_24h": gaps.get(sid, []),
        "qc_notes": QC_NOTES.get(sid, []),
        "sources": [{k: f[k] for k in ("file", "year", "sheet", "column")} for f in fs],
    })

_new_ids = {s["series_id"] for s in out_series}
out_series += [s for s in prev.get("series", []) if s["series_id"] not in _new_ids]   # serie non riscansionate
out_series.sort(key=lambda s: s["series_id"])

def dist_km(a, b):
    p1, p2 = math.radians(a["lat"]), math.radians(b["lat"])
    dp, dl = p2 - p1, math.radians(b["lon"] - a["lon"])
    h = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return round(2 * 6371 * math.asin(math.sqrt(h)), 2)

stations_out = []
for k, st in STATIONS.items():
    # distanze solo verso le stazioni in situ (le stazioni "griglia" non si elencano tra loro)
    near = {o: dist_km(st, ost) for o, ost in STATIONS.items()
            if o != k and st["lat"] is not None and ost["lat"] is not None and ost.get("kind") != "grid"}
    stations_out.append({"station_id": k, **st, "distance_km": near or None})

catalog = {
    "catalog_version": prev.get("catalog_version", "0.2"),
    "generated": datetime.now(TZ).isoformat(timespec="seconds"),
    "conventions": {**prev.get("conventions", {}),
        "time_reference": "UTC+1 (ora solare italiana / CET, senza ora legale) per tutte le misure. "
                          "Timestamp in ISO 8601 con offset '+01:00'. Sorgente: colonna 'Time UTC+1' "
                          "dei file Excel, arrotondata al minuto. La colonna 'Time' (ora locale con "
                          "ora legale) non è usata.",
        "series_id": "<STAZIONE>_<GRANDEZZA>[_<PROFONDITÀ cm, 3 cifre>]",
        "depth_m": "Profondità sotto il piano campagna, positiva verso il basso, in metri.",
        "coordinates": "WGS84 (EPSG:4326), gradi decimali.",
    },
    "variables": {**prev.get("variables", {}), **VARIABLES},
    "stations": stations_out,
    "instruments": [{"instrument_id": k, **{kk: vv for kk, vv in v.items() if kk != "file_pattern"}}
                    for k, v in INSTRUMENTS.items()],
    "series": out_series,
    "source_files": files,
}

if "qc_rules" in prev:
    catalog["qc_rules"] = prev["qc_rules"]
os.makedirs(OUT_DIR, exist_ok=True)
with open(os.path.join(OUT_DIR, "catalog.json"), "w", encoding="utf-8") as fh:
    json.dump(catalog, fh, ensure_ascii=False, indent=2)

st_lookup = {s["station_id"]: s for s in stations_out}
flat = pd.DataFrame([{k: s[k] for k in ("series_id", "station_id", "instrument_id", "variable",
                       "label", "unit", "depth_m", "sensor_model", "logger_port", "aggregation",
                       "step_min", "start", "end", "n_values", "min", "max", "status")}
                     | {"lat": st_lookup[s["station_id"]]["lat"],
                        "lon": st_lookup[s["station_id"]]["lon"],
                        "n_gaps_over_24h": len(s["gaps_over_24h"]),
                        "qc_notes": " | ".join(s["qc_notes"]),
                        "source_files": "; ".join(x["file"] for x in s["sources"])}
                     for s in out_series])
flat.to_csv(os.path.join(OUT_DIR, "series.csv"), index=False, encoding="utf-8-sig")
print(flat[["series_id", "variable", "unit", "step_min", "start", "end", "n_values", "n_gaps_over_24h"]].to_string())
