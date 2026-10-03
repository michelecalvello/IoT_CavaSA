"""
ingest_gee.py — passo 2b: prodotti a griglia scaricati da Google Earth Engine (ERA5-Land, SMAP L4).

Legge   : OUT_DIR/catalog.json, qc_rules.json, CSV GEE in INPUT_DIR
            cavasa_ERA5Land[_<anno>].csv  time,lon,lat,volumetric_soil_water_layer_1,volumetric_soil_water_layer_2,
                                        total_precipitation_hourly,temperature_2m
            cavasa_SMAP_L4[_<anno>].csv   time,lon,lat,sm_surface,sm_rootzone,sm_profile
Scrive  : data/raw/<serie>/<anno>.csv, data/agg/<serie>_{1h,1d}.csv, data/qc_log.csv,
          catalog.json e series.csv (stazioni "griglia", strumenti, serie)

Ogni pixel della griglia è una STAZIONE (kind = "grid"):
  ERA5-Land  ERA5_E<lon*10>_N<lat*10>        pixel 0,1° (coordinate agganciate a 0,1°: GEE restituisce
                                              centri con scarti di ~0,0009°)
  SMAP L4    SMAP_E<lon*100>_N<lat*100>      pixel ~9 km (coordinate = centro pixel)
Serie: <STAZIONE>_<GRANDEZZA>_<STRATO>; per i prodotti a strati il codice a 3 cifre è la profondità
INFERIORE dello strato in cm (ERA5: 007 = 0–7 cm, 028 = 7–28 cm; SMAP: 005 = 0–5 cm, 100 = 0–100 cm,
PRF = intero profilo fino al substrato); la profondità esatta è in `layer_label` / `depth_range_m`.

CONVERSIONI
- tempo: GEE fornisce UTC -> +1 h = UTC+1 (ora solare). Il dato UTC delle 23:00 del 31/12 diventa 00:00 del 1/1
  e finisce nel file dell'anno successivo.
- ERA5 total_precipitation_hourly: m -> mm (x1000). Pioggia nell'ora che TERMINA al timestamp GEE (convenzione
  ECMWF, ipotesi da verificare con i pluviometri); valori negativi di arrotondamento numerico azzerati.
- ERA5 temperature_2m: K -> °C.
- VWC (ERA5 layer 1/2, SMAP surface/rootzone/profile): m3/m3, nessuna conversione.

MODALITÀ INCREMENTALE: i CSV annuali già elaborati non servono più. Per ogni serie si uniscono i dati nuovi
a quelli già presenti in data/raw (a parità di istante vince il dato nuovo), si rifanno QC e aggregati.
Aggiungere un periodo = passare il solo nuovo CSV (anche multi-anno, es. cavasa_ERA5Land.csv).
REPLACE_PRODUCTS=ERA5,SMAP sostituisce integralmente i dati del prodotto (nuova esportazione con altri punti). Dopo l'ingestione rieseguire derive_cumulative.py.
"""
import glob, json, math, os, re
import numpy as np
import pandas as pd

INPUT_DIR = os.environ.get("INPUT_DIR", ".")
OUT_DIR = os.environ.get("OUT_DIR", ".")
CATALOG = os.path.join(OUT_DIR, "catalog.json")
RULES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "qc_rules.json")
TFMT = "%Y-%m-%d %H:%M"
ISO = "%Y-%m-%dT%H:%M+01:00"
UTC_TO_LOCAL = pd.Timedelta(hours=1)          # UTC -> UTC+1

# ----------------------------------------------------------------------------
# prodotti e colonne
#   colonna GEE: (codice grandezza, variabile, codice strato, (top, bottom) m, etichetta strato, conversione)
# ----------------------------------------------------------------------------
PRODUCTS = {
    "ERA5": {
        "pattern": r"^cavasa_ERA5Land(?:_\d{4})?\.csv$",
        "name": "ERA5-Land",
        "owner": "ECMWF – Copernicus Climate Change Service (ERA5-Land)",
        "manufacturer": "ECMWF / Copernicus C3S",
        "type": "Rianalisi ERA5-Land oraria (collezione GEE ECMWF/ERA5_LAND/HOURLY), pixel 0,1° (~9 km)",
        "sensor_model": "ERA5-Land",
        "step_min": 60.0,
        "station_notes": "Pixel ERA5-Land 0,1° (~9 km): media sull'area del pixel, non misura puntuale. "
                         "Dati scaricati da Google Earth Engine (collezione ECMWF/ERA5_LAND/HOURLY).",
        "columns": {
            "volumetric_soil_water_layer_1": ("VWC", "vwc", "007", (0.0, 0.07), "0–7 cm", lambda v: v),
            "volumetric_soil_water_layer_2": ("VWC", "vwc", "028", (0.07, 0.28), "7–28 cm", lambda v: v),
            "total_precipitation_hourly": ("RAIN", "rain", None, None, None, lambda v: (v * 1000).clip(lower=0)),
            "temperature_2m": ("TAIR", "t_air", None, None, None, lambda v: v - 273.15),
        },
        "notes": {
            "vwc": "Contenuto d'acqua stimato dal modello ERA5-Land (strato {layer}) come media del pixel di 0,1° (~9 km): "
                   "confronto con i sensori in situ solo indicativo (scala e profondità diverse). Valore istantaneo orario.",
            "rain": "Pioggia oraria ERA5-Land (total_precipitation_hourly, convertita da m a mm). Ipotesi: il valore si riferisce "
                    "all'ora che termina al timestamp GEE (convenzione ECMWF), coerente con i pluviometri; da verificare sul "
                    "confronto con i pluviometri del Centro Funzionale. Valori negativi di arrotondamento numerico azzerati.",
            "t_air": "Temperatura dell'aria a 2 m ERA5-Land (convertita da K a °C), valore istantaneo orario del pixel.",
        },
    },
    "SMAP": {
        "pattern": r"^cavasa_SMAP_L4(?:_\d{4})?\.csv$",
        "name": "SMAP L4",
        "owner": "NASA – SMAP Level 4 (GMAO)",
        "manufacturer": "NASA GMAO",
        "type": "Prodotto SMAP L4 Soil Moisture (collezione GEE NASA/SMAP/SPL4SMGP), pixel ~9 km, medie su 3 h",
        "sensor_model": "SMAP L4",
        "step_min": 180.0,
        "station_notes": "Pixel SMAP L4 (griglia EASE ~9 km): media sull'area del pixel, non misura puntuale. "
                         "Medie su 3 h con etichetta a metà intervallo (01:30 UTC = 02:30 UTC+1). "
                         "Dati scaricati da Google Earth Engine (collezione NASA/SMAP/SPL4SMGP).",
        "columns": {
            "sm_surface": ("VWC", "vwc", "005", (0.0, 0.05), "0–5 cm", lambda v: v),
            "sm_rootzone": ("VWC", "vwc", "100", (0.0, 1.0), "0–100 cm", lambda v: v),
            "sm_profile": ("VWC", "vwc", "PRF", None, "profilo intero", lambda v: v),
        },
        "notes": {
            "vwc": "Contenuto d'acqua SMAP L4 (strato {layer}) come media del pixel (~9 km) su intervalli di 3 h; "
                   "l'etichetta temporale è il punto medio dell'intervallo. Confronto con i sensori in situ solo indicativo. "
                   "Strato 'profilo intero': dalla superficie al substrato del modello (profondità variabile)."
        },
    },
}

NEW_VARIABLE = {"t_air": {"label_it": "Temperatura dell'aria (2 m)", "label_en": "Air temperature (2 m)",
                          "unit": "°C", "aggregation": "mean", "plot": "line"}}
DECIMALS = {"vwc": 4, "rain": 3, "t_air": 2}


def dist_km(a, b):
    p1, p2 = math.radians(a["lat"]), math.radians(b["lat"])
    dp, dl = p2 - p1, math.radians(b["lon"] - a["lon"])
    h = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return round(2 * 6371 * math.asin(math.sqrt(h)), 2)


catalog = json.load(open(CATALOG, encoding="utf-8"))
# REPLACE_PRODUCTS=ERA5,SMAP: elimina prima tutte le stazioni/serie/dati (anche derivati) di quei prodotti, poi reingerisce
for _p in [x for x in os.environ.get("REPLACE_PRODUCTS", "").split(",") if x]:
    _ids = {s["station_id"] for s in catalog["stations"] if s["station_id"].startswith(_p + "_")}
    for _s in [s for s in catalog["series"] if s["station_id"] in _ids]:
        for _d in glob.glob(os.path.join(OUT_DIR, "data", "raw", _s["series_id"])):
            for _f in glob.glob(os.path.join(_d, "*")):
                os.remove(_f)
            os.rmdir(_d)
        for _f in glob.glob(os.path.join(OUT_DIR, "data", "agg", _s["series_id"] + "_*.csv")):
            os.remove(_f)
    catalog["series"] = [s for s in catalog["series"] if s["station_id"] not in _ids]
    catalog["stations"] = [s for s in catalog["stations"] if s["station_id"] not in _ids]
    catalog["instruments"] = [i for i in catalog["instruments"] if i["station_id"] not in _ids]
    _qc = os.path.join(OUT_DIR, "data", "qc_log.csv")
    if os.path.exists(_qc):
        _l = pd.read_csv(_qc); _l[~_l["series_id"].str.startswith(_p + "_")].to_csv(_qc, index=False)
    _fl = os.path.join(OUT_DIR, "series.csv")
    if os.path.exists(_fl):
        _f = pd.read_csv(_fl, encoding="utf-8-sig"); _f[~_f["station_id"].isin(_ids)].to_csv(_fl, index=False, encoding="utf-8-sig")
    print(f"{_p}: rimossi {len(_ids)} pixel e relative serie")
rules = json.load(open(RULES, encoding="utf-8"))
catalog["variables"] = {**catalog["variables"], **{k: v for k, v in NEW_VARIABLE.items() if k not in catalog["variables"]}}
variables = catalog["variables"]

# ----------------------------------------------------------------------------
# lettura CSV GEE -> {(station_id, serie): DataFrame(time, value)}
# ----------------------------------------------------------------------------
def read_gee(path, key):
    d = pd.read_csv(path)
    d["time"] = (pd.to_datetime(d["time"]) + UTC_TO_LOCAL).dt.round("min")
    if key == "ERA5":
        d["glon"], d["glat"] = d["lon"].round(1), d["lat"].round(1)       # centri pixel 0,1° (scarti GEE ~0,0009° eliminati)
        d["members"] = None
    else:
        # I punti campionati da GEE sono più fitti dei pixel nativi SMAP (~9 km): punti diversi possono cadere nello
        # stesso pixel e hanno allora serie identiche. Si tiene un solo pixel per gruppo di punti identici,
        # posizionato nel baricentro dei punti (posizione indicativa: il centro nativo del pixel non è noto).
        d["glon"], d["glat"] = d["lon"].round(4), d["lat"].round(4)
        val = ["sm_surface", "sm_rootzone", "sm_profile"]
        sig = {}
        for pt, g in d.groupby(["glon", "glat"]):
            sig.setdefault(hash(g.sort_values("time")[val].round(6).to_numpy().tobytes()), []).append(pt)
        remap = {}
        for pts in sig.values():
            lon, lat = round(sum(p[0] for p in pts) / len(pts), 4), round(sum(p[1] for p in pts) / len(pts), 4)
            for p in pts:
                remap[p] = (lon, lat, pts)
        keep = {pts[0] for _, _, pts in remap.values()}
        d = d[[pt in keep for pt in zip(d["glon"], d["glat"])]].copy()
        m = [remap[pt] for pt in zip(d["glon"], d["glat"])]
        d["glon"], d["glat"], d["members"] = [x[0] for x in m], [x[1] for x in m], [x[2] for x in m]
    dec = 1 if key == "ERA5" else 2
    f = 10 if key == "ERA5" else 100
    d["station_id"] = (key + "_E" + (d["glon"] * f).round().astype(int).astype(str)
                       + "_N" + (d["glat"] * f).round().astype(int).astype(str))
    return d


new_data, new_stations, new_files = {}, {}, {}
for key, P in PRODUCTS.items():
    for path in sorted(glob.glob(os.path.join(INPUT_DIR, "*.csv"))):
        fname = re.sub(r"^[0-9a-f]{8}-", "", os.path.basename(path))       # rimuove prefisso di upload
        m = re.match(P["pattern"], fname)
        if not m:
            continue
        d = read_gee(path, key)
        print(f"{fname}: {len(d)} righe, {d['station_id'].nunique()} pixel, {d['time'].min()} → {d['time'].max()} (UTC+1)")
        for st, g in d.groupby("station_id"):
            new_stations[st] = {"product": key, "lat": round(float(g["glat"].iloc[0]), 4 if key == "SMAP" else 1),
                                "lon": round(float(g["glon"].iloc[0]), 4 if key == "SMAP" else 1),
                                "members": g["members"].iloc[0]}
            for col, (gcode, var, lcode, rng, llabel, conv) in P["columns"].items():
                sid = f"{st}_{gcode}" + (f"_{lcode}" if lcode else "")
                v = conv(pd.to_numeric(g[col], errors="coerce"))
                part = pd.DataFrame({"time": g["time"], "value": v}).dropna()
                new_data.setdefault(sid, []).append(part)
                for year in sorted(part["time"].dt.year.unique().tolist()):       # il file può coprire più anni
                    new_files.setdefault(sid, {})[year] = {"file": fname, "year": int(year), "sheet": None, "column": col}
                if col == "total_precipitation_hourly":
                    n_neg = int((pd.to_numeric(g[col], errors="coerce") < 0).sum())
                    if n_neg:
                        new_data.setdefault(("neg", sid), []).append(n_neg)
if not new_data:
    raise SystemExit("Nessun file GEE riconosciuto in INPUT_DIR (cavasa_ERA5Land_<anno>.csv, cavasa_SMAP_L4_<anno>.csv)")

# ----------------------------------------------------------------------------
# stazioni e strumenti "griglia"
# ----------------------------------------------------------------------------
stations = {s["station_id"]: s for s in catalog["stations"]}
instruments = {i["instrument_id"]: i for i in catalog["instruments"]}
insitu = [s for s in stations.values() if s.get("kind") != "grid"]
for st, info in new_stations.items():
    P = PRODUCTS[info["product"]]
    dec = 2 if info["product"] == "SMAP" else 1
    prev = stations.get(st, {})
    stations[st] = {
        "station_id": st, "name": f"{P['name']} · {info['lon']:.{dec}f}°E {info['lat']:.{dec}f}°N",
        "site": "Area di Salerno", "municipality": None,
        "lat": info["lat"], "lon": info["lon"], "elevation_m": None,
        "owner": P["owner"], "external_code": None,
        "kind": "grid", "product": P["name"],
        "notes": P["station_notes"] + (
            f" Il pixel comprende {len(info['members'])} punti di campionamento GEE con valori identici "
            f"(lon, lat: {'; '.join(f'{a:.4f}, {b:.4f}' for a, b in sorted(info['members']))}); "
            "posizione indicativa (baricentro dei punti)." if info["members"] and len(info["members"]) > 1 else ""),
        "sample_points": [list(p) for p in sorted(info["members"])] if info["members"] else None,
        "distance_km": {o["station_id"]: dist_km({"lat": info["lat"], "lon": info["lon"]}, o)
                        for o in insitu if o.get("lat") is not None} or None,
    }
    instruments[f"{st}_MODEL"] = {"instrument_id": f"{st}_MODEL", "station_id": st,
                                  "manufacturer": P["manufacturer"], "type": P["type"]}
catalog["stations"] = list(stations.values())
catalog["instruments"] = list(instruments.values())

# ----------------------------------------------------------------------------
# QC, scrittura, aggregati
# ----------------------------------------------------------------------------
def aggregate(d, agg, step, rule, var):
    expected = max(1.0, {"1h": 60 / step, "1d": 1440 / step}[rule])
    freq = {"1h": "1h", "1d": "1D"}[rule]
    g = d[d["flag"] < 2].set_index("time")["value"]
    if agg == "sum":
        rs = g.resample(freq, closed="right", label="left")
        out = pd.DataFrame({"value": rs.sum(min_count=1), "max": rs.max(), "n": rs.count()})
    else:
        rs = g.resample(freq)
        out = pd.DataFrame({"value": rs.mean(), "min": rs.min(), "max": rs.max(), "n": rs.count()})
    out = out[out["n"] > 0].copy()
    out["coverage"] = (out["n"] / expected).clip(upper=1).round(3)
    out["n"] = out["n"].astype(int)
    for c in ("value", "min", "max"):
        if c in out:
            out[c] = out[c].round(DECIMALS[var])
    out.index = out.index.strftime(TFMT)
    out.index.name = "time"
    return out


def col_info(sid):
    """(prodotto, colonna GEE, info colonna) di una serie."""
    for key, P in PRODUCTS.items():
        if sid.startswith(key + "_"):
            for col, info in P["columns"].items():
                gcode, var, lcode, rng, llabel, conv = info
                if sid.endswith(f"_{gcode}" + (f"_{lcode}" if lcode else "")):
                    return key, col, info


log, processed = [], []
old_series = {s["series_id"]: s for s in catalog["series"]}
out_series = []
for sid in sorted(k for k in new_data if isinstance(k, str)):
    key, col, (gcode, var, lcode, rng, llabel, conv) = col_info(sid)
    P = PRODUCTS[key]
    station_id = sid[: sid.index(f"_{gcode}")]
    st = stations[station_id]
    raw_dir = os.path.join(OUT_DIR, "data", "raw", sid)
    os.makedirs(raw_dir, exist_ok=True)

    # dati già presenti + dati nuovi (vince il nuovo)
    old = [pd.read_csv(p, usecols=["time", "value"]) for p in sorted(glob.glob(os.path.join(raw_dir, "*.csv")))]
    old = [o.assign(time=pd.to_datetime(o["time"])) for o in old]
    d = pd.concat(old + new_data[sid], ignore_index=True)
    d = d.drop_duplicates("time", keep="last").sort_values("time").reset_index(drop=True)

    flag = np.zeros(len(d), dtype=np.int8)
    r = rules["range_by_variable"].get(var)
    if r:
        bad = ((d["value"] < r["min"]) | (d["value"] > r["max"])).to_numpy()
        flag[bad] = 2
        if bad.any():
            log.append({"series_id": sid, "rule": f"fuori range [{r['min']}, {r['max']}]", "flag": 2, "n": int(bad.sum()),
                        "first": d["time"][bad].min().strftime(TFMT), "last": d["time"][bad].max().strftime(TFMT)})
    d["flag"] = flag
    n_neg = sum(new_data.get(("neg", sid), []))
    if n_neg:
        log.append({"series_id": sid, "rule": "valori negativi di arrotondamento numerico azzerati", "flag": 0,
                    "n": n_neg, "first": None, "last": None})

    years = sorted(d["time"].dt.year.unique().tolist())
    for y in years:
        dy = d[d["time"].dt.year == y].copy()
        dy["time"] = dy["time"].dt.strftime(TFMT)
        dy["value"] = dy["value"].round(DECIMALS[var])
        dy.to_csv(os.path.join(raw_dir, f"{y}.csv"), index=False)
    for rule in ("1h", "1d"):
        aggregate(d, variables[var]["aggregation"], P["step_min"], rule, var) \
            .to_csv(os.path.join(OUT_DIR, "data", "agg", f"{sid}_{rule}.csv"))

    dt = d["time"].diff()
    gaps = [{"from": d["time"][i - 1].strftime(ISO), "to": d["time"][i].strftime(ISO),
             "hours": round(dt[i].total_seconds() / 3600, 1)} for i in dt.index[dt > pd.Timedelta(hours=24)]]
    good = d[d["flag"] < 2]
    counts = d["flag"].value_counts().to_dict()
    notes = [P["notes"][var].format(layer=llabel)]
    if var == "rain":
        tot = good.groupby((good["time"] - pd.Timedelta(minutes=1)).dt.year)["value"].sum().round(0)   # il dato delle 00:00 del 1/1 chiude l'anno prima
        notes.append("Totali annui (dati disponibili): " + "; ".join(f"{y}: {v:.0f} mm" for y, v in tot.items()) + ".")
    prev = old_series.get(sid, {})
    srcs = {x["year"]: x for x in prev.get("sources", [])}
    srcs.update(new_files[sid])
    label = f"{variables[var]['label_it']}{' ' + llabel if llabel else ''} – {st['name']}"
    out_series.append({
        "series_id": sid, "instrument_id": f"{station_id}_MODEL", "station_id": station_id, "variable": var,
        "depth_m": rng[1] if rng else None, "depth_range_m": list(rng) if rng else None, "layer_label": llabel,
        "sensor_model": P["sensor_model"], "logger_port": None,
        "label": label, "unit": variables[var]["unit"], "aggregation": variables[var]["aggregation"],
        "step_min": P["step_min"],
        "start": d["time"].iloc[0].strftime(ISO), "end": d["time"].iloc[-1].strftime(ISO),
        "n_values": int(len(d)), "min": round(float(good["value"].min()), DECIMALS[var]),
        "max": round(float(good["value"].max()), DECIMALS[var]),
        "status": "active", "gaps_over_24h": gaps, "qc_notes": notes,
        "sources": [srcs[y] for y in sorted(srcs)],
        "data": {"format": "csv: time,value,flag", "raw": f"data/raw/{sid}/{{year}}.csv", "years": years,
                 "agg_1h": f"data/agg/{sid}_1h.csv", "agg_1d": f"data/agg/{sid}_1d.csv"},
        "qc_summary": {"n_total": int(len(d)), "n_valid": int(counts.get(0, 0)), "n_suspect": int(counts.get(1, 0)),
                       "n_bad": int(counts.get(2, 0)), "duplicates_removed": 0},
    })
    processed.append(sid)
    print(f"{sid:28s} n={len(d):6d} {d['time'].iloc[0]} → {d['time'].iloc[-1]}  min={out_series[-1]['min']}  max={out_series[-1]['max']}")

# ----------------------------------------------------------------------------
# catalogo, qc_log, series.csv
# ----------------------------------------------------------------------------
catalog["series"] = sorted([s for s in catalog["series"] if s["series_id"] not in set(processed)] + out_series,
                           key=lambda s: s["series_id"])
catalog["generated"] = pd.Timestamp.now(tz="Etc/GMT-1").strftime("%Y-%m-%dT%H:%M:%S+01:00")
catalog["conventions"]["gridded_products"] = (
    "Stazioni con kind = 'grid' sono pixel di prodotti a griglia scaricati da Google Earth Engine (pipeline/ingest_gee.py): "
    "ERA5-Land (0,1°, orario) e SMAP L4 (~9 km, medie su 3 h). I valori sono medie di pixel, non misure puntuali. "
    "I timestamp GEE (UTC) sono convertiti in UTC+1 (+1 h). Per i prodotti a strati il codice a 3 cifre della serie è la "
    "profondità inferiore dello strato in cm (PRF = intero profilo); vedi layer_label e depth_range_m.")
catalog["qc_rules"] = {k: v for k, v in rules.items() if not k.startswith("_")}
with open(CATALOG, "w", encoding="utf-8") as fh:
    json.dump(catalog, fh, ensure_ascii=False, indent=2)

qc_path = os.path.join(OUT_DIR, "data", "qc_log.csv")
_log = pd.DataFrame(log, columns=["series_id", "rule", "flag", "n", "first", "last"])
if os.path.exists(qc_path):
    _old = pd.read_csv(qc_path)
    _log = pd.concat([_old[~_old["series_id"].isin(processed)], _log], ignore_index=True)
    _log = _log.sort_values("series_id", kind="stable").reset_index(drop=True)
_log.to_csv(qc_path, index=False)

flat_path = os.path.join(OUT_DIR, "series.csv")
if os.path.exists(flat_path):
    flat = pd.read_csv(flat_path, encoding="utf-8-sig")
    flat = flat[~flat["series_id"].isin(processed)]
    rows = [{**{k: x[k] for k in ("series_id", "station_id", "instrument_id", "variable", "label", "unit",
                                  "depth_m", "sensor_model", "logger_port", "aggregation", "step_min",
                                  "start", "end", "n_values", "min", "max", "status")},
             "lat": stations[x["station_id"]]["lat"], "lon": stations[x["station_id"]]["lon"],
             "n_gaps_over_24h": len(x["gaps_over_24h"]), "qc_notes": " | ".join(x["qc_notes"]),
             "source_files": "; ".join(s["file"] for s in x["sources"])} for x in out_series]
    pd.concat([flat, pd.DataFrame(rows)], ignore_index=True)[flat.columns] \
        .sort_values("series_id").to_csv(flat_path, index=False, encoding="utf-8-sig")
print(f"\n{len(processed)} serie, {len(new_stations)} stazioni griglia aggiornate")
