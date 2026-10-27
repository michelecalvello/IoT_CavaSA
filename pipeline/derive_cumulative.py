"""
derive_cumulative.py — passo 3: serie derivate di pioggia cumulata (dal 1 gennaio e su 24 h).

Legge   : OUT_DIR/catalog.json (dopo standardize.py), qc_rules.json (blocco "cumulative"),
          OUT_DIR/data/raw/<pluviometro>/<anno>.csv
Scrive  : OUT_DIR/data/raw/<STAZIONE>_RAINCUM/<anno>.csv     time,value,flag
          OUT_DIR/data/agg/<STAZIONE>_RAINCUM_1h.csv / _1d.csv   time,value,n,coverage,flag
          OUT_DIR/data/raw/<STAZIONE>_RAINCUM24/<anno>.csv   time,value,flag
          OUT_DIR/data/agg/<STAZIONE>_RAINCUM24_1h.csv / _1d.csv   time,value,n,coverage,flag
Aggiorna: OUT_DIR/catalog.json e OUT_DIR/series.csv (le serie derivate hanno "derived_from").
Lo script è idempotente: a ogni esecuzione rigenera tutte le serie derivate.

PIOGGIA CUMULATA 24 h (RAINCUM24, variabile rain_cum24)
- value = pioggia (mm) caduta nelle 24 h che terminano al timestamp, cioè nella finestra (t − 24 h, t]
  (finestra mobile, non si azzera mai). Una riga per ogni dato di pioggia valido, sui soli dati con flag < 2.
  La finestra scorre sull'intera serie, anche a cavallo di fine anno.
- flag = 1 quando la cumulata 24 h è probabilmente sottostimata per dati mancanti: nella finestra mancano
  più di window24_missing_hours_threshold ore di dati (qc_rules.json -> "cumulative"). Succede anche nelle prime
  24 h della serie, quando la finestra comincia prima del primo dato.
- Aggregati (etichetta = inizio intervallo): value = MASSIMO della cumulata 24 h nell'intervallo (per l'aggregato
  giornaliero è il picco delle 24 h mobili del giorno, non il totale del giorno, che è già in <STAZIONE>_RAIN_1d);
  n = numero di dati di pioggia nell'intervallo; coverage = n / attesi; flag = massimo flag.

DEFINIZIONE
- Ogni dato di pioggia è riferito all'intervallo che TERMINA al suo timestamp. L'anno di
  cumulo di un dato è quindi l'anno di (timestamp - 1 min): il dato delle 00:00 del 1 gennaio
  chiude l'anno precedente (è l'ultimo valore del suo totale annuo) e la cumulata riparte da 0
  con il dato successivo. Stessa convenzione dei cumuli orari/giornalieri del catalogo.
- value = pioggia cumulata (mm) dal 1 gennaio ore 00:00 fino al timestamp, sui soli dati
  con flag < 2. Una riga per ogni dato di pioggia valido.
- flag = 1 quando la cumulata è probabilmente sottostimata per dati mancanti (vedi
  qc_rules.json -> "cumulative"); 0 altrimenti.

AGGREGATI (etichetta = inizio intervallo, come per la pioggia)
- value = ultimo valore della cumulata nell'intervallo (cumulata a fine intervallo)
- n = numero di dati di pioggia nell'intervallo; coverage = n / attesi; flag = massimo flag
"""
import json, os
import numpy as np
import pandas as pd

OUT_DIR = os.environ.get("OUT_DIR", ".")
CATALOG = os.path.join(OUT_DIR, "catalog.json")
RULES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "qc_rules.json")
TFMT = "%Y-%m-%d %H:%M"
ISO = "%Y-%m-%dT%H:%M+01:00"
SUFFIX = "RAINCUM"
SUFFIX24 = "RAINCUM24"
WINDOW_H = 24

VARIABLE24 = {
    "label_it": "Pioggia cumulata 24 h",
    "label_en": "24-hour rainfall",
    "unit": "mm", "aggregation": "max", "plot": "line",
    "notes": "Pioggia caduta nelle 24 h che terminano al timestamp (finestra mobile). Tratteggiata dove i dati mancanti "
             "nella finestra la rendono sottostimata (anche nelle prime 24 h della serie). Negli aggregati orario e "
             "giornaliero è il massimo della cumulata 24 h nell'intervallo.",
}

VARIABLE = {
    "label_it": "Pioggia cumulata dal 1° gennaio",
    "label_en": "Year-to-date rainfall",
    "unit": "mm", "aggregation": "last", "plot": "line",
    "notes": "Pioggia cumulata dal 1 gennaio ore 00:00 dell'anno corrente (si azzera ogni anno). "
             "Tratteggiata dove i dati mancanti rendono la cumulata sottostimata.",
}

catalog = json.load(open(CATALOG, encoding="utf-8"))
rules = json.load(open(RULES, encoding="utf-8"))["cumulative"]
GAP_MIN, MISS_THR_H = rules["gap_min_minutes"], rules["missing_hours_threshold"]
MISS_THR_24 = rules.get("window24_missing_hours_threshold", 2)
stations = {s["station_id"]: s for s in catalog["stations"]}

# idempotenza: via le derivate precedenti
catalog["series"] = [s for s in catalog["series"] if not s.get("derived_from")]
rain_series = [s for s in catalog["series"] if s["variable"] == "rain"]


def load_rain(s):
    parts = [pd.read_csv(os.path.join(OUT_DIR, s["data"]["raw"].replace("{year}", str(y))))
             for y in s["data"]["years"]]
    d = pd.concat(parts, ignore_index=True)
    d["time"] = pd.to_datetime(d["time"])
    d = d[d["flag"] < 2].sort_values("time").drop_duplicates("time").reset_index(drop=True)
    return d[["time", "value"]]


def cumulate(d, step_min):
    t = d["time"]
    cy = (t - pd.Timedelta(minutes=1)).dt.year          # anno di cumulo (dato 00:00 del 1/1 -> anno prima)
    d = d.assign(cy=cy.values)
    d["cum"] = np.round(d.groupby("cy")["value"].cumsum(), 2)
    # tempo mancante: buco rispetto al dato precedente (o all'inizio dell'anno per il primo dato)
    first = d["cy"] != d["cy"].shift(1)
    year_start = pd.to_datetime(d["cy"].astype(str) + "-01-01")
    prev = t.shift(1).where(~first, year_start)
    gap_min = (t - prev).dt.total_seconds() / 60
    miss_min = (gap_min - step_min).where(gap_min > GAP_MIN, 0.0).clip(lower=0)
    d["miss_h"] = miss_min.groupby(d["cy"]).cumsum() / 60
    d["flag"] = (d["miss_h"] > MISS_THR_H).astype(np.int8)
    return d


def rolling24(d, step_min):
    """Pioggia nelle 24 h che terminano a ogni timestamp: finestra (t - 24 h, t], sui dati validi."""
    n_exp = WINDOW_H * 60 / step_min                       # dati attesi nella finestra
    roll = d.set_index("time")["value"].rolling(f"{WINDOW_H}h", closed="right")
    out = d.copy()
    out["cum"] = np.round(np.clip(roll.sum().to_numpy(), 0, None), 2) + 0.0     # + 0.0: evita -0.0
    out["miss_h"] = np.maximum(0, n_exp - roll.count().to_numpy()) * step_min / 60
    out["flag"] = (out["miss_h"] > MISS_THR_24).astype(np.int8)
    return out


def aggregate(d, step_min, rule, how="last"):
    expected = {"1h": 60 / step_min, "1d": 1440 / step_min}[rule]
    freq = {"1h": "1h", "1d": "1D"}[rule]
    g = d.set_index("time")
    rs = g.resample(freq, closed="right", label="left")
    out = pd.DataFrame({"value": getattr(rs["cum"], how)(), "n": rs["cum"].count(), "flag": rs["flag"].max()})
    out = out[out["n"] > 0].copy()
    out["coverage"] = (out["n"] / expected).clip(upper=1).round(3)
    out["n"] = out["n"].astype(int)
    out["flag"] = out["flag"].astype(int)
    out["value"] = out["value"].round(2)
    out = out[["value", "n", "coverage", "flag"]]
    out.index = out.index.strftime(TFMT)
    out.index.name = "time"
    return out


def year_summary(d, step_min):
    rows, notes = [], []
    for y, g in d.groupby("cy"):
        miss = float(g["miss_h"].iloc[-1])
        flagged = g[g["flag"] == 1]
        rows.append({"year": int(y), "total_mm": float(g["cum"].iloc[-1]), "n": int(len(g)),
                     "first": g["time"].iloc[0].strftime(TFMT), "last": g["time"].iloc[-1].strftime(TFMT),
                     "missing_hours": round(miss, 1), "complete": bool(len(flagged) == 0)})
        if len(flagged):
            notes.append(f"{y}: circa {miss:.0f} h di dati mancanti nell'anno (prima del primo dato o in lacune > {GAP_MIN} min); "
                         f"dal {flagged['time'].iloc[0].strftime(TFMT)} la cumulata è sottostimata (flag 1, tratteggio). "
                         f"Totale {g['cum'].iloc[-1]:.1f} mm = valore minimo.")
    return rows, notes


os.makedirs(os.path.join(OUT_DIR, "data", "agg"), exist_ok=True)
new_series = []
for s in rain_series:
    sid = f"{s['station_id']}_{SUFFIX}"
    rain = load_rain(s)
    d = cumulate(rain, s["step_min"])
    years_summary, flag_notes = year_summary(d, s["step_min"])

    raw_dir = os.path.join(OUT_DIR, "data", "raw", sid)
    os.makedirs(raw_dir, exist_ok=True)
    years = sorted(d["time"].dt.year.unique().tolist())
    for y in years:
        dy = d[d["time"].dt.year == y]
        pd.DataFrame({"time": dy["time"].dt.strftime(TFMT), "value": dy["cum"], "flag": dy["flag"]}) \
            .to_csv(os.path.join(raw_dir, f"{y}.csv"), index=False)
    for rule in ("1h", "1d"):
        aggregate(d, s["step_min"], rule).to_csv(os.path.join(OUT_DIR, "data", "agg", f"{sid}_{rule}.csv"))

    st = stations[s["station_id"]]
    notes = [f"Serie derivata da {s['series_id']}: pioggia cumulata dal 1 gennaio ore 00:00 (si azzera ogni anno), "
             f"calcolata sui dati con flag < 2. Il dato delle 00:00 del 1 gennaio chiude l'anno precedente.",
             "Totali annui: " + "; ".join(f"{r['year']}: {r['total_mm']:.1f} mm" + ("" if r["complete"] else " (incompleto)")
                                           for r in years_summary) + "."] + flag_notes
    counts = d["flag"].value_counts().to_dict()
    new_series.append({
        "series_id": sid, "instrument_id": s["instrument_id"], "station_id": s["station_id"],
        "variable": "rain_cum", "depth_m": None, "sensor_model": None, "logger_port": None,
        "label": f"{VARIABLE['label_it']} – {st['name']}",
        "unit": "mm", "aggregation": "last", "step_min": s["step_min"],
        "start": d["time"].iloc[0].strftime(ISO), "end": d["time"].iloc[-1].strftime(ISO),
        "n_values": int(len(d)), "min": 0.0, "max": round(float(d["cum"].max()), 2),
        "status": "active",
        "gaps_over_24h": s["gaps_over_24h"],
        "qc_notes": notes, "sources": [],
        "derived_from": s["series_id"],
        "cumulative_years": years_summary,
        "data": {"format": "csv: time,value,flag", "raw": f"data/raw/{sid}/{{year}}.csv", "years": years,
                 "agg_1h": f"data/agg/{sid}_1h.csv", "agg_1d": f"data/agg/{sid}_1d.csv"},
        "qc_summary": {"n_total": int(len(d)), "n_valid": int(counts.get(0, 0)),
                       "n_suspect": int(counts.get(1, 0)), "n_bad": 0, "duplicates_removed": 0},
    })
    print(f"{sid:20s} righe={len(d):7d}  flag1={counts.get(1, 0):7d}  " +
          "  ".join(f"{r['year']}: {r['total_mm']:.1f} mm (lacune {r['missing_hours']} h)" for r in years_summary))

    # ---- pioggia cumulata 24 h (finestra mobile)
    sid24 = f"{s['station_id']}_{SUFFIX24}"
    d24 = rolling24(rain, s["step_min"])
    raw_dir24 = os.path.join(OUT_DIR, "data", "raw", sid24)
    os.makedirs(raw_dir24, exist_ok=True)
    years24 = sorted(d24["time"].dt.year.unique().tolist())
    for y in years24:
        dy = d24[d24["time"].dt.year == y]
        pd.DataFrame({"time": dy["time"].dt.strftime(TFMT), "value": dy["cum"], "flag": dy["flag"]}) \
            .to_csv(os.path.join(raw_dir24, f"{y}.csv"), index=False)
    for rule in ("1h", "1d"):
        aggregate(d24, s["step_min"], rule, how="max").to_csv(
            os.path.join(OUT_DIR, "data", "agg", f"{sid24}_{rule}.csv"))
    # massimo di 24 h per anno (solo dati non sottostimati, se ce ne sono)
    max_rows = []
    for y, g in d24.groupby(d24["time"].dt.year):
        ok = g[g["flag"] == 0]
        gg = ok if len(ok) else g
        i = gg["cum"].idxmax()
        max_rows.append({"year": int(y), "max_mm": float(gg.loc[i, "cum"]), "time": gg.loc[i, "time"].strftime(TFMT),
                         "n": int(len(g)), "n_underestimated": int((g["flag"] == 1).sum())})
    c24 = d24["flag"].value_counts().to_dict()
    notes24 = [f"Serie derivata da {s['series_id']}: pioggia caduta nelle {WINDOW_H} h che terminano al timestamp "
               f"(finestra mobile (t − {WINDOW_H} h, t]), calcolata sui dati con flag < 2.",
               f"Flag 1 (tratteggio) dove nella finestra mancano più di {MISS_THR_24:g} h di dati, anche nelle prime "
               f"{WINDOW_H} h della serie: la cumulata è un valore minimo.",
               f"Aggregati orario e giornaliero: massimo della cumulata {WINDOW_H} h nell'intervallo.",
               f"Massimo {WINDOW_H} h: " + "; ".join(f"{r['year']}: {r['max_mm']:.1f} mm ({r['time']})" for r in max_rows) + "."]
    new_series.append({
        "series_id": sid24, "instrument_id": s["instrument_id"], "station_id": s["station_id"],
        "variable": "rain_cum24", "depth_m": None, "sensor_model": None, "logger_port": None,
        "label": f"{VARIABLE24['label_it']} – {st['name']}",
        "unit": "mm", "aggregation": "max", "step_min": s["step_min"], "window_hours": WINDOW_H,
        "start": d24["time"].iloc[0].strftime(ISO), "end": d24["time"].iloc[-1].strftime(ISO),
        "n_values": int(len(d24)), "min": 0.0, "max": round(float(d24["cum"].max()), 2),
        "status": "active",
        "gaps_over_24h": s["gaps_over_24h"],
        "qc_notes": notes24, "sources": [],
        "derived_from": s["series_id"],
        "max_24h_by_year": max_rows,
        "data": {"format": "csv: time,value,flag", "raw": f"data/raw/{sid24}/{{year}}.csv", "years": years24,
                 "agg_1h": f"data/agg/{sid24}_1h.csv", "agg_1d": f"data/agg/{sid24}_1d.csv"},
        "qc_summary": {"n_total": int(len(d24)), "n_valid": int(c24.get(0, 0)),
                       "n_suspect": int(c24.get(1, 0)), "n_bad": 0, "duplicates_removed": 0},
    })
    print(f"{sid24:20s} righe={len(d24):7d}  flag1={c24.get(1, 0):7d}  " +
          "  ".join(f"{r['year']}: max {r['max_mm']:.1f} mm ({r['time']})" for r in max_rows))

catalog["series"] = sorted(catalog["series"] + new_series, key=lambda x: x["series_id"])
catalog["variables"]["rain_cum"] = VARIABLE
catalog["variables"]["rain_cum24"] = VARIABLE24
catalog["catalog_version"] = "0.5"
catalog["conventions"]["derived_series"] = (
    "Serie con 'derived_from' sono calcolate da altre serie (pipeline/derive_cumulative.py). "
    "Pioggia cumulata (variabile rain_cum): cumulo dal 1 gennaio ore 00:00; il dato di pioggia delle 00:00 "
    "del 1 gennaio chiude l'anno precedente. Aggregati (data/agg): value = cumulata a fine intervallo, "
    "etichetta = inizio intervallo; colonne time,value,n,coverage,flag. flag 1 = cumulata sottostimata "
    "per dati mancanti (soglia in qc_rules.json -> cumulative). "
    "Pioggia cumulata 24 h (variabile rain_cum24): pioggia nelle 24 h che terminano al timestamp (finestra mobile "
    "(t-24 h, t], non si azzera); aggregati orario e giornaliero = massimo della cumulata 24 h nell'intervallo, "
    "etichetta = inizio intervallo; flag 1 = finestra con più di window24_missing_hours_threshold ore di dati mancanti "
    "(anche nelle prime 24 h della serie).")
catalog["generated"] = pd.Timestamp.now(tz="Etc/GMT-1").strftime("%Y-%m-%dT%H:%M:%S+01:00")
with open(CATALOG, "w", encoding="utf-8") as fh:
    json.dump(catalog, fh, ensure_ascii=False, indent=2)

# series.csv: vista piatta (righe delle derivate rigenerate)
flat_path = os.path.join(OUT_DIR, "series.csv")
if os.path.exists(flat_path):
    flat = pd.read_csv(flat_path, encoding="utf-8-sig")
    flat = flat[~flat["series_id"].isin([x["series_id"] for x in new_series])]
    rows = [{**{k: x[k] for k in ("series_id", "station_id", "instrument_id", "variable", "label", "unit",
                                  "depth_m", "sensor_model", "logger_port", "aggregation", "step_min",
                                  "start", "end", "n_values", "min", "max", "status")},
             "lat": stations[x["station_id"]]["lat"], "lon": stations[x["station_id"]]["lon"],
             "n_gaps_over_24h": len(x["gaps_over_24h"]),
             "qc_notes": " | ".join(x["qc_notes"]), "source_files": ""} for x in new_series]
    flat = pd.concat([flat, pd.DataFrame(rows)], ignore_index=True)[flat.columns]
    flat.sort_values("series_id").to_csv(flat_path, index=False, encoding="utf-8-sig")
