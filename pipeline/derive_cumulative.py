"""
derive_cumulative.py — passo 3: serie derivate di pioggia cumulata dal 1 gennaio.

Legge   : OUT_DIR/catalog.json (dopo standardize.py), qc_rules.json (blocco "cumulative"),
          OUT_DIR/data/raw/<pluviometro>/<anno>.csv
Scrive  : OUT_DIR/data/raw/<STAZIONE>_RAINCUM/<anno>.csv     time,value,flag
          OUT_DIR/data/agg/<STAZIONE>_RAINCUM_1h.csv / _1d.csv   time,value,n,coverage,flag
Aggiorna: OUT_DIR/catalog.json e OUT_DIR/series.csv (le serie derivate hanno "derived_from").
Lo script è idempotente: a ogni esecuzione rigenera tutte le serie derivate.

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


def aggregate(d, step_min, rule):
    expected = {"1h": 60 / step_min, "1d": 1440 / step_min}[rule]
    freq = {"1h": "1h", "1d": "1D"}[rule]
    g = d.set_index("time")
    rs = g.resample(freq, closed="right", label="left")
    out = pd.DataFrame({"value": rs["cum"].last(), "n": rs["cum"].count(), "flag": rs["flag"].max()})
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
    d = cumulate(load_rain(s), s["step_min"])
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

catalog["series"] = sorted(catalog["series"] + new_series, key=lambda x: x["series_id"])
catalog["variables"]["rain_cum"] = VARIABLE
catalog["catalog_version"] = "0.4"
catalog["conventions"]["derived_series"] = (
    "Serie con 'derived_from' sono calcolate da altre serie (pipeline/derive_cumulative.py). "
    "Pioggia cumulata (variabile rain_cum): cumulo dal 1 gennaio ore 00:00; il dato di pioggia delle 00:00 "
    "del 1 gennaio chiude l'anno precedente. Aggregati (data/agg): value = cumulata a fine intervallo, "
    "etichetta = inizio intervallo; colonne time,value,n,coverage,flag. flag 1 = cumulata sottostimata "
    "per dati mancanti (soglia in qc_rules.json -> cumulative).")
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
