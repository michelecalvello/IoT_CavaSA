"""
standardize.py — passo 2: conversione delle serie nel formato dati standard.

Legge   : catalog.json (prodotto da build_catalog.py), qc_rules.json, file Excel in INPUT_DIR
Scrive  : OUT_DIR/data/raw/<series_id>/<anno>.csv   time,value,flag      (dato originale)
          OUT_DIR/data/agg/<series_id>_1h.csv        aggregati orari
          OUT_DIR/data/agg/<series_id>_1d.csv        aggregati giornalieri
          OUT_DIR/data/qc_log.csv                    riepilogo dei valori segnalati
          OUT_DIR/catalog.json                       arricchito con percorsi file e statistiche QC

FORMATO STANDARD
- time  : "YYYY-MM-DD HH:MM" in UTC+1 (ora solare, senza ora legale), senza offset nel testo:
          così il browser non applica conversioni di fuso e mostra l'ora così com'è.
- value : valore numerico nell'unità indicata dal catalogo
- flag  : 0 valido, 1 sospetto, 2 errato (vedi qc_rules.json)

AGGREGATI (solo valori con flag < 2)
- grandezze istantanee (aggregation = mean): media, min, max; etichetta = inizio intervallo
  (es. "2024-05-28 10:00" = 10:00–10:59)
- pioggia (aggregation = sum): somma; ogni dato rappresenta l'intervallo che TERMINA al suo
  timestamp, quindi l'ora "10:00" comprende i dati da 10:00 (escluso) a 11:00 (incluso)
- n = numero di valori usati; coverage = n / valori attesi (in base al passo nominale)
"""
import json, os, re, glob
import numpy as np
import pandas as pd

INPUT_DIR = os.environ.get("INPUT_DIR", ".")
OUT_DIR = os.environ.get("OUT_DIR", ".")
CATALOG = os.path.join(OUT_DIR, "catalog.json")
RULES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "qc_rules.json")
TFMT = "%Y-%m-%d %H:%M"

catalog = json.load(open(CATALOG, encoding="utf-8"))
rules = json.load(open(RULES, encoding="utf-8"))

# mappa nome file "pulito" -> percorso reale (gestisce eventuali prefissi di upload)
paths = {re.sub(r"^[0-9a-f]{8}-", "", os.path.basename(p)): p
         for p in glob.glob(os.path.join(INPUT_DIR, "*.xlsx"))}
_sheet_cache = {}
def read_sheet(fname, sheet):
    key = (fname, sheet)
    if key not in _sheet_cache:
        _sheet_cache[key] = pd.read_excel(paths[fname], sheet_name=sheet, engine="openpyxl")
    return _sheet_cache[key]

def load_series(s):
    parts = []
    for src in s["sources"]:
        df = read_sheet(src["file"], src["sheet"])
        part = pd.DataFrame({
            "time": pd.to_datetime(df.iloc[:, 1], errors="coerce").dt.round("min"),   # colonna B = 'Time UTC+1' (o 'Time (UTC+1)')
            "value": pd.to_numeric(df[src["column"]], errors="coerce")})
        parts.append(part.dropna())
    d = pd.concat(parts).sort_values("time")
    n_dup = int(d["time"].duplicated().sum())
    d = d.drop_duplicates("time", keep="first").reset_index(drop=True)
    return d, n_dup

def apply_qc(s, d, log):
    sid, var = s["series_id"], s["variable"]
    flag = np.zeros(len(d), dtype=np.int8)

    def mark(mask, f, rule):
        mask = np.asarray(mask) & (flag < f)
        if mask.any():
            flag[mask] = f
            t = d["time"][mask]
            log.append({"series_id": sid, "rule": rule, "flag": f, "n": int(mask.sum()),
                        "first": t.min().strftime(TFMT), "last": t.max().strftime(TFMT)})

    r = rules["range_by_variable"].get(var)
    if r:
        mark((d["value"] < r["min"]) | (d["value"] > r["max"]), 2, f"fuori range [{r['min']}, {r['max']}]")
    codes = rules["error_codes"].get(sid)
    if codes:
        mark(d["value"].isin(codes), 2, f"codice di errore {codes}")
    for w in rules["windows"]:
        if sid in w["series"]:
            m = (d["time"] >= pd.Timestamp(w["from"])) & (d["time"] <= pd.Timestamp(w["to"]))
            mark(m, w["flag"], w["reason"])
    sp = rules["spike_by_series"].get(sid)
    if sp:
        # picco isolato: il valore si discosta da ENTRAMBI i vicini oltre la soglia,
        # mentre i due vicini sono tra loro coerenti (una piena reale sale/scende con continuità)
        v = d["value"].where(flag < 2)
        prev, nxt = v.shift(1), v.shift(-1)
        isolated = ((v - prev).abs() > sp["threshold"]) & ((v - nxt).abs() > sp["threshold"]) \
                   & ((prev - nxt).abs() <= sp["neighbor_tolerance"])
        mark(isolated, 1, f"picco isolato > {sp['threshold']} rispetto ai valori adiacenti")
    d["flag"] = flag
    return d

def aggregate(s, d, rule):
    step = s["step_min"]
    expected = {"1h": 60 / step, "1d": 1440 / step}[rule]
    freq = {"1h": "1h", "1d": "1D"}[rule]
    g = d[d["flag"] < 2].set_index("time")["value"]
    if s["aggregation"] == "sum":
        rs = g.resample(freq, closed="right", label="left")
        out = pd.DataFrame({"value": rs.sum(min_count=1), "max": rs.max(), "n": rs.count()})
    else:
        rs = g.resample(freq)
        out = pd.DataFrame({"value": rs.mean(), "min": rs.min(), "max": rs.max(), "n": rs.count()})
    out = out[out["n"] > 0]
    out["coverage"] = (out["n"] / expected).clip(upper=1).round(3)
    out["n"] = out["n"].astype(int)
    out.index = out.index.strftime(TFMT)
    out.index.name = "time"
    return out

def fmt_value(x, var):
    return x.round({"vwc": 4, "rain": 2}.get(var, 2))

log, processed = [], []
os.makedirs(os.path.join(OUT_DIR, "data", "agg"), exist_ok=True)
for s in catalog["series"]:
    if s.get("derived_from"):       # serie derivate (derive_cumulative.py): non hanno file Excel sorgente
        continue
    missing = sorted({x["file"] for x in s["sources"] if x["file"] not in paths})
    if missing:                     # modalità incrementale: serie già elaborate, Excel non presenti in INPUT_DIR
        if "data" not in s:
            raise SystemExit(f"{s['series_id']}: file Excel mancanti in INPUT_DIR ({', '.join(missing)}) e serie non ancora elaborata")
        print(f"{s['series_id']:18s} saltata (Excel non presenti, dati già elaborati)")
        continue
    sid, var = s["series_id"], s["variable"]
    d, n_dup = load_series(s)
    if n_dup:
        log.append({"series_id": sid, "rule": "timestamp duplicati rimossi", "flag": None,
                    "n": n_dup, "first": None, "last": None})
    d = apply_qc(s, d, log)

    # lacune > 24 h sulla serie completa (anche tra un file annuale e il successivo)
    dt = d["time"].diff()
    s["gaps_over_24h"] = [{"from": d["time"][i - 1].strftime("%Y-%m-%dT%H:%M+01:00"),
                           "to": d["time"][i].strftime("%Y-%m-%dT%H:%M+01:00"),
                           "hours": round(dt[i].total_seconds() / 3600, 1)}
                          for i in dt.index[dt > pd.Timedelta(hours=24)]]

    # dati originali, un file per anno
    raw_dir = os.path.join(OUT_DIR, "data", "raw", sid)
    os.makedirs(raw_dir, exist_ok=True)
    years = sorted(d["time"].dt.year.unique().tolist())
    for y in years:
        dy = d[d["time"].dt.year == y].copy()
        dy["time"] = dy["time"].dt.strftime(TFMT)
        dy["value"] = fmt_value(dy["value"], var)
        dy.to_csv(os.path.join(raw_dir, f"{y}.csv"), index=False)

    # aggregati
    for rule in ("1h", "1d"):
        a = aggregate(s, d, rule)
        for c in ("value", "min", "max"):
            if c in a:
                a[c] = fmt_value(a[c], var)
        a.to_csv(os.path.join(OUT_DIR, "data", "agg", f"{sid}_{rule}.csv"))

    processed.append(sid)
    counts = d["flag"].value_counts().to_dict()
    s["data"] = {
        "format": "csv: time,value,flag",
        "raw": f"data/raw/{sid}/{{year}}.csv",
        "years": years,
        "agg_1h": f"data/agg/{sid}_1h.csv",
        "agg_1d": f"data/agg/{sid}_1d.csv",
    }
    s["qc_summary"] = {"n_total": int(len(d)),
                       "n_valid": int(counts.get(0, 0)),
                       "n_suspect": int(counts.get(1, 0)),
                       "n_bad": int(counts.get(2, 0)),
                       "duplicates_removed": n_dup}
    print(f"{sid:18s} n={len(d):7d}  sospetti={counts.get(1,0):4d}  errati={counts.get(2,0):4d}  dup={n_dup}  anni={years}")

catalog["catalog_version"] = "0.3"
catalog["conventions"]["data_format"] = (
    "File CSV con intestazione. Dati originali: time,value,flag (un file per serie e per anno). "
    "Aggregati: time,value,[min],max,n,coverage. time = 'YYYY-MM-DD HH:MM' in UTC+1. "
    "flag: 0 valido, 1 sospetto, 2 errato. Aggregati calcolati sui soli flag < 2.")
catalog["conventions"]["aggregation_labels"] = (
    "Grandezze istantanee: etichetta = inizio intervallo. Pioggia: ogni dato è riferito "
    "all'intervallo che termina al timestamp; l'etichetta oraria/giornaliera indica l'inizio "
    "del periodo di cumulo (es. giorno '2024-10-01' = pioggia da 2024-10-01 00:00 esclusa "
    "a 2024-10-02 00:00 inclusa).")
catalog["qc_rules"] = {k: v for k, v in rules.items() if not k.startswith("_")}
with open(CATALOG, "w", encoding="utf-8") as fh:
    json.dump(catalog, fh, ensure_ascii=False, indent=2)
_qc_path = os.path.join(OUT_DIR, "data", "qc_log.csv")
_log = pd.DataFrame(log)
if os.path.exists(_qc_path):       # modalità incrementale: si conservano le righe delle serie non rielaborate
    _old = pd.read_csv(_qc_path)
    _log = pd.concat([_old[~_old["series_id"].isin(processed)], _log], ignore_index=True)
    _log = _log.sort_values("series_id", kind="stable").reset_index(drop=True)
_log.to_csv(_qc_path, index=False)
# allinea series.csv (numero di lacune ricalcolato sulla serie completa)
_flat_path = os.path.join(OUT_DIR, "series.csv")
if os.path.exists(_flat_path):
    _flat = pd.read_csv(_flat_path, encoding="utf-8-sig")
    _ng = {s["series_id"]: len(s["gaps_over_24h"]) for s in catalog["series"]}
    _flat["n_gaps_over_24h"] = _flat["series_id"].map(_ng)
    _flat.to_csv(_flat_path, index=False, encoding="utf-8-sig")
print("\nQC log:"); print(pd.DataFrame(log).to_string())
