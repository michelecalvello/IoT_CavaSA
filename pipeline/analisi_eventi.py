"""
analisi_eventi.py — analisi eventi pluviometrici / immagazzinamento idrico 0–75 cm (sito V e NV).

Legge   : data/raw/<serie>/<anno>.csv (solo flag < 2)
Scrive  : data/analisi_eventi.json   (letto da analisi.html)

Tempi: UTC+1 come nel resto del portale.
Pioggia di riferimento: media dei pluviometri CF 21521, 18925, 18957 (passo 10 min; dato riferito all'intervallo
che termina al timestamp). Evento = pioggia separata dalla precedente/successiva da >= GAP_H ore asciutte.
Stoccaggio S(t) [mm] = somma_i VWC_i * spessore_i * 10, spessori da punti medi tra sensori, fino a 75 cm:
    NV (15, 30, 60 cm): 0–22.5, 22.5–45, 45–75       V (15, 35, 60 cm): 0–25, 25–47.5, 47.5–75
"""
import os, json, glob
import numpy as np, pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
YEAR = "2025"
GAP_H, MIN_P = 6, 5.0
POST_H = 48
GAUGES = ["CF_21521", "CF_18925", "CF_18957"]
SITES = {
    "NV": {"station": "SAL_NV", "depths": ["015", "030", "060"], "bounds": [0, 22.5, 45, 75]},
    "V":  {"station": "SAL_V",  "depths": ["015", "035", "060"], "bounds": [0, 25, 47.5, 75]},
}

def rawv(sid):
    fs = sorted(glob.glob(os.path.join(ROOT, "data", "raw", sid, "*.csv")))
    d = pd.concat([pd.read_csv(f) for f in fs])
    d["time"] = pd.to_datetime(d["time"])
    return d[d["flag"] < 2].set_index("time")["value"].sort_index()

# ---------------- pioggia ----------------
R = pd.DataFrame({g: rawv(g + "_RAIN") for g in GAUGES}).loc[YEAR].resample("10min").sum(min_count=1)
Rm = R.mean(axis=1)
wet_t = Rm.index[(Rm.fillna(0) > 0).values]
ev, s, last = [], wet_t[0], wet_t[0]
for x in wet_t[1:]:
    if x - last > pd.Timedelta(hours=GAP_H):
        ev.append((s, last)); s = x
    last = x
ev.append((s, last))
events = pd.DataFrame(ev, columns=["t_first", "t_last"])
events["start"] = events.t_first - pd.Timedelta(minutes=10)      # inizio del primo intervallo bagnato
events["end"] = events.t_last
events["next_start"] = events.start.shift(-1)
events["P"] = [Rm[a:b].sum() for a, b in zip(events.t_first, events.t_last)]
events["P_near"] = [R["CF_18957"][a:b].sum() for a, b in zip(events.t_first, events.t_last)]
events["P_min"] = [R[list(GAUGES)][a:b].sum().min() for a, b in zip(events.t_first, events.t_last)]
events["P_max"] = [R[list(GAUGES)][a:b].sum().max() for a, b in zip(events.t_first, events.t_last)]
events["dur_h"] = (events.end - events.start).dt.total_seconds() / 3600
events["I10"] = [Rm[a:b].max() * 6 for a, b in zip(events.t_first, events.t_last)]                  # mm/h su 10 min
events["I30"] = [Rm[a:b].rolling(3, min_periods=1).sum().max() * 2 for a, b in zip(events.t_first, events.t_last)]
events["I60"] = [Rm[a:b].rolling(6, min_periods=1).sum().max() for a, b in zip(events.t_first, events.t_last)]
events["Imean"] = events.P / events.dur_h
events["API7"] = [Rm[a - pd.Timedelta(days=7):a - pd.Timedelta(minutes=10)].sum() for a in events.start]
events["API3"] = [Rm[a - pd.Timedelta(days=3):a - pd.Timedelta(minutes=10)].sum() for a in events.start]
events = events[events.P >= MIN_P].reset_index(drop=True)

# ---------------- suolo ----------------
VW, S5 = {}, {}
for site, c in SITES.items():
    cols = {}
    for dp in c["depths"]:
        cols[dp] = rawv(f"{c['station']}_VWC_{dp}").loc[YEAR].resample("5min").mean().interpolate(limit=3)
    V = pd.DataFrame(cols)
    thick = np.diff(c["bounds"])
    S = sum(V[dp] * t * 10 for dp, t in zip(c["depths"], thick))
    VW[site] = V; S5[site] = S
SMAX = {k: float(v.quantile(0.995)) for k, v in S5.items()}
SMIN = {k: float(v.quantile(0.005)) for k, v in S5.items()}

def med(x, a, b):
    w = x[a:b].dropna(); return float(w.median()) if len(w) >= 2 else np.nan

rows = []
for _, e in events.iterrows():
    nxt = e.next_start if pd.notna(e.next_start) else e.end + pd.Timedelta(hours=POST_H)
    we = min(e.end + pd.Timedelta(hours=POST_H), nxt)
    post_h = (we - e.end).total_seconds() / 3600
    for site, c in SITES.items():
        S, V = S5[site], VW[site]
        s0 = med(S, e.start - pd.Timedelta(minutes=30), e.start)
        w = S[e.start:we].dropna()
        if np.isnan(s0) or len(w) < 20: continue
        sp, tp = float(w.max()), w.idxmax()
        send = med(S, we - pd.Timedelta(minutes=30), we)
        th = np.diff(c["bounds"])
        lay = {}
        for dp, t in zip(c["depths"], th):
            b = med(V[dp], e.start - pd.Timedelta(minutes=30), e.start)
            m = float(V[dp][e.start:we].max())
            lay[dp] = round((m - b) * t * 10, 2)
        rows.append(dict(
            ev=e.start.strftime("%Y-%m-%d %H:%M"), site=site, P=round(e.P, 1), P_near=round(e.P_near, 1),
            P_min=round(e.P_min, 1), P_max=round(e.P_max, 1),
            dur_h=round(e.dur_h, 1), I10=round(e.I10, 1), I30=round(e.I30, 1), I60=round(e.I60, 1),
            Imean=round(e.Imean, 2), API7=round(e.API7, 1), API3=round(e.API3, 1),
            S0=round(s0, 1), D0=round(SMAX[site] - s0, 1), f0=round(s0 / SMAX[site], 3),
            dS=round(sp - s0, 1), dS_end=round(send - s0, 1) if not np.isnan(send) else None,
            t_peak_h=round((tp - e.start).total_seconds() / 3600, 1), post_h=round(post_h, 1),
            eff=round((sp - s0) / e.P, 3), loss=round(e.P - (sp - s0), 1),
            layers=lay, month=int(e.start.month)))
EV = pd.DataFrame(rows)

# ---------------- bilancio mensile e svuotamento in tempo asciutto ----------------
Rd = Rm.resample("1D", closed="right", label="left").sum(min_count=100)
monthly = {}
for site in SITES:
    Sd = S5[site].resample("1D").mean()
    out = []
    for m in range(1, 13):
        a, b = pd.Timestamp(f"{YEAR}-{m:02d}-01"), pd.Timestamp(f"{YEAR}-{m:02d}-01") + pd.offsets.MonthEnd(0)
        s_a = float(Sd[a - pd.Timedelta(days=1):a + pd.Timedelta(days=1)].mean()) if m > 1 else np.nan
        s_b = float(Sd[b - pd.Timedelta(days=1):b].mean())
        s_a = s_a if not np.isnan(s_a) else float(Sd[a:a + pd.Timedelta(days=2)].mean())
        P = float(Rd[a:b].sum())
        out.append(dict(month=m, P=round(P, 1), S_start=round(s_a, 1), S_end=round(s_b, 1),
                        dS=round(s_b - s_a, 1), loss=round(P - (s_b - s_a), 1)))
    monthly[site] = out

# perdita di stoccaggio in tempo asciutto: giorno d con pioggia ~0 in d-2..d e anche in d+1
# (variazione giornaliera dS = S(d+1) - S(d), medie giornaliere); include drenaggio + evapotraspirazione
dry = (Rd.rolling(3).sum() < 0.3)
dry_ok = dry & dry.shift(-1).fillna(False)
depl, dry_pts = {}, {}
for site in SITES:
    Sd = S5[site].resample("1D").mean()
    rate = (Sd.shift(-1) - Sd)[dry_ok].dropna()
    f = ((Sd - SMIN[site]) / (SMAX[site] - SMIN[site])).reindex(rate.index)
    dry_pts[site] = [[d.strftime("%Y-%m-%d"), round(float(Sd[d]), 1), round(float(r), 2)] for d, r in rate.items()]
    cls = [(0, .25), (.25, .5), (.5, .75), (.75, 1.5)]
    out = []
    for lo, hi in cls:
        q = rate[(f >= lo) & (f < hi)]
        if len(q) >= 3:
            out.append(dict(f_lo=lo, f_hi=min(hi, 1.0), n=int(len(q)), rate=round(float(q.median()), 2),
                            p25=round(float(q.quantile(.25)), 2), p75=round(float(q.quantile(.75)), 2)))
    depl[site] = out

# ---------------- serie giornaliere ----------------
idx = pd.date_range(f"{YEAR}-02-01", f"{YEAR}-12-31", freq="D")
daily = {"t": [d.strftime("%Y-%m-%d") for d in idx],
         "rain": [None if np.isnan(v) else round(float(v), 1) for v in Rd.reindex(idx)]}
for site in SITES:
    Sd = S5[site].resample("1D").mean().reindex(idx)
    daily["S_" + site] = [None if np.isnan(v) else round(float(v), 1) for v in Sd]

# ---------------- finestre orarie per eventi >= 10 mm ----------------
hourly = {}
Rh = Rm.resample("1h", closed="right", label="left").sum(min_count=1)
for _, e in events[events.P >= 10].iterrows():
    nxt = e.next_start if pd.notna(e.next_start) else e.end + pd.Timedelta(hours=POST_H)
    a, b = e.start.floor("h") - pd.Timedelta(hours=6), min(e.end + pd.Timedelta(hours=POST_H), nxt).ceil("h")
    h = pd.date_range(a, b, freq="1h")
    rec = {"t": [t.strftime("%Y-%m-%d %H:%M") for t in h],
           "rain": [None if np.isnan(v) else round(float(v), 2) for v in Rh.reindex(h)]}
    for site, c in SITES.items():
        for dp in c["depths"]:
            v = VW[site][dp].resample("1h").mean().reindex(h)
            rec[f"{site}_{dp}"] = [None if np.isnan(x) else round(float(x), 3) for x in v]
        sh = S5[site].resample("1h").mean().reindex(h)
        rec["S_" + site] = [None if np.isnan(x) else round(float(x), 1) for x in sh]
    hourly[e.start.strftime("%Y-%m-%d %H:%M")] = rec

# ---------------- statistiche ----------------
stats = {}
for site in SITES:
    d = EV[EV.site == site]
    st = {"n": int(len(d)), "Smax": round(SMAX[site], 1), "Smin": round(SMIN[site], 1)}
    st["spearman_dS"] = {k: round(float(d[k].corr(d["dS"], method="spearman")), 2)
                         for k in ["P", "P_near", "I10", "I30", "I60", "dur_h", "Imean", "API7", "S0", "D0"]}
    st["spearman_eff"] = {k: round(float(d[k].corr(d["eff"], method="spearman")), 2)
                          for k in ["P", "I10", "I60", "dur_h", "API7", "S0", "D0"]}
    big = d[d.P >= 10]
    st["eff_median_P10"] = round(float(big.eff.median()), 2)
    st["eff_sum"] = round(float(d.dS.sum() / d.P.sum()), 2)
    # modelli: (1) dS = k*P ; (2) secchio dS = min(k*P, D0)
    def r2(y, yh): return round(float(1 - ((y - yh) ** 2).sum() / ((y - y.mean()) ** 2).sum()), 2)
    k1 = float((d.P * d.dS).sum() / (d.P ** 2).sum())
    best = None
    for k in np.arange(0.2, 1.51, 0.01):
        yh = np.minimum(k * d.P, d.D0)
        sse = ((d.dS - yh) ** 2).sum()
        if best is None or sse < best[0]: best = (sse, k)
    k2 = float(best[1])
    st["model_lin"] = {"k": round(k1, 2), "r2": r2(d.dS, k1 * d.P)}
    st["model_bucket"] = {"k": round(k2, 2), "r2": r2(d.dS, np.minimum(k2 * d.P, d.D0))}
    # regressione multipla su dS: P, I60, D0 (standardizzati) -- indicativa (n piccolo)
    X = d[["P", "I60", "D0"]].astype(float); y = d.dS.astype(float)
    Z = (X - X.mean()) / X.std(); A = np.column_stack([np.ones(len(Z)), Z.values])
    beta = np.linalg.lstsq(A, (y - y.mean()) / y.std(), rcond=None)[0]
    st["beta_std"] = {"P": round(float(beta[1]), 2), "I60": round(float(beta[2]), 2), "D0": round(float(beta[3]), 2)}
    st["r2_multi"] = r2((y - y.mean()) / y.std(), A @ beta)
    # residuo del modello a secchio vs caratteristiche della pioggia (a parita' di P e D0)
    res = d.dS - np.minimum(k2 * d.P, d.D0)
    st["resid_corr"] = {k: round(float(d[k].corr(res, method="spearman")), 2) for k in ["I10", "I30", "I60", "dur_h", "API7", "S0"]}
    mid = d[(d.P >= 10) & (d.P <= 40) & (d.D0 > 25)]
    st["mid_events"] = {"n": int(len(mid)),
                        "eff_I60_spearman": round(float(mid.eff.corr(mid.I60, method="spearman")), 2),
                        "eff_dur_spearman": round(float(mid.eff.corr(mid.dur_h, method="spearman")), 2),
                        "eff_S0_spearman": round(float(mid.eff.corr(mid.S0, method="spearman")), 2)}
    resp = d[d.dS >= 5]
    st["resp_n"] = int(len(resp)); st["resp_P_min"] = float(resp.P.min()) if len(resp) else None
    st["tpeak_med_h"] = round(float(resp.t_peak_h.median()), 1) if len(resp) else None
    big2 = d[d.P >= 20]
    st["P_ge20"] = {"n": int(len(big2)), "dS_med": round(float(big2.dS.median()), 1), "dS_max": round(float(big2.dS.max()), 1),
                    "eff_med": round(float(big2.eff.median()), 2)}
    # coefficiente di variazione dell'incremento nei grandi eventi e quota dell'incremento per strato
    lay = pd.DataFrame(list(d[d.dS >= 10].layers)); 
    if len(lay): st["layer_share"] = {c: round(float(lay[c].sum() / lay.sum().sum()), 2) for c in lay.columns}
    # quanto resta dopo il drenaggio: eventi P >= 10 mm con almeno 24 h di finestra dopo la fine della pioggia
    ret = d[(d.P >= 10) & (d.post_h >= 24) & d.dS_end.notna()]
    big_r = ret[ret.dS > 3]
    st["retained"] = {"n": int(len(ret)), "eff_end_sum": round(float(ret.dS_end.sum() / ret.P.sum()), 2),
                      "eff_peak_sum": round(float(ret.dS.sum() / ret.P.sum()), 2),
                      "frac_kept_med": round(float((big_r.dS_end / big_r.dS).median()), 2) if len(big_r) else None,
                      "n_frac": int(len(big_r))}
    # per stagione
    seas = {"Feb–Apr": [2, 3, 4], "Mag–Ago": [5, 6, 7, 8], "Set–Nov": [9, 10, 11], "Dic": [12]}
    ss = {}
    for nm, ms in seas.items():
        q = d[d.month.isin(ms)]
        if len(q): ss[nm] = {"n": int(len(q)), "P": round(float(q.P.sum()), 0), "dS": round(float(q.dS.sum()), 0),
                              "eff": round(float(q.dS.sum() / q.P.sum()), 2), "S0_med": round(float(q.S0.median()), 0)}
    st["seasons"] = ss
    stats[site] = st

out = {"generated": pd.Timestamp.now(tz="Etc/GMT-1").strftime("%Y-%m-%dT%H:%M:%S+01:00"),
       "params": {"gap_h": GAP_H, "min_P": MIN_P, "post_h": POST_H, "year": int(YEAR),
                  "bounds": {k: v["bounds"] for k, v in SITES.items()},
                  "depths": {k: v["depths"] for k, v in SITES.items()}},
       "events": EV.to_dict(orient="records"), "monthly": monthly, "depletion": depl, "dry_pts": dry_pts,
       "daily": daily, "hourly": hourly, "stats": stats}
def clean(o):                                   # NaN non e' JSON valido -> null
    if isinstance(o, dict): return {k: clean(v) for k, v in o.items()}
    if isinstance(o, list): return [clean(v) for v in o]
    if isinstance(o, float) and o != o: return None
    return o
with open(os.path.join(ROOT, "data", "analisi_eventi.json"), "w", encoding="utf-8") as fh:
    json.dump(clean(out), fh, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
print("eventi", len(events), "righe evento×sito", len(EV), "| finestre orarie", len(hourly))
print(json.dumps(stats, ensure_ascii=False, indent=1))
