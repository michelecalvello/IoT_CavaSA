"""
curve_caratteristiche.py — curve caratteristiche θ–ψ in situ nelle fasi di drying e di wetting (siti Meter NV e V).

Legge   : data/raw/SAL_{NV,V}_{VWC,PSI}_<prof>/*.csv  (6 coppie: 2 siti × 3 profondità, solo flag < 2)
          data/raw/CF_{21521,18925,18957}_RAIN/*.csv  (pioggia: media dei 3 pluviometri CF, come analisi_eventi.py)
          data/richards4_{NV,V}.json                  (solo per il confronto con la calibrazione del modello a 4 strati)
Scrive  : data/curve_caratteristiche.json              (letto da curvecaratteristiche.html)

Tempi: UTC+1 come nel resto del portale. Suzione s = −ψ (kPa, positiva).

METODO (dettagli e letteratura in curvecaratteristiche.html)
1. Serie orarie (mediana dei dati a 5 min) di θ e di log10 s per ciascuna coppia.
2. Segmentazione in fasi con un filtro ZigZag su θ (mediana mobile 3 h): una fase di WETTING va da un minimo locale al massimo
   successivo, una fase di DRYING dal massimo al minimo successivo; la soglia di inversione (0,005–0,012 m3/m3 secondo il sensore)
   è tarata sul rumore di θ (a 15 cm di NV il ciclo termico diurno estivo arriva a ~0,009).
3. Controllo di coerenza θ–ψ ora per ora (il TEROS 21 risponde con ritardo, specie a suzione alta): un punto di wetting è "valido"
   se θ sale e la suzione scende; un punto di drying è valido se θ non sale e la suzione non scende. Gli altri sono "transitori"
   (ψ non in equilibrio con θ) e non entrano nei fit. Nei fit si escludono anche s < 3 kPa (fuori dal campo di specifica del sensore).
4. Curve: mediane per classi di 0,1 decadi di log s; van Genuchten unimodale (θr, θs, α, n; m = 1 − 1/n) sul ramo di drying;
   ramo di wetting dai punti validi più gli stati di fine wetting (θ al picco, suzione minima raggiunta dopo il picco);
   curva di wetting derivata da quella di drying alla Kool & Parker (θr e n invariati, αw = αd·10^δ con δ ≥ 0, θs,w = κ·θs,d)
   con bootstrap per eventi.
5. Isteresi senza modello: spostamento orizzontale Δlog10 s tra gli stati di fine wetting e la curva di drying a pari θ.
6. Wetting e pioggia: per ogni fase di wetting pioggia nella finestra che precede il picco, intensità, pioggia antecedente,
   ritardi pioggia→θ e θ→ψ; probabilità di risposta per classe di pioggia.
"""
import os, glob, json, math
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy.stats import spearmanr, linregress

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
OUT = os.path.join(ROOT, "data", "curve_caratteristiche.json")

SITES = {"NV": ("SAL_NV", ["015", "030", "060"]), "V": ("SAL_V", ["015", "035", "060"])}
GAUGES = ["CF_21521", "CF_18925", "CF_18957"]
THR = {"NV_015": 0.012, "NV_030": 0.005, "NV_060": 0.006, "V_015": 0.005, "V_035": 0.005, "V_060": 0.005}   # soglia ZigZag su θ
PLATEAU = 0.002        # tolleranza per collocare inizio/fine dell'ascesa sul plateau
TOL_S = 0.003          # decadi/h: velocità di log s oltre la quale si considera una variazione di ψ
TOL_TH = 0.0005        # m3/m3/h
S_FIT_MIN = 3.0        # kPa: i dati con s < 3 kPa non entrano nei fit
BIN_W = 0.1            # decadi
DELAY_H = {"015": 24, "030": 36, "035": 36, "060": 48}      # finestra a monte del wetting per attribuire la pioggia (ore)
NBOOT = 200
RNG = np.random.default_rng(20261120)
S_REF = [10, 30, 100, 300, 1000]      # kPa, per i confronti Δθ a pari suzione


# --------------------------------------------------------------------------- dati
def rawv(sid):
    fs = sorted(glob.glob(os.path.join(ROOT, "data", "raw", sid, "*.csv")))
    d = pd.concat([pd.read_csv(f) for f in fs])
    d["time"] = pd.to_datetime(d["time"])
    return d[d["flag"] < 2].set_index("time")["value"].sort_index()


def hourly_pair(site, dp, idx):
    st = SITES[site][0]
    th, psi = rawv(f"{st}_VWC_{dp}"), rawv(f"{st}_PSI_{dp}")
    d = pd.concat([th.rename("th"), (-psi).rename("s")], axis=1).resample("1h").median().reindex(idx)
    d["s"] = d["s"].clip(lower=0.1)                       # risoluzione del sensore: 0,1 kPa
    d.loc[d["th"].isna() | d["s"].isna()] = np.nan
    return d


def rain_series(idx):
    R = pd.DataFrame({g: rawv(g + "_RAIN") for g in GAUGES}).resample("10min").sum(min_count=1)
    Rm = R.mean(axis=1)                                   # mm per 10 min: dato riferito all'intervallo che termina al timestamp
    return Rm


def rain_events(Rm, gap_h=6):
    wet = Rm.index[(Rm.fillna(0) > 0).values]
    ev, s, last = [], wet[0], wet[0]
    for x in wet[1:]:
        if x - last > pd.Timedelta(hours=gap_h):
            ev.append((s, last)); s = x
        last = x
    ev.append((s, last))
    E = pd.DataFrame(ev, columns=["t_first", "t_last"])
    E["start"] = E.t_first - pd.Timedelta(minutes=10)
    E["end"] = E.t_last
    E["P"] = [Rm[a:b].sum() for a, b in zip(E.t_first, E.t_last)]
    return E


# --------------------------------------------------------------------------- segmentazione
def zigzag(x, thr):
    """Pivot (indice, +1 massimo / −1 minimo) di una serie senza NaN."""
    n = len(x); piv = []
    if n < 2:
        return piv
    mode = 0; lo_i = hi_i = 0; lo = hi = x[0]; ext_i, ext_v = 0, x[0]
    for i in range(1, n):
        v = x[i]
        if mode == 0:
            if v > hi: hi, hi_i = v, i
            if v < lo: lo, lo_i = v, i
            if hi - lo >= thr:
                if hi_i > lo_i: piv.append((lo_i, -1)); mode = 1; ext_i, ext_v = hi_i, hi
                else:           piv.append((hi_i, +1)); mode = -1; ext_i, ext_v = lo_i, lo
        elif mode == 1:
            if v > ext_v: ext_i, ext_v = i, v
            elif ext_v - v >= thr:
                piv.append((ext_i, +1)); mode = -1; ext_i, ext_v = i, v
        else:
            if v < ext_v: ext_i, ext_v = i, v
            elif v - ext_v >= thr:
                piv.append((ext_i, -1)); mode = 1; ext_i, ext_v = i, v
    piv.append((ext_i, +1 if mode == 1 else -1))
    return piv


def refine(x, piv, tol=PLATEAU, frac=0.1):
    """Minimo → ultimo istante del plateau prima della risalita, cioè entro max(tol, frac × ampiezza della risalita) dal minimo
    (il ciclo termico giornaliero o una deriva lenta non devono anticipare l'inizio della risalita);
    massimo → primo istante in cui si raggiunge il plateau. Restituisce (indice raffinato, tipo, indice del pivot originale)."""
    out = []
    for k, (i, kind) in enumerate(piv):
        if kind < 0:
            j = i; lim = piv[k + 1][0] if k + 1 < len(piv) else len(x) - 1
            amp = (x[lim] - x[i]) if k + 1 < len(piv) else 0.0
            tl = max(tol, frac * amp)
            while j + 1 < lim and x[j + 1] <= x[i] + tl: j += 1
        else:
            j = i; lo = piv[k - 1][0] if k > 0 else 0
            while j - 1 > lo and x[j - 1] >= x[i] - tol: j -= 1
        out.append((j, kind, i))
    return out


def prepare(d):
    h = d.copy()
    h["ls"] = np.log10(h["s"])
    h["ths"] = h["th"].interpolate(limit=3, limit_area="inside").rolling(3, center=True, min_periods=2).median()
    h["lss"] = h["ls"].interpolate(limit=3, limit_area="inside").rolling(3, center=True, min_periods=2).median()
    h.loc[h["th"].isna() & h["ths"].isna(), ["ths", "lss"]] = np.nan
    h["rth"] = (h["ths"].shift(-2) - h["ths"].shift(2)) / 4
    h["rs"] = (h["lss"].shift(-2) - h["lss"].shift(2)) / 4
    return h


def segment_blocks(h, thr):
    """Divide in blocchi continui, applica lo ZigZag e restituisce id di fase, tipo (W/D), elenco segmenti."""
    ok = h["ths"].notna().values
    n = len(h)
    seg_id = np.full(n, -1, dtype=int); seg_ty = np.full(n, "", dtype="U1")
    segs = []
    i = 0
    while i < n:
        if not ok[i]: i += 1; continue
        j = i
        while j + 1 < n and ok[j + 1]: j += 1
        if j - i >= 12:
            x = h["ths"].values[i:j + 1]
            piv = refine(x, zigzag(x, thr))
            b = [0] + [p[0] for p in piv] + [len(x) - 1]
            kinds = [None] + [p[1] for p in piv] + [None]
            iorig = [None] + [p[2] for p in piv] + [None]
            for m in range(len(b) - 1):
                a_, c_ = b[m], b[m + 1]
                if c_ <= a_: continue
                if m == 0:   ty = "D" if (piv and piv[0][1] < 0) else "W"
                elif m == len(b) - 2: ty = "W" if kinds[m] < 0 else "D"
                else: ty = "W" if kinds[m + 1] > 0 else "D"
                partial = (m == 0) or (m == len(b) - 2)
                sid = len(segs)
                iv = i + (iorig[m] if (ty == "W" and m > 0) else a_)      # minimo vero di θ (prima della correzione del plateau)
                th0 = float(x[iorig[m]]) if (ty == "W" and m > 0) else float(x[a_])
                segs.append(dict(id=sid, ty=ty, a=i + a_, b=i + c_, partial=partial, iv=int(iv),
                                 th0=th0, th1=float(x[c_])))
        i = j + 1
    for s in segs:      # l'ora di pivot appartiene alla fase che inizia
        seg_id[s["a"]:s["b"] + 1] = s["id"]; seg_ty[s["a"]:s["b"] + 1] = s["ty"]
    return seg_id, seg_ty, segs


def classify(h, ty):
    """W valido, D valido, T transitorio (ψ non in equilibrio con θ), '-' non utilizzabile."""
    rs, rth = h["rs"].values, h["rth"].values
    ph = np.full(len(h), "-", dtype="U1")
    okn = np.isfinite(rs) & np.isfinite(rth) & (ty != "")
    w = okn & (ty == "W"); d = okn & (ty == "D")
    ph[w] = np.where((rs[w] < -TOL_S) & (rth[w] > 0), "W", "T")
    ph[d] = np.where((rs[d] > -TOL_S) & (rth[d] <= TOL_TH), "D", "T")
    return ph


def reversal_points(h, segs, post_max=14 * 24):
    """Stati di inversione di ogni wetting: minimo (θ_v, s_v) e massimo (θ_p al picco, suzione minima s_p dopo il picco)."""
    out = []
    ls = h["lss"].values; ths = h["ths"].values; n = len(h)
    for k, s in enumerate(segs):
        if s["ty"] != "W": continue
        i0, i1 = s["a"], s["b"]
        nxt = n - 1
        for s2 in segs[k + 1:]:
            if s2["ty"] == "W" and s2["a"] > i1: nxt = s2["a"]; break
        j_end = min(i1 + post_max, nxt, n - 1)
        # il segmento D che segue può finire in un buco di dati: ci si ferma al primo NaN
        seg_slice = ls[max(0, i1 - 6):j_end + 1]
        if len(seg_slice) < 3 or not np.isfinite(seg_slice).any(): continue
        jm = int(np.nanargmin(seg_slice)) + max(0, i1 - 6)
        out.append(dict(k=k, id=s["id"], i_v=i0, i_p=i1, i_smin=jm, j_end=j_end,
                        th_v=float(np.nanmedian(ths[max(0, s["iv"] - 1):s["iv"] + 2])), ls_v=float(np.nanmedian(ls[max(0, i0 - 1):i0 + 2])),
                        th_p=float(s["th1"]), ls_p=float(ls[jm]), th_at_smin=float(ths[jm]),
                        interior=bool(jm < j_end - 2)))
    return out


# --------------------------------------------------------------------------- curve
def vg(ls, p):
    thr, ths, la, n = p
    return thr + (ths - thr) * (1 + (10 ** la * 10 ** np.asarray(ls)) ** n) ** (-(1 - 1 / n))


def vg_ls_of_theta(th, p):
    thr, ths, la, n = p; m = 1 - 1 / n
    Se = np.clip((np.asarray(th) - thr) / (ths - thr), 1e-6, 1 - 1e-9)
    return np.log10(Se ** (-1 / m) - 1) / n - la


def bins(ls, th, width=BIN_W, minn=3):
    ls = np.asarray(ls); th = np.asarray(th)
    if len(ls) == 0: return np.zeros((0, 6))
    b = np.floor((ls - (-1.0)) / width).astype(int)
    out = []
    for k in np.unique(b):
        m = b == k
        if m.sum() >= minn:
            out.append((ls[m].mean(), np.median(th[m]), np.quantile(th[m], .1), np.quantile(th[m], .9), m.sum(), th[m].std()))
    return np.array(out)


def _bounds_vg(y, th_all):
    lo = [0.0, max(float(np.max(y)), float(np.max(th_all)) * 0.98), -3.5, 1.05]
    hi = [max(1e-6, min(float(np.min(y)), float(np.min(th_all)))), 0.7, 1.5, 8.0]
    return lo, hi


def fit_vg(ls, th, width=BIN_W, minn=3):
    B = bins(ls, th, width, minn)
    if len(B) < 6: return None
    x, y = B[:, 0], B[:, 1]
    lo, hi = _bounds_vg(y, th)
    best = None
    for thr0 in (0.0, hi[0] * 0.5):
        for ths0 in (lo[1] * 1.02, lo[1] * 1.3):
            for la0 in (-2.0, -1.0, 0.0):
                for n0 in (1.3, 2.0, 4.0):
                    p0 = np.clip([thr0, ths0, la0, n0], lo, hi)
                    try: r = least_squares(lambda p: vg(x, p) - y, p0, bounds=(lo, hi))
                    except Exception: continue
                    if best is None or r.cost < best.cost: best = r
    p = best.x
    return dict(p=[float(v) for v in p], rmse=float(np.sqrt(np.mean((vg(x, p) - y) ** 2))), nbins=int(len(B)),
                npts=int(len(ls)), at_bound=[bool(abs(p[i] - lo[i]) < 1e-9 or abs(p[i] - hi[i]) < 1e-9) for i in range(4)])


def wet_fit(pD, BW):
    """Curva di wetting derivata da quella di drying (Kool e Parker): θr e n invariati, αw = αd·10^δ con δ ∈ [0; 1,2]
    (in equilibrio il wetting non può stare a destra del drying), θs,w = κ·θs,d con κ ∈ [0,5; 1]. Restituisce [δ, κ]."""
    thr, ths, lad, n = pD
    xw, yw = BW[:, 0], BW[:, 1]
    lo, hi = [0.0, 0.5], [1.2, 1.0]
    res = lambda q: vg(xw, [thr, q[1] * ths, lad + q[0], n]) - yw
    best = None
    for d0 in (0.0, 0.3, 0.7):
        for k0 in (0.7, 0.9, 1.0):
            try: r = least_squares(res, [d0, k0], bounds=(lo, hi))
            except Exception: continue
            if best is None or r.cost < best.cost: best = r
    return best.x


def monotone_curve(B):
    """Curva θ(log s) mediana, resa monotona decrescente: restituisce (ls, th)."""
    ls, th = B[:, 0], B[:, 1]
    o = np.argsort(ls); ls, th = ls[o], th[o]
    th = np.minimum.accumulate(th)
    return ls, th


def shift_to_curve(ls_c, th_c, th_pts, ls_pts):
    """Δlog10 s = log s della curva di drying al θ del punto − log s del punto (positivo: il punto è più umido/a minor suzione)."""
    # th_c decrescente → invertire per np.interp
    tc = th_c[::-1]; lc = ls_c[::-1]
    # tolleranza: elimina valori uguali
    ok = np.concatenate([[True], np.diff(tc) > 1e-9])
    tc, lc = tc[ok], lc[ok]
    inside = (th_pts > tc.min()) & (th_pts < tc.max())
    return np.where(inside, np.interp(th_pts, tc, lc) - ls_pts, np.nan), inside


def boot_ci(vals, B=1000, stat=np.median, q=(2.5, 97.5)):
    vals = np.asarray([v for v in vals if np.isfinite(v)])
    if len(vals) < 3: return [None, None]
    bs = [stat(RNG.choice(vals, len(vals), replace=True)) for _ in range(B)]
    return [float(np.percentile(bs, q[0])), float(np.percentile(bs, q[1]))]


def interp_th_at(B, s_ref):
    ls, th = monotone_curve(B)
    if len(ls) < 3: return None
    if np.log10(s_ref) < ls.min() or np.log10(s_ref) > ls.max(): return None
    return float(np.interp(np.log10(s_ref), ls, th))


# --------------------------------------------------------------------------- principale
def main():
    idx = pd.date_range("2025-01-31 12:00", "2026-09-30 23:00", freq="1h")
    Rm = rain_series(idx); Rcum = Rm.fillna(0).cumsum()
    rain_in = lambda a, b: float(Rcum.asof(b) - Rcum.asof(a))
    E = rain_events(Rm)
    Rh = Rm.fillna(0).resample("1h", closed="right", label="left").sum().reindex(idx).fillna(0)

    rich = {k: json.load(open(os.path.join(ROOT, "data", f"richards4_{k}.json"))) for k in SITES}
    KPA_M = 1 / 9.81

    pairs, events_all, resp_rows = {}, [], []
    for site, (st, dps) in SITES.items():
        for li, dp in enumerate(dps):
            key = f"{site}_{dp}"
            d = hourly_pair(site, dp, idx)
            h = prepare(d)
            seg_id, seg_ty, segs = segment_blocks(h, THR[key])
            ph = classify(h, seg_ty)
            rev = reversal_points(h, segs)
            th, ls = h["th"].values, h["ls"].values
            thr_ = THR[key]

            # dati per le curve
            fit_ok = np.isfinite(h["ths"].values) & np.isfinite(h["lss"].values) & (10 ** h["lss"].values >= S_FIT_MIN)
            mD = (ph == "D") & fit_ok; mW = (ph == "W") & fit_ok
            ix = np.arange(len(h)); year = h.index.year.values
            d25 = mD & (h.index < "2025-10-01"); d26 = mD & (h.index >= "2026-01-01")
            lsS, thS = h["lss"].values, h["ths"].values
            BD = bins(lsS[mD], thS[mD], BIN_W, 3)
            BD25 = bins(lsS[d25], thS[d25], BIN_W, 3); BD26 = bins(lsS[d26], thS[d26], BIN_W, 3)

            # eventi di wetting con stati di fine wetting
            segW = {s["id"]: s for s in segs if s["ty"] == "W"}
            thD = h["ths"].values
            evrows = []
            prev_peak_t = None
            for r in rev:
                s = segs[r["k"]]
                t0, t1 = h.index[s["a"]], h.index[s["b"]]
                ev = dict(seg=s["id"], partial=bool(s["partial"]), i0=int(s["a"]), i1=int(s["b"]),
                          th0=r["th_v"], th1=r["th_p"], ls0=r["ls_v"], ls_end=r["ls_p"], th_end=r["th_at_smin"],
                          interior=r["interior"], i_smin=int(r["i_smin"]))
                ev["dth"] = ev["th1"] - ev["th0"]; ev["dls"] = ev["ls0"] - ev["ls_end"]
                ev["dur_h"] = float((t1 - t0).total_seconds() / 3600)
                # lag di risposta di ψ rispetto a θ: primo istante in cui log s scende di 0,05 sotto il valore iniziale
                lw = h["lss"].values[max(0, s["a"] - 6):r["j_end"] + 1]
                below = np.where(lw <= ev["ls0"] - 0.05)[0]
                ev["lag_on_h"] = float(below[0] - (s["a"] - max(0, s["a"] - 6))) if (len(below) and 10 ** ev["ls0"] >= S_FIT_MIN) else None
                # (a suzioni iniziali < 3 kPa, cioè vicino alla saturazione, il calo di 0,05 decadi non è misurabile con il TEROS 21)
                # pioggia
                ws = t0 - pd.Timedelta(hours=DELAY_H[dp])
                if prev_peak_t is not None: ws = max(ws, prev_peak_t)
                P = rain_in(ws, t1)
                seg_r = Rm[ws:t1]
                ev["P"] = P
                ev["I10"] = float(seg_r.max() * 6) if len(seg_r) and np.isfinite(seg_r.max()) else None
                ev["I60"] = float(seg_r.rolling(6, min_periods=1).sum().max()) if len(seg_r) else None
                first_wet = seg_r[seg_r > 0]
                ev["rain_to_theta_h"] = float((t0 - first_wet.index[0]).total_seconds() / 3600) if len(first_wet) else None
                ev["rain_dur_h"] = float((first_wet.index[-1] - first_wet.index[0]).total_seconds() / 3600) if len(first_wet) else None
                ev["A7"] = rain_in(ws - pd.Timedelta(days=7), ws)
                ev["A3"] = rain_in(ws - pd.Timedelta(days=3), ws)
                prev_peak_t = t1
                evrows.append(ev)
                s["_ev"] = ev
            # intervallo di θ del sensore e tipo (principale/scansione)
            th_lo, th_hi = float(np.nanpercentile(th, 1)), float(np.nanpercentile(th, 99))
            for ev in evrows:
                ev["se0"] = (ev["th0"] - th_lo) / (th_hi - th_lo); ev["se1"] = (ev["th1"] - th_lo) / (th_hi - th_lo)
                ev["kind"] = "principale" if (ev["se0"] < 0.2 and ev["se1"] > 0.7) else "scansione"

            # punti di fine wetting utilizzabili e dati W
            peaks = [e for e in evrows if (not e["partial"]) and e["interior"] and e["dls"] > 0.1 and 10 ** e["ls_end"] >= S_FIT_MIN]
            # W: punti validi con id evento; D: punti validi con id fase
            Wl = [(ix_, seg_id[ix_]) for ix_ in np.where(mW)[0]]
            thW = np.concatenate([thS[mW], [e["th1"] for e in peaks]]); lsW = np.concatenate([lsS[mW], [e["ls_end"] for e in peaks]])
            BW = bins(lsW, thW, 0.2, 2)

            fits = {}
            fD = fit_vg(lsS[mD], thS[mD]) if mD.sum() > 200 else None
            fD25 = fit_vg(lsS[d25], thS[d25]) if d25.sum() > 200 else None
            fD26 = fit_vg(lsS[d26], thS[d26]) if d26.sum() > 200 else None
            fW = fit_vg(lsW, thW, 0.2, 2) if len(BW) >= 6 else None
            fits["D"], fits["D2025"], fits["D2026"], fits["W_free"] = fD, fD25, fD26, fW

            # curva di wetting dal drying + bootstrap per eventi
            joint = None
            if fD and len(BW) >= 4 and len(BD) >= 8:
                pD = fD["p"]
                q = wet_fit(pD, BW)
                wsegs = sorted(set(seg_id[mW]) | {e["seg"] for e in peaks})
                idxW = {sid: np.where(mW & (seg_id == sid))[0] for sid in wsegs}
                pk = {e["seg"]: e for e in peaks}
                boots = []
                if len(wsegs) >= 3:
                    for _ in range(NBOOT):
                        sw = RNG.choice(wsegs, len(wsegs))
                        iW = np.concatenate([idxW[s_] for s_ in sw])
                        lw_ = np.concatenate([lsS[iW], [pk[s_]["ls_end"] for s_ in sw if s_ in pk]])
                        tw_ = np.concatenate([thS[iW], [pk[s_]["th1"] for s_ in sw if s_ in pk]])
                        bW_ = bins(lw_, tw_, 0.2, 2)
                        if len(bW_) < 4: continue
                        try: boots.append(wet_fit(pD, bW_))
                        except Exception: continue
                boots = np.array(boots) if len(boots) else None
                def ci(v): return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] if v is not None and len(v) >= 20 else [None, None]
                ratio = 10 ** q[0]
                joint = dict(thr=float(pD[0]), ths_d=float(pD[1]), alpha_d=float(10 ** pD[2]), n=float(pD[3]), alpha_w=float(10 ** (pD[2] + q[0])), kappa=float(q[1]),
                             ratio_alpha=float(ratio), ths_w=float(q[1] * pD[1]),
                             ratio_alpha_ci=ci(10 ** boots[:, 0]) if boots is not None else [None, None],
                             ratio_alpha_p95=float(np.percentile(10 ** boots[:, 0], 95)) if boots is not None and len(boots) >= 20 else None,
                             kappa_ci=ci(boots[:, 1]) if boots is not None else [None, None],
                             n_boot=int(0 if boots is None else len(boots)),
                             rmse_d=float(fD["rmse"]),
                             rmse_w=float(np.sqrt(np.mean((vg(BW[:, 0], [pD[0], q[1] * pD[1], pD[2] + q[0], pD[3]]) - BW[:, 1]) ** 2))),
                             at_bound=dict(alpha_w=bool(q[0] <= 1e-6 or q[0] >= 1.199), kappa=bool(q[1] >= 0.999 or q[1] <= 0.501)))
            fits["joint"] = joint

            # spostamento orizzontale (senza modello) rispetto alla curva di drying
            shifts = dict(peaks=None, wvalid=None)
            if len(BD) >= 8:
                lc, tc = monotone_curve(BD)
                if peaks:
                    sh, ins = shift_to_curve(lc, tc, np.array([e["th1"] for e in peaks]), np.array([e["ls_end"] for e in peaks]))
                    for e, v in zip(peaks, sh): e["shift"] = None if not np.isfinite(v) else float(v)
                    shifts["peaks"] = dict(n=int(np.isfinite(sh).sum()), median=float(np.nanmedian(sh)) if np.isfinite(sh).any() else None,
                                           q25=float(np.nanpercentile(sh, 25)) if np.isfinite(sh).any() else None,
                                           q75=float(np.nanpercentile(sh, 75)) if np.isfinite(sh).any() else None,
                                           ci=boot_ci(sh[np.isfinite(sh)]))
                # wetting valido: mediana per evento, poi mediana tra eventi
                per_ev = []
                for sid in sorted(set(seg_id[mW])):
                    m = mW & (seg_id == sid)
                    sh, ins = shift_to_curve(lc, tc, thS[m], lsS[m])
                    if np.isfinite(sh).sum() >= 3: per_ev.append(float(np.nanmedian(sh)))
                if len(per_ev) >= 3:
                    shifts["wvalid"] = dict(n=len(per_ev), median=float(np.median(per_ev)), q25=float(np.percentile(per_ev, 25)),
                                            q75=float(np.percentile(per_ev, 75)), ci=boot_ci(per_ev))
            # deviazione orizzontale dalla curva di drying (decadi di log s; positiva = a sinistra, cioè a suzione minore: lato isteresi;
            # negativa = a destra, cioè ψ in ritardo rispetto a θ: lato dinamico)
            devstats = None
            if len(BD) >= 8:
                lc, tc = monotone_curve(BD)
                tcr, lcr = tc[::-1], lc[::-1]
                okk = np.concatenate([[True], np.diff(tcr) > 1e-9]); tcr, lcr = tcr[okk], lcr[okk]
                def dev_of(thv, lsv):
                    return np.where((thv > tcr.min()) & (thv < tcr.max()), np.interp(thv, tcr, lcr) - lsv, np.nan)
                dev_all = dev_of(thS, lsS)
                inW = (seg_ty == "W") & fit_ok
                def st(m):
                    v = dev_all[m]; v = v[np.isfinite(v)]
                    if len(v) < 10: return None
                    return dict(n=int(len(v)), median=float(np.median(v)), left=float(np.mean(v > 0.1)), right=float(np.mean(v < -0.1)))
                devstats = dict(W_valid=st(mW), W_transient=st(inW & (ph == "T")), D_valid=st(mD))
                lag_sens = []
                for tau in (0, 3, 6, 12):
                    lsc = pd.Series(h["lss"].values).shift(-tau).values
                    rsc = (pd.Series(lsc).shift(-2).values - pd.Series(lsc).shift(2).values) / 4
                    m = (seg_ty == "W") & (rsc < -TOL_S) & (h["rth"].values > 0) & np.isfinite(lsc) & (10 ** np.where(np.isfinite(lsc), lsc, 0) >= S_FIT_MIN)
                    v = dev_of(thS, lsc)[m]; v = v[np.isfinite(v)]
                    lag_sens.append(dict(tau=tau, n=int(len(v)), median=float(np.median(v)) if len(v) else None,
                                         left=float(np.mean(v > 0.1)) if len(v) else None))
                devstats["lag_sens"] = lag_sens
                # per evento: deviazione massima a destra (dinamica) e a sinistra (isteresi) e velocità massima di risalita di θ
                for ev in evrows:
                    a_, b_ = ev["i0"], max(ev["i_smin"], ev["i1"])
                    v = dev_all[a_:b_ + 1]; v = v[np.isfinite(v)]
                    ev["dyn_dev"] = float(max(0.0, -np.min(v))) if len(v) else None
                    ev["hys_dev"] = float(max(0.0, np.max(v))) if len(v) else None
                    rr = h["rth"].values[ev["i0"]:ev["i1"] + 1]; ev["rate_max"] = float(np.nanmax(rr)) if np.isfinite(rr).any() else None
            # Δθ a pari suzione: drying vs wetting (classi) e tra cicli
            dth_ref = {}
            pe_ = None if not fD else [fD["p"][0], 0.85 * fD["p"][1], fD["p"][2] + math.log10(2.0), fD["p"][3]]
            for sr in S_REF:
                a = interp_th_at(BD, sr)
                a25, a26 = interp_th_at(BD25, sr) if len(BD25) > 4 else None, interp_th_at(BD26, sr) if len(BD26) > 4 else None
                w = interp_th_at(BW, sr) if len(BW) > 4 else None
                dth_ref[str(sr)] = dict(D=a, D2025=a25, D2026=a26, W=w,
                                        E=None if pe_ is None else float(vg(math.log10(sr), pe_)),
                                        Dfit=None if not fD else float(vg(math.log10(sr), fD["p"])),
                                        Wj=None if not joint else float(vg(math.log10(sr), [joint["thr"], joint["ths_w"], math.log10(joint["alpha_w"]), joint["n"]])))

            # confronto con la calibrazione Richards 4 strati (strato corrispondente)
            rp = rich[site]["params"][li]
            rich_curve = dict(thr=rp["tr"], ths=rp["ts"], alpha=(10 ** rp["la"]) * KPA_M, n=rp["n"])

            # curve da tracciare
            lgrid = np.round(np.arange(-0.5, 3.61, 0.05), 2)
            curves = dict(ls=lgrid.tolist())
            if fD: curves["D"] = [round(float(v), 4) for v in vg(lgrid, fD["p"])]
            if fD25: curves["D2025"] = [round(float(v), 4) for v in vg(lgrid, fD25["p"])]
            if fD26: curves["D2026"] = [round(float(v), 4) for v in vg(lgrid, fD26["p"])]
            if joint:
                curves["W"] = [round(float(v), 4) for v in vg(lgrid, [joint["thr"], joint["ths_w"], math.log10(joint["alpha_w"]), joint["n"]])]
            if fD:
                pe = [fD["p"][0], 0.85 * fD["p"][1], fD["p"][2] + math.log10(2.0), fD["p"][3]]     # Kool & Parker (1987)
                curves["E"] = [round(float(v), 4) for v in vg(lgrid, pe)]
            curves["R"] = [round(float(v), 4) for v in vg(lgrid, [rich_curve["thr"], rich_curve["ths"], math.log10(rich_curve["alpha"]), rich_curve["n"]])]
            rng_ls = [float(np.nanmin(lsS[fit_ok])), float(np.nanmax(lsS[fit_ok]))]

            # drying: statistiche sulle fasi
            dstats = []
            for s in segs:
                if s["ty"] != "D" or s["partial"]: continue
                dur = (s["b"] - s["a"]) / 24
                dstats.append(dict(id=s["id"], i0=int(s["a"]), i1=int(s["b"]), th0=s["th0"], th1=s["th1"], dur_d=float(dur),
                                   rate=float((s["th1"] - s["th0"]) / max(dur, 1e-6))))

            # probabilità di risposta alla pioggia (episodi separati da ≥ 24 h)
            ep = E.copy().sort_values("start").reset_index(drop=True)
            cl, cur = [], None
            for _, e_ in ep.iterrows():
                if cur is None or e_.start - cur["end"] > pd.Timedelta(hours=24):
                    if cur is not None: cl.append(cur)
                    cur = dict(start=e_.start, end=e_.end, P=e_.P)
                else:
                    cur["end"] = e_.end; cur["P"] += e_.P
            cl.append(cur)
            for c in cl:
                if c["P"] < 2: continue
                a = c["start"]; b = c["end"] + pd.Timedelta(hours=DELAY_H[dp])
                if a < h.index[0] or b > h.index[-1]: continue
                base = h["ths"].asof(a)
                w = h["ths"].loc[a:b]
                if not np.isfinite(base) or w.notna().sum() < 6: continue
                resp = float(np.nanmax(w.values) - base)
                resp_rows.append(dict(key=key, P=float(c["P"]), rise=resp, resp=bool(resp >= thr_)))

            # serializzazione compatta
            def pack(a, scale, nd):
                return [None if not np.isfinite(v) else int(round(v * scale)) for v in a]
            pairs[key] = dict(
                site=site, depth=int(dp), label=f"{site} · {int(dp)} cm", thr=thr_,
                th=pack(th, 1000, 3), ls=pack(ls, 100, 2), ph="".join(ph.tolist()),
                n=int(len(h)), th_range=[th_lo, th_hi], ls_range=rng_ls,
                counts=dict(W=int((ph == "W").sum()), D=int((ph == "D").sum()), T=int((ph == "T").sum()),
                            valid_fit_W=int(mW.sum()), valid_fit_D=int(mD.sum())),
                segments=[dict(id=s["id"], ty=s["ty"], a=int(s["a"]), b=int(s["b"]), partial=bool(s["partial"]),
                               th0=round(s["th0"], 4), th1=round(s["th1"], 4)) for s in segs],
                events=[{k: (round(v, 4) if isinstance(v, float) else v) for k, v in e.items()} for e in evrows],
                dry=dstats,
                fits={k: v for k, v in fits.items()},
                shifts=shifts, dth_ref=dth_ref, curves=curves, devstats=devstats,
                richards=rich_curve,
                bins_D=[[round(float(v), 4) for v in row[:5]] for row in BD],
                bins_W=[[round(float(v), 4) for v in row[:5]] for row in BW],
            )
            print(f"{key}: seg W={sum(1 for s in segs if s['ty']=='W')} D={sum(1 for s in segs if s['ty']=='D')}  "
                  f"W/D/T h = {(ph=='W').sum()}/{(ph=='D').sum()}/{(ph=='T').sum()}  fit D rmse={fD['rmse'] if fD else None}  "
                  f"joint ratio={None if not joint else round(joint['ratio_alpha'],2)} ci={None if not joint else joint['ratio_alpha_ci']}  "
                  f"shift peaks median={None if not shifts['peaks'] else shifts['peaks']['median']}")
            for e in evrows:
                e["key"] = key
            events_all += evrows

    # correlazioni wetting–pioggia
    ev = pd.DataFrame([{**{k: v for k, v in e.items() if not isinstance(v, (list, dict))}} for e in events_all])
    corr = {}
    for key, g in ev.groupby("key"):
        g = g[(~g.partial) & (g.P >= 1)]
        c = dict(n=int(len(g)))
        for name, col in (("P", "P"), ("I60", "I60"), ("I10", "I10"), ("A7", "A7")):
            gg = g.dropna(subset=[col])
            if len(gg) >= 5:
                r = spearmanr(gg.dth, gg[col]); c["rho_" + name] = float(r[0]); c["p_" + name] = float(r[1])
        gg = g[g.P > 0]
        if len(gg) >= 5:
            lr = linregress(np.log10(gg.P), np.log10(gg.dth.clip(lower=1e-4)))
            c["loglog"] = dict(slope=float(lr.slope), intercept=float(lr.intercept), r2=float(lr.rvalue ** 2), p=float(lr.pvalue))
        c["delay_rain_theta_med"] = float(np.nanmedian(g.rain_to_theta_h.astype(float))) if g.rain_to_theta_h.notna().any() else None
        c["delay_rain_theta_q25"] = float(np.nanpercentile(g.rain_to_theta_h.dropna(), 25)) if g.rain_to_theta_h.notna().any() else None
        c["delay_rain_theta_q75"] = float(np.nanpercentile(g.rain_to_theta_h.dropna(), 75)) if g.rain_to_theta_h.notna().any() else None
        c["lag_psi_med"] = float(np.nanmedian(g.lag_on_h.astype(float))) if g.lag_on_h.notna().any() else None
        c["no_rain"] = int((ev[(ev.key == key) & (~ev.partial)].P < 1).sum())
        gd = g.dropna(subset=["dyn_dev", "rate_max"]) if "dyn_dev" in g else g.iloc[0:0]
        if len(gd) >= 5:
            r = spearmanr(gd.rate_max, gd.dyn_dev); c["rho_dyn_rate"] = float(r[0]); c["p_dyn_rate"] = float(r[1]); c["n_dyn"] = int(len(gd))
        gl = g.dropna(subset=["lag_on_h", "ls0"])
        if len(gl) >= 5:
            r = spearmanr(gl.ls0, gl.lag_on_h.astype(float)); c["rho_lag_s0"] = float(r[0]); c["p_lag_s0"] = float(r[1]); c["n_lag"] = int(len(gl))
        corr[key] = c

    resp = pd.DataFrame(resp_rows)
    classes = [(2, 5), (5, 10), (10, 20), (20, 40), (40, 1000)]
    resp_tab = {}
    for key, g in resp.groupby("key"):
        rows = []
        for a, b in classes:
            m = (g.P >= a) & (g.P < b)
            rows.append(dict(lo=a, hi=b if b < 1000 else None, n=int(m.sum()), resp=int(g.resp[m].sum()),
                             rise_med=float(g.rise[m].median()) if m.any() else None))
        # P50 da regressione logistica su log10 P
        P50 = None
        if g.resp.nunique() == 2 and len(g) >= 15:
            x = np.log10(g.P.values); yv = g.resp.values.astype(float)
            def nll(p):
                z = p[0] + p[1] * x; return np.sum(np.log1p(np.exp(-np.abs(z))) + np.maximum(z, 0) - yv * z)
            from scipy.optimize import minimize
            r = minimize(nll, [0.0, 1.0], method="Nelder-Mead")
            if r.x[1] > 0.05: P50 = float(10 ** (-r.x[0] / r.x[1]))
        resp_tab[key] = dict(classes=rows, P50=P50, n=int(len(g)), resp=int(g.resp.sum()))

    out = dict(
        generated=pd.Timestamp.now(tz="Europe/Rome").isoformat(timespec="seconds"),
        t0=idx[0].strftime("%Y-%m-%dT%H:%M"), n=int(len(idx)),
        params=dict(thr=THR, plateau=PLATEAU, tol_s=TOL_S, tol_th=TOL_TH, s_fit_min=S_FIT_MIN, bin_w=BIN_W, delay_h=DELAY_H,
                    nboot=NBOOT, s_ref=S_REF, rain_gap_h=6),
        rain=dict(mm_h=[int(round(v * 10)) for v in Rh.values]),
        pairs=pairs, correlations=corr, response=resp_tab,
        rain_events=[dict(t=e.start.strftime("%Y-%m-%dT%H:%M"), te=e.end.strftime("%Y-%m-%dT%H:%M"), P=round(float(e.P), 1)) for _, e in E[E.P >= 5].iterrows()],
        rain_total=dict(mm=round(float(Rm.sum()), 1)),
    )
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, separators=(",", ":"), default=lambda o: None if (isinstance(o, float) and not np.isfinite(o)) else (float(o) if isinstance(o, (np.floating,)) else (int(o) if isinstance(o, np.integer) else (bool(o) if isinstance(o, np.bool_) else str(o)))))
    print("scritto", OUT, round(os.path.getsize(OUT) / 1e6, 2), "MB")


if __name__ == "__main__":
    main()
