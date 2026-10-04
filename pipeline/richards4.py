"""
richards4.py — modello di Richards 1D a 4 strati per un pendio indefinito (siti NV e V).

Strati (profondità verticali): NV 0–22.5 / 22.5–45 / 45–75 / 75–100 cm ; V 0–25 / 25–47.5 / 47.5–75 / 75–100 cm.
I primi tre sono monitorati (VWC e potenziale matriciale a 15, 30|35, 60 cm), il quarto no.
Il flusso è 1D normale al pendio: coordinata n = z·cosβ, gravità cosβ, flussi per unità di area di pendio (pioggia·cosβ).
Condizione al fondo: λ·K(h)·cosβ (λ=0 roccia impermeabile; λ=1 drenaggio libero) — λ è calibrato e assorbe anche
il drenaggio laterale/fratture. Condizione in alto: pioggia (media 3 pluviometri CF); se la superficie satura il surplus è ruscellamento.
Evapotraspirazione: PET Hargreaves-Samani da ERA5 Tair × Kc, prelevata nei primi strati con profilo esponenziale e riduzione di Feddes.
Tempi UTC+1.
"""
import os, glob, math
import numpy as np, pandas as pd
from numba import njit

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
BETA = math.radians(30.0)
COSB = math.cos(BETA)
DZ = 0.025                       # m (verticale)
NCELL = 40                       # 1 m
KPA_M = 1.0 / 9.81               # 1 kPa = 0.1019 m di colonna d'acqua
SITES = {
    "NV": dict(st="SAL_NV", dep=["015", "030", "060"], z=[0.15, 0.30, 0.60], cuts=[9, 18, 30], zr=0.05),
    "V":  dict(st="SAL_V",  dep=["015", "035", "060"], z=[0.15, 0.35, 0.60], cuts=[10, 19, 30], zr=0.25),
}
GAUGES = ["CF_21521", "CF_18925", "CF_18957"]


def rawv(sid):
    fs = sorted(glob.glob(os.path.join(ROOT, "data", "raw", sid, "*.csv")))
    d = pd.concat([pd.read_csv(f) for f in fs]); d["time"] = pd.to_datetime(d["time"])
    return d[d["flag"] < 2].set_index("time")["value"].sort_index()


def pet_daily():
    t = rawv("ERA5_E147_N407_TAIR").loc["2024-12-25":"2025-12-31"]
    tmax, tmin, tm = t.resample("1D").max(), t.resample("1D").min(), t.resample("1D").mean()
    lat = math.radians(40.69)
    def Ra(doy):
        dr = 1 + 0.033 * math.cos(2 * math.pi * doy / 365); d = 0.409 * math.sin(2 * math.pi * doy / 365 - 1.39)
        ws = math.acos(-math.tan(lat) * math.tan(d))
        return 24 * 60 / math.pi * 0.0820 * dr * (ws * math.sin(lat) * math.sin(d) + math.cos(lat) * math.cos(d) * math.sin(ws))
    ra = pd.Series([Ra(d.dayofyear) for d in tm.index], index=tm.index)
    return (0.0023 * (tm + 17.8) * np.sqrt((tmax - tmin).clip(lower=0)) * 0.408 * ra).loc["2025"]   # mm/g


def forcing(dt_min=30):
    """pioggia (mm per passo) e PET (mm per passo) su griglia regolare dal 2025-01-01 00:00 al 2025-12-31 23:30"""
    R = pd.DataFrame({g: rawv(g + "_RAIN") for g in GAUGES}).loc["2025"].resample("10min").sum(min_count=1)
    r = R.mean(axis=1).fillna(0.0).resample(f"{dt_min}min", closed="right", label="right").sum()
    r = r.loc["2025-01-01 00:00":"2025-12-31 23:30"]
    pet = pet_daily().reindex(r.index.floor("1D")).values / (24 * 60 / dt_min)
    return r.index, r.values.astype(np.float64), pet.astype(np.float64)


@njit(cache=True, error_model='numpy')
def _vg(h, ts, tr, a, n, ks):
    m = 1.0 - 1.0 / n
    if h >= 0.0:
        return ts + 1e-4 * h, ks, 1e-4, 1.0
    x = (a * -h) ** n
    se = (1.0 + x) ** (-m)
    th = tr + (ts - tr) * se
    t1 = se ** (1.0 / m)
    k = ks * math.sqrt(se) * (1.0 - (1.0 - t1) ** m) ** 2
    c = a * m * n / (1.0 - m) * (ts - tr) * t1 * (1.0 - t1) ** m
    if k < 1e-18: k = 1e-18
    if c < 1e-6: c = 1e-6
    return th, k, c, se


@njit(cache=True, error_model='numpy')
def _thomas(a, b, c, d, n):
    cp = np.empty(n); dp = np.empty(n)
    cp[0] = c[0] / b[0]; dp[0] = d[0] / b[0]
    for i in range(1, n):
        den = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / den; dp[i] = (d[i] - a[i] * dp[i - 1]) / den
    x = np.empty(n); x[n - 1] = dp[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


@njit(cache=True, error_model='numpy')
def _run(h0, ts, tr, al, nn, ks, lam, kc, w, rain, pet, dt_s, cosb, dz, hwp, ha):
    N = h0.shape[0]; dn = dz * cosb
    h = h0.copy()
    nst = rain.shape[0]
    TH = np.empty((nst, N)); HH = np.empty((nst, N)); RUN = np.zeros(nst); AET = np.zeros(nst)
    thn = np.empty(N)
    for i in range(N):
        thn[i] = _vg(h[i], ts[i], tr[i], al[i], nn[i], ks[i])[0]
    for s in range(nst):
        tleft = dt_s; dts = dt_s
        qtop = rain[s] * 1e-3 / dt_s * cosb          # m/s (per area di pendio)
        run_acc = 0.0; aet_acc = 0.0
        while tleft > 1e-6:
            if dts > tleft: dts = tleft
            ok = False
            for attempt in range(12):
                hk = h.copy(); conv = False
                dirich = False
                for rep in range(2):
                    for it in range(25):
                        th = np.empty(N); K = np.empty(N); C = np.empty(N)
                        for i in range(N):
                            th[i], K[i], C[i], _ = _vg(hk[i], ts[i], tr[i], al[i], nn[i], ks[i])
                        A = np.zeros(N); B = np.zeros(N); Cc = np.zeros(N); D = np.zeros(N)
                        # sink evapotraspirazione (m/s per cella, volume/area di pendio -> flusso)
                        sink = np.zeros(N)
                        pe = pet[s] * 1e-3 / dt_s * cosb * kc
                        for i in range(N):
                            if w[i] > 0.0:
                                f = 1.0 if hk[i] > ha else (0.0 if hk[i] < hwp else (hk[i] - hwp) / (ha - hwp))
                                sink[i] = pe * w[i] * f
                        for i in range(N):
                            cap = dn / dts
                            B[i] = C[i] * cap
                            D[i] = C[i] * cap * hk[i] - (th[i] - thn[i]) * cap - sink[i]
                            if i > 0:
                                kf = math.sqrt(K[i] * K[i - 1]) if (K[i] > 0 and K[i - 1] > 0) else 0.0
                                g = kf / dn
                                A[i] = -g; B[i] += g
                                D[i] += -kf * cosb * (-1.0) * 0.0  # placeholder (gravità sotto)
                                D[i] += kf * cosb
                            if i < N - 1:
                                kf = math.sqrt(K[i] * K[i + 1]) if (K[i] > 0 and K[i + 1] > 0) else 0.0
                                g = kf / dn
                                Cc[i] = -g; B[i] += g
                                D[i] += -kf * cosb
                        # top
                        if dirich:
                            g0 = K[0] * 2.0 / dn
                            B[0] += g0; D[0] += K[0] * cosb
                        else:
                            D[0] += qtop
                        # bottom: flusso uscente lam*K*cosb (esplicito)
                        D[N - 1] += -lam * K[N - 1] * cosb
                        hn = _thomas(A, B, Cc, D, N)
                        err = 0.0
                        for i in range(N):
                            if hn[i] != hn[i]: err = 1e9; hn[i] = hk[i]
                            if hn[i] < -1e5: hn[i] = -1e5
                            e = abs(hn[i] - hk[i]) / (1.0 + abs(hn[i]) * 0.02); err = max(err, e)
                            hk[i] = hn[i]
                        if err < 2e-4:
                            conv = True; break
                    if not conv: break
                    if (not dirich) and hk[0] > 0.0:
                        dirich = True; continue
                    break
                if conv:
                    ok = True; break
                dts *= 0.5
                if dts < 20.0: break
            if not ok:
                # rinuncia: accetta ultimo iterato
                pass
            # aggiorna
            qin = qtop
            if dirich:
                K0 = _vg(hk[0], ts[0], tr[0], al[0], nn[0], ks[0])[1]
                qin = K0 * (cosb - 2.0 * hk[0] / dn)
                if qin > qtop: qin = qtop
                if qin < 0.0: qin = 0.0
                run_acc += (qtop - qin) * dts / cosb * 1e3
            for i in range(N):
                h[i] = hk[i] if hk[i] < 0.0 else hk[i]
                thn[i] = _vg(h[i], ts[i], tr[i], al[i], nn[i], ks[i])[0]
            pe = pet[s] * 1e-3 / dt_s * cosb * kc
            for i in range(N):
                if w[i] > 0.0:
                    f = 1.0 if h[i] > ha else (0.0 if h[i] < hwp else (h[i] - hwp) / (ha - hwp))
                    aet_acc += pe * w[i] * f * dts / cosb * 1e3
            tleft -= dts
            if dts < dt_s: dts = min(dts * 2.0, dt_s)
        for i in range(N):
            TH[s, i] = thn[i]; HH[s, i] = h[i]
        RUN[s] = run_acc; AET[s] = aet_acc
    return TH, HH, RUN, AET


DT_MIN = 60
LPAR = ["lks", "la", "n", "ts", "tr"]          # parametri per strato


def cell_arrays(site, P):
    """P: array (4,5) [log10 Ks m/s, log10 alpha 1/m, n, theta_s, theta_r] per strato -> array per cella"""
    cuts = [0] + SITES[site]["cuts"] + [NCELL]
    ks = np.empty(NCELL); al = np.empty(NCELL); nn = np.empty(NCELL); ts = np.empty(NCELL); tr = np.empty(NCELL)
    for L in range(4):
        sl = slice(cuts[L], cuts[L + 1])
        ks[sl] = 10 ** P[L, 0]; al[sl] = 10 ** P[L, 1]; nn[sl] = P[L, 2]; ts[sl] = P[L, 3]; tr[sl] = P[L, 4]
    return ks, al, nn, ts, tr


def root_weights(site):
    zc = (np.arange(NCELL) + 0.5) * DZ
    w = np.exp(-zc / SITES[site]["zr"]) * (zc <= 0.6)
    return w / w.sum()


_OBS = {}
def observations(site):
    if site in _OBS: return _OBS[site]
    S = SITES[site]; o = {}
    for d in S["dep"]:
        o["th" + d] = rawv(f"{S['st']}_VWC_{d}").loc["2025"].resample(f"{DT_MIN}min").mean()
        o["ps" + d] = rawv(f"{S['st']}_PSI_{d}").loc["2025"].resample(f"{DT_MIN}min").median()
    _OBS[site] = pd.DataFrame(o); return _OBS[site]


def initial_h(site, t0):
    S = SITES[site]; ob = observations(site)
    row = ob.loc[t0:].iloc[:96].median()                      # primi 2 giorni
    s = np.array([max(-row["ps" + d], 1.0) for d in S["dep"]])   # kPa, almeno 1
    zc = (np.arange(NCELL) + 0.5) * DZ
    ls = np.interp(zc, S["z"], np.log10(s))
    return -(10 ** ls) * KPA_M


def simulate(site, P, lam, kc, t0="2025-02-01", t1="2025-12-31 23:30", fcg=None):
    idx, rain, pet = fcg if fcg is not None else forcing(DT_MIN)
    m = (idx >= pd.Timestamp(t0)) & (idx <= pd.Timestamp(t1))
    ks, al, nn, ts, tr = cell_arrays(site, P)
    h0 = initial_h(site, t0)
    TH, HH, RUN, AET = _run(h0, ts, tr, al, nn, ks, lam, kc, root_weights(site), rain[m], pet[m], DT_MIN * 60.0, COSB, DZ, -150.0, -1.0)
    return idx[m], TH, HH, RUN, AET


def sensor_series(site, TH, HH):
    """θ e suzione (kPa) ai sensori: media delle 2 celle adiacenti al confine"""
    S = SITES[site]; out = {}
    for d, z in zip(S["dep"], S["z"]):
        j = int(round(z / DZ))                # confine tra cella j-1 e j
        out["th" + d] = 0.5 * (TH[:, j - 1] + TH[:, j])
        out["ps" + d] = 0.5 * (HH[:, j - 1] + HH[:, j]) / KPA_M       # kPa (negativo)
    return out


def storage(site, TH):
    return (TH[:, :30].sum(axis=1) * DZ * 1e3, TH.sum(axis=1) * DZ * 1e3)   # mm per area orizzontale (0–75, 0–100)
