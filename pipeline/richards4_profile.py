"""
richards4_profile.py — profilo di identificabilità di log10 Ks per strato.

uso: python pipeline/richards4_profile.py NV|V
Per ogni strato L e per una griglia di log10 Ks fissati, ri-ottimizza tutti gli altri parametri (partendo dal miglior
adattamento di richards4_calib.py) e registra il residuo normalizzato. Un profilo piatto = Ks non identificabile
(solo limite inferiore); un minimo netto = Ks vincolato. Scrive data/richards4_profile_<sito>.json
"""
import sys, os, json, time
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import richards4 as R

site = sys.argv[1]
S = R.SITES[site]
FCG = R.forcing(R.DT_MIN); OB = R.observations(site)
T0, T1 = "2025-02-01", "2025-11-15 23:00"
J = json.load(open(os.path.join(R.ROOT, "data", f"richards4_{site}.json")))
best = min(J["starts"], key=lambda s: s["rms"]); x0 = np.array(best["x"])
lo, hi = np.array(J["bounds_lo"]), np.array(J["bounds_hi"])


def resid(x):
    P = x[:20].reshape(4, 5)
    idx, TH, HH, RUN, AET = R.simulate(site, P, x[20], x[21], t0=T0, t1=T1, fcg=FCG)
    sp = R.sensor_series(site, TH, HH); ob = OB.reindex(idx); res = []
    for d in S["dep"]:
        res.append(((sp["th" + d] - ob["th" + d].values) / 0.02)[::3])
        sm = np.log10(np.clip(-sp["ps" + d], 1.0, 100.0)); so = np.log10(np.clip(-ob["ps" + d].values, 1.0, 100.0))
        res.append(((sm - so) / 0.3)[::3])
    r = np.concatenate(res); r = np.where(np.isfinite(r), r, 0.0)
    Pm = P; pr = []
    for L in range(1, 4):
        pr += [(Pm[L, 0] - Pm[L - 1, 0]) / 2.0 * 3, (Pm[L, 1] - Pm[L - 1, 1]) / 1.0 * 3, (Pm[L, 2] - Pm[L - 1, 2]) / 1.0 * 3]
    return np.concatenate([r, pr]), len(r)


out = {"site": site, "best_rms": best["rms"], "grid": {}}
nobs = resid(x0)[1]
for L in range(3):                         # strati monitorati; il 4° è confuso con λ per costruzione
    jk = L * 5; pts = []
    for v in [-8.0, -7.0, -6.0, -5.0, -4.0]:
        free = np.ones(len(x0), bool); free[jk] = False
        xs = x0.copy(); xs[jk] = v
        def f(z):
            xx = xs.copy(); xx[free] = z; return resid(xx)[0]
        t = time.time()
        sol = least_squares(f, x0[free], bounds=(lo[free], hi[free]), x_scale=((hi - lo) / 4)[free], diff_step=2e-3, max_nfev=12)
        rms = float(np.sqrt(np.mean(sol.fun[:nobs] ** 2)))
        pts.append({"lks": v, "rms": rms})
        print(site, "strato", L + 1, "lks", v, "rms %.3f" % rms, "%.0fs" % (time.time() - t), flush=True)
        out["grid"]["L%d" % (L + 1)] = pts
        json.dump(out, open(os.path.join(R.ROOT, "data", f"richards4_profile_{site}.json"), "w"))
