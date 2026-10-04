"""
richards4_calib.py — calibrazione inversa del modello a 4 strati (richards4.py) su θ e ψ osservati.

uso: python pipeline/richards4_calib.py NV|V [n_starts]
Scrive: data/richards4_<sito>.json  (parametri, incertezza, serie simulate/osservate, bilancio)
Finestra di calibrazione 2025-02-01 → 2025-11-15; validazione 2025-11-15 → 2025-12-31 (stato iniziale dalle osservazioni).
"""
import sys, os, json, time
import numpy as np, pandas as pd
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import richards4 as R

site = sys.argv[1]; nstart = int(sys.argv[2]) if len(sys.argv) > 2 else 3
S = R.SITES[site]
FCG = R.forcing(R.DT_MIN)
OB = R.observations(site)
T0, T1, TV0, TV1 = "2025-02-01", "2025-11-15 23:00", "2025-11-15", "2025-12-31 23:00"
thmax = [float(OB["th" + d].quantile(0.995)) for d in S["dep"]]
thmin = [float(OB["th" + d].quantile(0.005)) for d in S["dep"]]
NAMES = ["lks", "la", "n", "ts", "tr"]
lo, hi = [], []
for L in range(4):
    tsm = thmax[min(L, 2)]
    lo += [-9.0, -0.5, 1.1, tsm if L < 3 else 0.15, 0.0]
    hi += [-4.0, 2.0, 3.5, tsm + 0.10 if L < 3 else 0.45, min(0.10, 0.9 * thmin[min(L, 2)])]
lo += [0.005, 0.2]; hi += [0.5, 1.5]          # lambda, Kc
lo, hi = np.array(lo), np.array(hi)


def unpack(x):
    return x[:20].reshape(4, 5), x[20], x[21]


def start(lk):
    x = []
    for L in range(4):
        x += [lk, np.log10(0.35 * 9.81), 1.5, thmax[min(L, 2)] + 0.03 if L < 3 else thmax[2], 0.5 * min(0.1, thmin[min(L, 2)])]
    x += [0.05, 0.8]
    return np.clip(np.array(x), lo + 1e-6, hi - 1e-6)


def residuals_for(x, t0, t1, full=False):
    P, lam, kc = unpack(x)
    idx, TH, HH, RUN, AET = R.simulate(site, P, lam, kc, t0=t0, t1=t1, fcg=FCG)
    sp = R.sensor_series(site, TH, HH)
    ob = OB.reindex(idx)
    res = []
    for d in S["dep"]:
        a = (sp["th" + d] - ob["th" + d].values) / 0.02
        res.append(a[::3])
        sm = np.log10(np.clip(-sp["ps" + d], 1.0, 100.0)); so_raw = -ob["ps" + d].values
        so = np.log10(np.clip(so_raw, 1.0, 100.0))
        b = (sm - so) / 0.3
        res.append(b[::3])
    r = np.concatenate(res)
    r = np.where(np.isfinite(r), r, 0.0)
    if full: return r, (idx, TH, HH, RUN, AET, sp, ob)
    return r


def prior(x):
    P, lam, kc = unpack(x)
    out = []
    for L in range(1, 4):
        out += [(P[L, 0] - P[L - 1, 0]) / 2.0 * 3, (P[L, 1] - P[L - 1, 1]) / 1.0 * 3, (P[L, 2] - P[L - 1, 2]) / 1.0 * 3]
    return np.array(out)


def fun(x):
    r = np.concatenate([residuals_for(x, T0, T1), prior(x)])
    return r


best = None; runs = []
for k, lk in enumerate([-6.0, -5.0, -7.0, -5.5][:nstart]):
    t = time.time()
    try:
        sol = least_squares(fun, start(lk), bounds=(lo, hi), x_scale=(hi - lo) / 4, diff_step=2e-3, max_nfev=45, verbose=0)
    except Exception as e:
        print("start", k, "fallito", e, flush=True); continue
    cost = float(np.sqrt(np.mean(sol.fun[:-9] ** 2)))
    print(f"start {k} lk0={lk}: rms residuo normalizzato {cost:.3f}  nfev {sol.nfev}  {time.time()-t:.0f}s", flush=True)
    runs.append((cost, sol.x))
    if best is None or sol.cost < best.cost: best = sol

x = best.x
# incertezza (matrice di covarianza approssimata dalla Jacobiana)
J = best.jac; dof = max(len(best.fun) - len(x), 1)
s2 = 2 * best.cost / dof
try:
    cov = np.linalg.pinv(J.T @ J) * s2; sd = np.sqrt(np.clip(np.diag(cov), 0, None))
except Exception:
    sd = np.full(len(x), np.nan)

# serie complete + validazione
rc, (idx, TH, HH, RUN, AET, sp, ob) = residuals_for(x, T0, TV1, full=True)
P, lam, kc = unpack(x)
rv, (idv, THv, HHv, RUNv, AETv, spv, obv) = residuals_for(x, TV0, TV1, full=True)


def rmse_theta(sp_, ob_, a, b):
    out = {}
    for d in S["dep"]:
        m = np.isfinite(ob_["th" + d].values)
        out[d] = float(np.sqrt(np.mean((sp_["th" + d][m] - ob_["th" + d].values[m]) ** 2)))
    return out


cal = (idx >= pd.Timestamp(T0)) & (idx <= pd.Timestamp(T1))
spc = {k: v[cal] for k, v in sp.items()}; obc = ob[cal]
ncal = int(cal.sum())
S75, S100 = R.storage(site, TH)
rain = FCG[1][(FCG[0] >= pd.Timestamp(T0)) & (FCG[0] <= pd.Timestamp(TV1))]
bal = {"P_mm": float(rain.sum()), "runoff_mm": float(RUN.sum()), "aet_mm": float(AET.sum()),
       "dS100_mm": float(S100[-1] - S100[0]), "leak_mm": float(rain.sum() - RUN.sum() - AET.sum() - (S100[-1] - S100[0]))}
def cc(a, b):
    m = np.isfinite(a) & np.isfinite(b); return float(np.corrcoef(a[m], b[m])[0, 1]) if m.sum() > 5 else None
out = {
    "site": site, "dt_min": R.DT_MIN, "calib": [T0, T1], "valid": [TV0, TV1], "beta_deg": 30, "layers_cm": [0] + [c * 2.5 for c in S["cuts"]] + [100],
    "params": [{"layer": L + 1, **{n: float(P[L, j]) for j, n in enumerate(NAMES)},
                **{n + "_sd": float(sd[L * 5 + j]) for j, n in enumerate(NAMES)}} for L in range(4)],
    "lambda": float(lam), "lambda_sd": float(sd[20]), "kc": float(kc), "kc_sd": float(sd[21]),
    "rms_norm": float(np.sqrt(np.mean(best.fun[:-9] ** 2))), "starts": [{"rms": c, "x": list(map(float, xx))} for c, xx in runs],
    "rmse_theta_cal": {d: float(np.sqrt(np.nanmean((spc["th" + d] - obc["th" + d].values) ** 2))) for d in S["dep"]},
    "rmse_theta_val": rmse_theta(spv, obv, 0, 0),
    "r_theta_val": {d: cc(spv["th" + d], obv["th" + d].values) for d in S["dep"]},
    "balance": bal, "bounds_lo": list(map(float, lo)), "bounds_hi": list(map(float, hi)),
}
dd = pd.DatetimeIndex(idx)
daily = pd.DataFrame({"S75": S75, "S100": S100}, index=dd).resample("1D").mean()
ser = {"t": [x_.strftime("%Y-%m-%dT%H:%M") for x_ in dd[::6]]}
for d in S["dep"]:
    ser["th_m_" + d] = [None if not np.isfinite(v) else round(float(v), 4) for v in sp["th" + d][::6]]
    ser["th_o_" + d] = [None if not np.isfinite(v) else round(float(v), 4) for v in ob["th" + d].values[::6]]
    ser["ps_m_" + d] = [round(float(max(v, -2000)), 1) for v in sp["ps" + d][::6]]
    ser["ps_o_" + d] = [None if not np.isfinite(v) else round(float(v), 1) for v in ob["ps" + d].values[::6]]
ser["S75"] = [round(float(v), 1) for v in S75[::6]]; ser["S100"] = [round(float(v), 1) for v in S100[::6]]
ser["runoff_mm_d"] = [round(float(v), 2) for v in pd.Series(RUN, index=dd).resample("1D").sum().values]
ser["runoff_t"] = [x_.strftime("%Y-%m-%d") for x_ in pd.Series(RUN, index=dd).resample("1D").sum().index]
out["series"] = ser
with open(os.path.join(R.ROOT, "data", f"richards4_{site}.json"), "w") as fh:
    json.dump(out, fh, allow_nan=False)
print(json.dumps({k: v for k, v in out.items() if k not in ("series", "starts", "bounds_lo", "bounds_hi")}, indent=1))
