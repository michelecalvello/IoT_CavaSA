"""
richards4_fs.py — FS a 75 cm (e 100 cm) con la pressione dei pori simulata dal modello a 4 strati, per ciascuna delle
soluzioni di calibrazione (le soluzioni con residuo simile ma Ks e drenaggio basale diversi danno un intervallo).
Aggiunge a data/richards4_<sito>.json il campo "fs" (serie giornaliere del minimo di FS e statistiche per scenario).
FS = [c' + (γ z cos²β + σs − u) tanφ'] / (γ z sinβ cosβ), σs = min(Se·s, 20 kPa) se h<0, u = γw·h se h>0.
"""
import sys, os, json, math
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import richards4 as R
fcg = R.forcing(60); c, s = math.cos(R.BETA), math.sin(R.BETA); phi = math.radians(35)
for site in ("NV", "V"):
    fn = os.path.join(R.ROOT, "data", f"richards4_{site}.json"); j = json.load(open(fn)); res = []
    for k, st in enumerate(j["starts"]):
        x = np.array(st["x"]); P = x[:20].reshape(4, 5)
        idx, TH, HH, RUN, AET = R.simulate(site, P, x[20], x[21], fcg=fcg)
        ks, al, nn, ts, tr = R.cell_arrays(site, P)
        o = {"rms": st["rms"], "lambda": float(x[20]), "kc": float(x[21])}
        for zz, nm in ((0.75, "z75"), (1.0, "z100")):
            i = int(round(zz / R.DZ)) - 1; h = HH[:, i]; m = 1 - 1 / nn[i]
            Se = np.where(h < 0, (1 + (al[i] * np.abs(h)) ** nn[i]) ** (-m), 1.0)
            sig = np.where(h < 0, np.minimum(Se * (-h) * 9.81, 20.0), -9.81 * np.maximum(h, 0))
            FS = (1 + (17 * zz * c * c + sig) * math.tan(phi)) / (17 * zz * s * c)
            sr = pd.Series(FS, index=idx); d = sr.resample("1D").min()
            o[nm] = {"fs_min": float(FS.min()), "t_min": str(idx[FS.argmin()]), "h_max_m": float(h.max()), "hours_h_pos": int((h > 0).sum()),
                     "hours_fs_lt1": int((FS < 1).sum()), "daily_min": [round(float(v), 2) for v in d.clip(upper=8).values]}
            if nm == "z75": o["t"] = [t.strftime("%Y-%m-%d") for t in d.index]
        res.append(o)
    j["fs"] = res
    json.dump(j, open(fn, "w"), allow_nan=False)
    print(site, [(round(r["rms"], 3), round(r["z75"]["fs_min"], 2), r["z75"]["hours_h_pos"]) for r in res])
