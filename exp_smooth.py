"""Experiment 3: smoothed depletion threshold and capacity noise.
(a) Replace (A-A_d)_+ by w*log(1+exp((A-A_d)/w)) and recompute the
    folds, the slow return rate and the normalized linear arousal
    variance on the regulated branch near onset.
(b) Nonsmooth model with additional noise sigR*R(1-R)dB2 in capacity:
    normalized linear arousal variance on both branches."""
import numpy as np, json, warnings
from model import Model, refine_fold, A_D
warnings.filterwarnings("ignore")

dgrid = np.logspace(-6, np.log10(0.1), 120)
out = {"d": dgrid.tolist(), "smooth": {}, "capnoise": {}}


def reg_root(m, s):
    r = m.equilibria(s, n=20001)
    return r[0] if len(r) else None


for w in [0.0, 0.0025, 0.005, 0.01, 0.02]:
    m = Model("A", w=w)
    if w == 0.0:
        s_on = m.s_on()
    else:
        _, s_on = refine_fold(m, m.s_on())
    _, s_rec = refine_fold(m, m.s_rec())
    rate, var = [], []
    for d in dgrid:
        s = s_on - d
        A = reg_root(m, s)
        if A is None or (w == 0.0 and A >= A_D):
            rate.append(np.nan); var.append(np.nan); continue
        J = m.jac(A, m.Rnull(A, s), s)
        ev = np.linalg.eigvals(J)
        rate.append(float(min(abs(ev.real))))
        var.append(float(m.lin_var(A, s, 1.0)[0, 0]))
    out["smooth"][str(w)] = dict(s_on=s_on, s_rec=s_rec, rate=rate, var=var)
    rate = np.array(rate); var = np.array(var)
    print(f"w={w:<7} s_on={s_on:.4f} s_rec={s_rec:.4f}  "
          f"rate(d=0.05)={np.interp(0.05, dgrid, rate):.4f}  "
          f"var(d=0.05)={np.interp(0.05, dgrid, var):.4f}")

# distance at which the smoothed rate falls to half of the nonsmooth one,
# and the variance doubles relative to the nonsmooth one
ref = out["smooth"]["0.0"]
for w in [0.0025, 0.005, 0.01, 0.02]:
    r = np.array(out["smooth"][str(w)]["rate"])
    v = np.array(out["smooth"][str(w)]["var"])
    rr = np.array(ref["rate"]); vr = np.array(ref["var"])
    ok = ~np.isnan(r)
    half = dgrid[ok][np.argmax(r[ok] >= 0.5 * rr[ok])]
    dbl = dgrid[ok][np.argmax(v[ok] <= 2 * vr[ok])]
    out["smooth"][str(w)].update(d_half=float(half), d_double=float(dbl))
    print(f"w={w}: rate < half of nonsmooth only for d < {half:.2e}"
          f" ({half / w:.2f} w);  variance > 2x only for d < {dbl:.2e} ({dbl / w:.2f} w)")

# (b) capacity noise in the nonsmooth model
m = Model("A")
s_on = m.s_on(); _, s_rec = refine_fold(m, m.s_rec())
for sigR in [0.0, 0.5, 1.0]:
    res = {}
    for br in ["reg", "dys"]:
        vv = []
        for d in dgrid:
            s = s_on - d if br == "reg" else s_rec + d
            r = m.equilibria(s, n=20001)
            if len(r) < 3:
                vv.append(np.nan); continue
            A = r[0] if br == "reg" else r[-1]
            # arousal noise 1 (normalization), capacity noise sigR/sig
            vv.append(float(m.lin_var(A, s, 1.0, sigR)[0, 0]))
        res[br] = vv
    out["capnoise"][str(sigR)] = res
    vr = np.array(res["reg"]); vd = np.array(res["dys"])
    print(f"sigR/sig={sigR}: reg var range [{np.nanmin(vr):.4f},{np.nanmax(vr):.4f}]"
          f"  dys var at d=1e-4: {np.interp(1e-4, dgrid, vd):.3f}")
json.dump(out, open("exp_smooth.json", "w"))
