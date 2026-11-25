"""Experiment 1: nonlinear SDE at fixed aversive input on both branches.
Compares the simulated stationary arousal variance and lag-1
autocorrelation with the linear-noise (Lyapunov) prediction, and
records the fraction of paths leaving the basin. Child A, a_E = 0.07."""
import numpy as np, json, time, warnings
from model import Model, refine_fold, A_D
from sde import simulate
warnings.filterwarnings("ignore")

m = Model("A")
s_on = m.s_on()
A_rec, s_rec = refine_fold(m, m.s_rec())
d = np.array([0.002, 0.004, 0.007, 0.01, 0.015, 0.02, 0.03, 0.05, 0.08, 0.12])
N, T, dt, every, burn = 100, 2000.0, 0.01, 50, 200.0
L = 10                    # lag in samples
lag = L * every * dt      # 5 dimensionless = 1 min

out = {"s_on": s_on, "s_rec": s_rec, "d": d.tolist(), "N": N, "T": T,
       "lag": lag, "res": {}}

for sig in [0.02, 0.05]:
    for branch in ["reg", "dys"]:
        s_vals = s_on - d if branch == "reg" else s_rec + d
        eq = [m.equilibria(s) for s in s_vals]
        A_st = np.array([e[0] if branch == "reg" else e[-1] for e in eq])
        A_sd = np.array([e[1] for e in eq])           # saddle
        R_st = m.Rnull(A_st, s_vals)
        R_sd = m.Rnull(A_sd, s_vals)
        S = np.repeat(s_vals[:, None], N, 1)
        A0 = np.repeat(A_st[:, None], N, 1)
        R0 = np.repeat(R_st[:, None], N, 1)
        t0 = time.time()
        tau, A, R = simulate(m, lambda t: S, A0, R0, sig, T, dt, every,
                             seed=11 if branch == "reg" else 23)
        # first escape: entering the opposite quadrant of the saddle
        if branch == "reg":
            esc = (A > A_sd[None, :, None]) & (R < R_sd[None, :, None])
        else:
            esc = (A < A_sd[None, :, None]) & (R > R_sd[None, :, None])
        k_esc = np.where(esc.any(0), esc.argmax(0), len(tau))
        valid = (np.arange(len(tau))[:, None, None] < k_esc[None]) & \
                (tau[:, None, None] >= burn)
        V, AC, frac, tesc = [], [], [], []
        for j in range(len(d)):
            x = A[:, j, :]
            v = valid[:, j, :]
            mu = x[v].mean()
            y = np.where(v, x - mu, np.nan)
            V.append(float(np.nanmean(y ** 2)))
            pair = v[L:] & v[:-L]
            AC.append(float(np.sum((y[L:] * y[:-L])[pair]) /
                            np.sqrt(np.sum(y[L:][pair] ** 2) *
                                    np.sum(y[:-L][pair] ** 2))))
            frac.append(float(np.mean(k_esc[j] < len(tau))))
            ke = k_esc[j][k_esc[j] < len(tau)]
            tesc.append(float(np.median(tau[ke]) / 5.0) if len(ke) else None)
        # linear predictions
        Vlin, AClin = [], []
        for j, s in enumerate(s_vals):
            Q = m.lin_var(A_st[j], s, sig)
            Rj = R_st[j]
            J = m.jac(A_st[j], Rj, s)
            from scipy.linalg import expm
            Vlin.append(float(Q[0, 0]))
            AClin.append(float((expm(J * lag) @ Q)[0, 0] / Q[0, 0]))
        out["res"][f"{branch}_{sig}"] = dict(
            s=s_vals.tolist(), A=A_st.tolist(), V=V, AC=AC, esc=frac,
            Vlin=Vlin, AClin=AClin, tesc_min=tesc, gap=(np.abs(A_sd - A_st)).tolist())
        print(branch, sig, "%.0fs" % (time.time() - t0))
        for j in range(len(d)):
            print("  d=%.3f  V/Vlin=%.3f  AC=%.3f (lin %.3f)  esc=%.2f  t_esc=%s" %
                  (d[j], V[j] / Vlin[j], AC[j], AClin[j], frac[j], tesc[j]))

json.dump(out, open("exp_stationary.json", "w"), indent=1)
