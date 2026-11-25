"""Experiment 2: nonlinear SDE under slow ramps of the aversive input
(Child A, a_E = 0.07), as in the hysteresis experiment of the paper:
s increases from 0 to 0.8 (or decreases from 0.8 to 0) over 400 min.
Rolling variance and lag-1 autocorrelation of detrended arousal are
computed before each transition and summarized by Kendall's tau."""
import numpy as np, json, warnings
from scipy.ndimage import gaussian_filter1d
from scipy.stats import kendalltau
from model import Model, BETA
from sde import simulate
warnings.filterwarnings("ignore")

m = Model("A")
N, dt, every = 200, 0.01, 50                 # record every 0.5 tau = 0.1 min
T_min = 400.0
T = T_min * BETA
rate = 0.8 / T                               # per tau
SAMPLE_MIN = every * dt / BETA               # 0.1 min

A_UP, A_DOWN = 0.70, 0.62                    # transition detection on arousal
SMOOTH_MIN = 1.0                             # 1-min moving average for detection
LAG_MIN = 1.0                                # autocorrelation lag (min)
WIN_MIN, PRE_MIN = 20.0, 100.0               # rolling window, analysed segment
BW_LIST = [5.0, 20.0]                        # Gaussian detrending bandwidth (min)


def det_run(direction):
    """Noise-free ramp: transition input."""
    if direction == "up":
        A0, R0 = m.equilibria(0.0)[0], 1.0
        sf = lambda t: rate * t
    else:
        e = m.equilibria(0.8)[-1]
        A0, R0 = e, m.Rnull(e, 0.8)
        sf = lambda t: 0.8 - rate * t
    tau, A, R = simulate(m, sf, np.array([A0]), np.array([R0]), 0.0, T, dt, every)
    A = A[:, 0]
    k = np.argmax(A > A_UP) if direction == "up" else np.argmax(A < A_DOWN)
    return float(sf(tau[k]))


def rolling(x, w, L):
    """Rolling variance and lag-L autocorrelation, window w samples."""
    n = len(x) - w + 1
    V = np.empty(n); AC = np.empty(n)
    for i in range(n):
        y = x[i:i + w]
        y = y - y.mean()
        V[i] = np.mean(y ** 2)
        AC[i] = np.sum(y[L:] * y[:-L]) / np.sqrt(np.sum(y[L:] ** 2) * np.sum(y[:-L] ** 2))
    return V, AC


out = {"N": N, "T_min": T_min, "win_min": WIN_MIN, "pre_min": PRE_MIN,
       "s_det": {}, "res": {}}
for direction in ["up", "down"]:
    out["s_det"][direction] = det_run(direction)

w = int(round(WIN_MIN / SAMPLE_MIN))
L = int(round(LAG_MIN / SAMPLE_MIN))
npre = int(round(PRE_MIN / SAMPLE_MIN))
for sig in [0.02, 0.05]:
    for direction in ["up", "down"]:
        if direction == "up":
            A0 = np.full(N, m.equilibria(0.0)[0]); R0 = np.ones(N)
            sf = lambda t: rate * t
        else:
            e = m.equilibria(0.8)[-1]
            A0 = np.full(N, e); R0 = np.full(N, m.Rnull(e, 0.8))
            sf = lambda t: 0.8 - rate * t
        tau, A, R = simulate(m, sf, A0, R0, sig, T, dt, every,
                             seed=101 if direction == "up" else 202)
        s_path = sf(tau)
        k_s = int(round(SMOOTH_MIN / SAMPLE_MIN))
        As = np.apply_along_axis(lambda c: np.convolve(c, np.ones(k_s) / k_s, 'same'), 0, A)
        crossed = (As > A_UP) if direction == "up" else (As < A_DOWN)
        crossed[:k_s] = False
        kc = np.where(crossed.any(0), crossed.argmax(0), -1)
        res = {"s_trans": [], "tau": {}, "ratio": {}, "curves": {}}
        for bw in BW_LIST:
            res["tau"][str(bw)] = {"V": [], "AC": []}
            res["ratio"][str(bw)] = {"V": [], "AC": []}
            res["curves"][str(bw)] = {"s": [], "V": [], "AC": []}
        for i in range(N):
            if kc[i] < 0:
                continue
            res["s_trans"].append(float(s_path[kc[i]]))
            # analysed segment: PRE_MIN before the transition, ending 5 min
            # before it (avoid the transition itself)
            end = kc[i] - int(round(5.0 / SAMPLE_MIN))
            start = end - npre - w
            if start < 0:
                continue
            x = A[start:end, i]
            for bw in BW_LIST:
                trend = gaussian_filter1d(x, bw / SAMPLE_MIN, mode="nearest")
                V, AC = rolling(x - trend, w, L)
                tt = np.arange(len(V))
                res["tau"][str(bw)]["V"].append(kendalltau(tt, V)[0])
                res["tau"][str(bw)]["AC"].append(kendalltau(tt, AC)[0])
                res["ratio"][str(bw)]["V"].append(float(V[-1] / V[0]))
                res["ratio"][str(bw)]["AC"].append(float(AC[-1] - AC[0]))
                # indicators against input level at the window end
                se = s_path[start + w - 1:end]
                res["curves"][str(bw)]["s"].append(se[::50].tolist())
                res["curves"][str(bw)]["V"].append((V[::50] / sig ** 2).tolist())
                res["curves"][str(bw)]["AC"].append(AC[::50].tolist())
        st = np.array(res["s_trans"])
        print(f"{direction} sig={sig}: transitions {len(st)}/{N}, "
              f"s_trans median {np.median(st):.3f} [{np.percentile(st,5):.3f},"
              f"{np.percentile(st,95):.3f}]  (deterministic {out['s_det'][direction]:.3f})")
        for bw in BW_LIST:
            tV = np.array(res["tau"][str(bw)]["V"]); tA = np.array(res["tau"][str(bw)]["AC"])
            rV = np.array(res["ratio"][str(bw)]["V"]); rA = np.array(res["ratio"][str(bw)]["AC"])
            print(f"   bw={bw:>4} min  n={len(tV)}  Kendall V: median {np.median(tV):+.2f}"
                  f" (frac>0.5 {np.mean(tV>0.5):.2f})  AC: median {np.median(tA):+.2f}"
                  f" (frac>0.5 {np.mean(tA>0.5):.2f})  V_end/V_start median {np.median(rV):.2f}"
                  f"  dAC median {np.median(rA):+.3f}")
        out["res"][f"{direction}_{sig}"] = res
json.dump(out, open("exp_ramp.json", "w"))
