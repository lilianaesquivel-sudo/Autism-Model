"""Figures for the stochastic and smoothing experiments."""
import numpy as np, json, warnings
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import Model, BETA
from sde import simulate
warnings.filterwarnings("ignore")

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
    "axes.titlesize": 9, "axes.labelsize": 9, "legend.fontsize": 7.5,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "lines.linewidth": 1.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#e6e6e6", "grid.linewidth": 0.6,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02})
REG, DYS, INK, MUTED = "#2a78d6", "#eb6834", "#222222", "#8a8a8a"
W = 6.5


def label(ax, s):
    ax.text(-0.16, 1.04, s, transform=ax.transAxes, fontweight="bold",
            fontsize=10, va="bottom")


# ---------------------------------------------------------------- Fig. A
st = json.load(open("exp_stationary.json"))
d = np.array(st["d"])
fig, axs = plt.subplots(1, 3, figsize=(W, 2.2))
for br, col, name in [("reg", REG, "regulated"), ("dys", DYS, "dysregulated")]:
    for sig, mk in [(0.02, "o"), (0.05, "^")]:
        r = st["res"][f"{br}_{sig}"]
        V = np.array(r["V"]) / sig ** 2
        AC = np.array(r["AC"]); esc = np.array(r["esc"])
        ok = esc <= 0.5
        if sig == 0.02:
            axs[0].plot(d, np.array(r["Vlin"]) / sig ** 2, color=col, lw=1.2,
                        label=f"{name}, linear")
            axs[1].plot(d, r["AClin"], color=col, lw=1.2)
        axs[0].plot(d[ok], V[ok], mk, ms=5, mfc="white", mec=col, mew=1.1,
                    ls="none", label=fr"{name}, SDE $\hat\varsigma={sig}$")
        axs[1].plot(d[ok], AC[ok], mk, ms=5, mfc="white", mec=col, mew=1.1, ls="none")
        axs[2].plot(d, esc, marker=mk, ms=4.5, color=col, lw=1.0, mfc="white",
                    mew=1.1, ls="-" if sig == 0.02 else "--")
for ax in axs:
    ax.set_xscale("log")
    ax.set_xlabel("distance to the fold")
axs[0].set_ylabel(r"$V_{\mathcal{A}}/\hat\varsigma^{2}$")
axs[1].set_ylabel("autocorrelation, lag 1 min")
axs[2].set_ylabel("fraction escaped in 400 min")
axs[1].set_ylim(0, 0.85)
axs[2].set_ylim(-0.03, 1.05)
for a, s in zip(axs, "abc"):
    label(a, f"({s})")
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, frameon=False,
           bbox_to_anchor=(0.5, -0.17))
# legend for panel (c) line styles
axs[2].plot([], [], color=MUTED, ls="-", label=r"$\hat\varsigma=0.02$")
axs[2].plot([], [], color=MUTED, ls="--", label=r"$\hat\varsigma=0.05$")
axs[2].legend(frameon=False, loc="lower left")
fig.tight_layout(w_pad=1.4)
fig.savefig("fig_sde_stationary.pdf")
plt.close(fig)

# ---------------------------------------------------------------- Fig. B
rp = json.load(open("exp_ramp.json"))
m = Model("A")
fig, axs = plt.subplots(1, 3, figsize=(W, 2.3), gridspec_kw={"width_ratios": [1.1, 1.1, 1]})
# (a) equilibria + sample paths
ss = np.linspace(0.0, 0.8, 801)
lo, sad, hi = [], [], []
for s in ss:
    e = m.equilibria(s)
    if len(e) == 3:
        lo.append((s, e[0])); sad.append((s, e[1])); hi.append((s, e[2]))
    elif e[0] < 0.65:
        lo.append((s, e[0]))
    else:
        hi.append((s, e[0]))
for pts, ls in [(lo, "-"), (hi, "-"), (sad, ":")]:
    p = np.array(pts)
    axs[0].plot(p[:, 0], p[:, 1], color=MUTED, lw=1.0, ls=ls, zorder=1)
T = 400.0 * BETA; rate = 0.8 / T
n = 3
tau, Au, _ = simulate(m, lambda t: rate * t, np.full(n, m.equilibria(0.0)[0]),
                      np.ones(n), 0.02, T, 0.01, 50, seed=5)
e8 = m.equilibria(0.8)[-1]
tau, Ad, _ = simulate(m, lambda t: 0.8 - rate * t, np.full(n, e8),
                      np.full(n, m.Rnull(e8, 0.8)), 0.02, T, 0.01, 50, seed=6)
for i in range(n):
    axs[0].plot(rate * tau, Au[:, i], color=REG, lw=0.5, alpha=0.8, zorder=2)
    axs[0].plot(0.8 - rate * tau, Ad[:, i], color=DYS, lw=0.5, alpha=0.8, zorder=2)
axs[0].set_xlabel(r"aversive input $s$")
axs[0].set_ylabel(r"arousal $\mathcal{A}$")
axs[0].set_xlim(0, 0.6); axs[0].set_ylim(0.45, 0.85)
axs[0].text(0.22, 0.52, "increasing", color=REG, fontsize=7.5)
axs[0].text(0.43, 0.80, "decreasing", color=DYS, fontsize=7.5)

# (b) rolling variance vs s, ensemble median and IQR
bins = np.linspace(0.12, 0.40, 29)
for direction, col in [("up", REG), ("down", DYS)]:
    c = rp["res"][f"{direction}_0.02"]["curves"]["5.0"]
    S = np.concatenate([np.array(x) for x in c["s"]])
    V = np.concatenate([np.array(x) for x in c["V"]])
    idx = np.digitize(S, bins)
    mid, med, q1, q3 = [], [], [], []
    for k in range(1, len(bins)):
        sel = idx == k
        if sel.sum() > 50:
            mid.append(0.5 * (bins[k - 1] + bins[k]))
            med.append(np.median(V[sel])); q1.append(np.percentile(V[sel], 25))
            q3.append(np.percentile(V[sel], 75))
    axs[1].fill_between(mid, q1, q3, color=col, alpha=0.18, lw=0)
    axs[1].plot(mid, med, color=col, lw=1.5,
                label="increasing" if direction == "up" else "decreasing")
s_on = rp["s_det"]["up"]
axs[1].set_xlabel(r"aversive input $s$ at window end")
axs[1].set_ylabel(r"rolling variance $/\hat\varsigma^{2}$")
axs[1].legend(frameon=False, loc="upper center")
axs[1].set_xlim(0.12, 0.40)

# (c) Kendall tau
groups, data, cols = [], [], []
for ind, nm in [("V", "variance"), ("AC", "autocorr.")]:
    for direction, col in [("up", REG), ("down", DYS)]:
        for sig in [0.02, 0.05]:
            data.append(rp["res"][f"{direction}_{sig}"]["tau"]["5.0"][ind])
            cols.append(col)
pos = [0.7, 1.05, 1.6, 1.95, 3.05, 3.4, 3.95, 4.3]
bp = axs[2].boxplot(data, positions=pos, widths=0.3, patch_artist=True,
                    showfliers=False, medianprops=dict(color=INK, lw=1.2),
                    whiskerprops=dict(color=MUTED), capprops=dict(color=MUTED))
for patch, col, k in zip(bp["boxes"], cols, range(8)):
    patch.set_facecolor(col); patch.set_alpha(0.35 if k % 2 else 0.75)
    patch.set_edgecolor(col)
axs[2].axhline(0, color=MUTED, lw=0.8)
axs[2].set_xticks([1.32, 3.68]); axs[2].set_xticklabels(["variance", "autocorr."])
from matplotlib.patches import Patch
axs[2].legend(handles=[Patch(fc=REG, alpha=0.75, label="increasing"), Patch(fc=DYS, alpha=0.75, label="decreasing"),
                       Patch(fc=MUTED, alpha=0.75, label=r"$\hat\varsigma=0.02$"), Patch(fc=MUTED, alpha=0.35, label=r"$\hat\varsigma=0.05$")],
              frameon=False, loc="lower center", ncol=2, fontsize=6.5, handlelength=1.0, columnspacing=0.8)
axs[2].set_ylabel(r"Kendall $\tau$")
axs[2].set_ylim(-1.25, 1)
axs[2].grid(axis="x", visible=False)
for a, s in zip(axs, "abc"):
    label(a, f"({s})")
fig.tight_layout(w_pad=1.2)
fig.savefig("fig_sde_ramps.pdf")
plt.close(fig)

# ---------------------------------------------------------------- Fig. C
sm = json.load(open("exp_smooth.json"))
dg = np.array(sm["d"])
ramp = ["#86b6ef", "#3987e5", "#256abf", "#104281"]
fig, axs = plt.subplots(1, 2, figsize=(W * 0.78, 2.2))
for (w, col) in zip(["0.0025", "0.005", "0.01", "0.02"], ramp):
    r = sm["smooth"][w]
    axs[0].plot(dg, r["rate"], color=col, lw=1.3, label=fr"$w={float(w):g}$")
    axs[1].plot(dg, r["var"], color=col, lw=1.3)
r0 = sm["smooth"]["0.0"]
axs[0].plot(dg, r0["rate"], color=INK, lw=1.5, ls="--", label="nonsmooth")
axs[1].plot(dg, r0["var"], color=INK, lw=1.5, ls="--")
for ax in axs:
    ax.set_xscale("log"); ax.set_xlabel(r"$s_{\mathrm{on}}(w)-s$")
    ax.set_xlim(1e-6, 0.1)
axs[0].set_yscale("log")
axs[0].set_ylabel("slow return rate")
axs[1].set_ylabel(r"$V_{\mathcal{A}}/\hat\varsigma^{2}$ (linear)")
axs[1].set_yscale("log"); axs[1].set_ylim(0.03, 1.0)
axs[0].legend(frameon=False, loc="lower right")
for a, s in zip(axs, "ab"):
    label(a, f"({s})")
fig.tight_layout(w_pad=1.6)
fig.savefig("fig_smooth_threshold.pdf")
plt.close(fig)
print("done")
