# Nonlinear stochastic simulations and smoothed-threshold checks

Python code for Sections 5.8-5.9 (Figures 9-11) of
"Bistability and Asymmetric Critical Slowing Down in a State-Dependent
Model of Sensory Regulation in Autism".

Requirements: Python 3.10+, numpy, scipy, matplotlib.

    python exp_stationary.py   # fixed inputs, both branches  -> exp_stationary.json  (~1 min)
    python exp_ramp.py         # slow ramps, rolling indicators -> exp_ramp.json      (~1 min)
    python exp_smooth.py       # softplus threshold, capacity noise -> exp_smooth.json (seconds)
    python figures.py          # fig_sde_stationary.pdf, fig_sde_ramps.pdf, fig_smooth_threshold.pdf

model.py  : nondimensional arousal-capacity subsystem (Child A/B parameters),
            equilibria, folds, Jacobian, Lyapunov covariance; w>0 smooths the threshold.
sde.py    : vectorized Euler-Maruyama integrator (dtau = 0.01 = 0.12 s).
Random seeds are fixed, so the numbers quoted in the paper are reproducible.
