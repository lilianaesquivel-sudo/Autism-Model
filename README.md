# Autism-Model

MATLAB / GNU Octave code for the paper

> L. Esquivel, N. Agudelo, D. E. Pico-Lache, M. Mesa-Mazo.
> *A state-dependent dynamical model of sensory regulation during autism therapy.*

The code reproduces every numerical result of the paper: the illustrative
example of Section 5.1 (Table 7, Figure 4), the simulations of Section 7
(Figures 5–9), the sensitivity analysis of Section 7.1, and the eigenvalues
quoted in the Discussion.

## Quick start

In MATLAB (R2016b or later) or GNU Octave (6 or later):

```matlab
cd src
run_all
```

`run_all` prints each computed value next to the value reported in the paper
and saves Figures 5–9 as PNG files in `figures/`. It takes about one minute in
MATLAB and a few minutes in Octave; most of the time is the sensitivity
analysis (2000 parameter samples). To skip the figures or reduce the number of
samples:

```matlab
SAVE_FIGS = false; NSAMP = 200; run_all
```

No toolboxes are needed (only `ode45`, `fzero`, `eig`, `trapz`).

## The model

State: arousal `A`, regulatory capacity `R`, participation `P`, all in [0, 1];
time in minutes. With inputs `SH` (aversive), `SZ` (seeking-related) and the
target level `Asigma`,

```
dA/dt = alphaH*SH*(1-A) + alphaZ*SZ*(Asigma-A) + eta*Zbar*(1-Gamma(A))*(Ast-A) + beta*R*(Ab-A)
dR/dt = mu*(1-R) - nuH*SH*R - nuA*max(A-Ad,0)*R
dP/dt = lambda*R*Gamma(A)*(1-P) - delta*(1-Gamma(A))*P
Gamma(A) = exp(-((A-Ast)/kappa)^2),   Ad = Ast + kappa
```

The inputs are piecewise constant, so the system is integrated separately on
each interval between switching times (`ode45`, RelTol 1e-8, AbsTol 1e-10),
the state at the end of one interval being the initial condition of the next
(the Carathéodory solution). Equilibria are the zeros of `F(A) = f(A, R(A))`.

## Files (`src/`)

| File | Content |
|---|---|
| `params_default.m` | Parameter values (Table 5) |
| `child_profile.m` | Profiles of Child A and Child B (Table 4), `Zbar`, `Ab` |
| `make_event.m`, `session_events.m` | Events of the illustrative example and of the session of Section 7.2 |
| `event_contributions.m` | Registration and individualized contributions `E^H`, `E^Z` |
| `inputs_from_events.m` | Aggregation, modulatory shares and target level: `[SH, SZ, Asigma]` |
| `gamma_band.m`, `capacity_nullcline.m` | `Gamma(A)` and the equilibrium capacity `R(A)` |
| `model_rhs.m` | Right-hand side of the model |
| `simulate_pw.m`, `sim_profile.m`, `sim_session.m` | Piecewise integration |
| `equilibria_const.m`, `resting_state.m` | Equilibria, their stability, resting state |
| `count_zeros.m`, `bistable_range.m` | Bistable range of the aversive input |
| `prop63_explicit.m` | Explicit criterion of Proposition 6.3 (`eta = 0`) |
| `example_inputs.m` | Section 5.1, Table 7, Figure 4 |
| `fig_session.m` | Section 7.2, Figure 5 |
| `fig_bifurcation.m` | Section 7.3, Figure 6 (bistable ranges, hysteresis) |
| `fig_memory.m` | Section 7.4, Figure 7 |
| `fig_timing.m`, `timing_setup.m` | Section 7.5, Figure 8 |
| `fig_window.m` | Section 7.5, Figure 9 (timing versus intensity of control) |
| `sensitivity.m` | Section 7.1, sensitivity of the bistable range |
| `eigen_discussion.m` | Eigenvalues quoted in the Discussion |
| `run_all.m` | Runs everything and compares with the paper |

Domains are numbered 1 auditory, 2 visual, 3 tactile, 4 gustatory/olfactory,
5 vestibular/proprioceptive.

## Notes

* The joint sensitivity analysis uses random samples (seed 1), so its
  percentages agree with the paper up to sampling variation.
* All parameter values are illustrative (see Table 5 of the paper).

## License

MIT
