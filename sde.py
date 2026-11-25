"""Euler--Maruyama integration of the nonlinear arousal--capacity SDE
    dA = f(A,R,s) dtau + sig A(1-A) dB1,
    dR = g(A,R,s) dtau + sigR R(1-R) dB2,
vectorized over any array shape (ensemble members x parameter values).
"""
import numpy as np


def simulate(m, s_fun, A0, R0, sig, T, dt=0.01, every=50, sigR=0.0,
             seed=0, stop=None):
    """Integrate on [0,T]; s_fun(tau) returns s (scalar or array matching A0).
    Records A, R every `every` steps. Returns (tau, A_rec, R_rec)."""
    rng = np.random.default_rng(seed)
    A = np.array(A0, dtype=float).copy()
    R = np.array(R0, dtype=float).copy()
    n = int(round(T / dt))
    nrec = n // every + 1
    Arec = np.empty((nrec,) + A.shape)
    Rrec = np.empty((nrec,) + A.shape)
    Arec[0], Rrec[0] = A, R
    sq = np.sqrt(dt)
    k = 1
    for i in range(1, n + 1):
        s = s_fun((i - 1) * dt)
        dW = rng.standard_normal(A.shape) * sq
        An = A + m.f(A, R, s) * dt + sig * A * (1 - A) * dW
        Rn = R + m.g(A, R, s) * dt
        if sigR > 0:
            Rn += sigR * R * (1 - R) * rng.standard_normal(A.shape) * sq
        A = np.clip(An, 1e-9, 1 - 1e-9)
        R = np.clip(Rn, 1e-9, 1 - 1e-9)
        if i % every == 0:
            Arec[k], Rrec[k] = A, R
            k += 1
    tau = np.arange(nrec) * every * dt
    return tau, Arec, Rrec
