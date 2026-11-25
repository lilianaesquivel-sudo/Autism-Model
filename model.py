"""Arousal--capacity subsystem of the sensory-regulation model
(nondimensional form, tau = beta * t). Child A defaults.

The depletion switch (A - A_d)_+ can be replaced by the softplus
w*log(1 + exp((A - A_d)/w)); w = 0 gives the original nonsmooth model.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.linalg import solve_continuous_lyapunov

# Physical parameters (Table of the paper), time in minutes
BETA, ALPHA_H, ALPHA_Z, ETA = 5.0, 2.5, 10.0, 2.5
MU, NU_H, NU_A = 0.05, 0.2, 16.0
A_STAR, KAPPA = 0.5, 0.15
A_D = A_STAR + KAPPA

# Nondimensional parameters
AH_HAT = ALPHA_H / BETA      # 0.5
AZ_HAT = ALPHA_Z / BETA      # 2
NUH_HAT = NU_H / MU          # 4
EPS = MU / BETA              # 0.01
THETA_A = NU_A / MU          # 320

CHILD = {"A": dict(A_b=0.5, Zbar=0.14), "B": dict(A_b=0.25, Zbar=0.60)}


def a_E_of(child):
    return ETA * CHILD[child]["Zbar"] / BETA


def Gamma(A):
    return np.exp(-((A - A_STAR) / KAPPA) ** 2)


def h(A):
    return (1.0 - Gamma(A)) * (A_STAR - A)


def dh(A):
    G = Gamma(A)
    return -(1.0 - G) - 2.0 * (A - A_STAR) ** 2 / KAPPA ** 2 * G


def sw(A, w=0.0):
    """Switch (A-A_d)_+ or its softplus smoothing."""
    x = A - A_D
    if w == 0.0:
        return np.maximum(x, 0.0)
    return w * np.logaddexp(0.0, x / w)


def dsw(A, w=0.0):
    x = A - A_D
    if w == 0.0:
        return (x > 0).astype(float) if np.ndim(x) else float(x > 0)
    return 1.0 / (1.0 + np.exp(-x / w))


class Model:
    def __init__(self, child="A", a_E=None, w=0.0, aZ=0.0, A_sigma=A_STAR):
        self.A_b = CHILD[child]["A_b"]
        self.a_E = a_E_of(child) if a_E is None else a_E
        self.w = w
        self.aZ = aZ
        self.A_sigma = A_sigma

    # drift of the planar subsystem
    def f(self, A, R, s):
        return (AH_HAT * s * (1 - A) + self.aZ * (self.A_sigma - A)
                + self.a_E * h(A) + R * (self.A_b - A))

    def g(self, A, R, s):
        return EPS * (1 - R - NUH_HAT * s * R - THETA_A * sw(A, self.w) * R)

    def Rnull(self, A, s):
        return 1.0 / (1.0 + NUH_HAT * s + THETA_A * sw(A, self.w))

    def F(self, A, s):
        return self.f(A, self.Rnull(A, s), s)

    def jac(self, A, R, s):
        fA = -AH_HAT * s - self.aZ - R + self.a_E * dh(A)
        fR = self.A_b - A
        gA = -EPS * THETA_A * dsw(A, self.w) * R
        gR = -EPS * (1 + NUH_HAT * s + THETA_A * sw(A, self.w))
        return np.array([[fA, fR], [gA, gR]])

    def equilibria(self, s, n=4001):
        """All zeros of F on (0,1), sorted."""
        A = np.linspace(1e-6, 1 - 1e-6, n)
        Fv = self.F(A, s)
        roots = []
        for i in np.nonzero(np.sign(Fv[:-1]) * np.sign(Fv[1:]) < 0)[0]:
            roots.append(brentq(self.F, A[i], A[i + 1], args=(s,), xtol=1e-14))
        return np.array(roots)

    def branches(self, s):
        """(regulated, dysregulated) arousal or None."""
        r = self.equilibria(s)
        if len(r) == 3:
            return r[0], r[2]
        if len(r) == 1:
            return (r[0], None) if r[0] < A_D else (None, r[0])
        return None, None

    def lin_var(self, A, s, sig, sigR=0.0):
        """Stationary covariance of the linear-noise approximation."""
        R = self.Rnull(A, s)
        J = self.jac(A, R, s)
        D = np.diag([(sig * A * (1 - A)) ** 2, (sigR * R * (1 - R)) ** 2])
        return solve_continuous_lyapunov(J, -D)

    # fold locations
    def s_on(self):
        """Onset fold. Nonsmooth: F(A_d; s)=0. Smooth: lower fold."""
        if self.w == 0.0:
            return brentq(lambda s: self.F(A_D, s), 1e-6, 1.0, xtol=1e-14)
        return self._smooth_fold(upper=False)

    def s_rec(self):
        return self._smooth_fold(upper=True)

    def _smooth_fold(self, upper):
        """Fold where F = F_A = 0, located by the number of equilibria."""
        def n_eq(s):
            return len(self.equilibria(s))
        ss = np.linspace(0.0, 1.0, 1001)
        ns = np.array([n_eq(s) for s in ss])
        idx = np.nonzero(ns == 3)[0]
        if len(idx) == 0:
            return np.nan
        lo, hi = (ss[idx[0] - 1], ss[idx[0]]) if upper else (ss[idx[-1]], ss[idx[-1] + 1])
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            three = n_eq(mid) == 3
            if upper:
                lo, hi = (lo, mid) if three else (mid, hi)
            else:
                lo, hi = (mid, hi) if three else (lo, mid)
        return 0.5 * (lo + hi)


def _dF(m, A, s):
    R = m.Rnull(A, s)
    J = m.jac(A, R, s)
    dR = -THETA_A * dsw(A, m.w) * R ** 2
    return J[0, 0] + J[0, 1] * dR


def refine_fold(m, s0):
    """Newton refinement of a smooth fold F = F_A = 0 near s0."""
    from scipy.optimize import fsolve
    r = m.equilibria(s0 + 1e-4) if True else None
    # initial A: midpoint of the two colliding roots closest together
    rr = m.equilibria(s0 + 1e-4)
    if len(rr) < 2:
        rr = m.equilibria(s0 - 1e-4)
    k = np.argmin(np.diff(rr))
    A0 = 0.5 * (rr[k] + rr[k + 1])
    sol = fsolve(lambda x: [m.F(x[0], x[1]), _dF(m, x[0], x[1])], [A0, s0],
                 xtol=1e-14)
    return sol[0], sol[1]
