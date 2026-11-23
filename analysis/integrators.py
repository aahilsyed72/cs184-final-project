"""Integrator accuracy study for light bending around a Schwarzschild black hole.

A photon's orbit in its orbital plane, with u = 1/r and units where the
Schwarzschild radius r_s = 1, obeys

    u'' = -u + 1.5 u^2        (' = d/dphi)

A ray coming in from infinity starts at u = 0, u' = 1/b, where b is the impact
parameter. If it escapes, u returns to 0 at some angle phi_exit, and the
deflection angle is phi_exit - pi. Rays with b below the critical value
b_c = 3*sqrt(3)/2 ~= 2.598 are captured.

This script compares three integrators written by hand (Euler, fixed-step RK4,
and adaptive Dormand-Prince RK45) against a ground truth computed from the
exact deflection integral. It produces:

  images/integrator_error.png        error vs. function evaluations, one panel per b
  analysis/results/convergence.csv   raw data behind that plot
  analysis/results/critical_b.md     smallest escaping b per method and setting

Run from the repo root:  analysis/.venv/bin/python analysis/integrators.py
"""

import csv
import math
import os
import time
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import IntegrationWarning, quad, solve_ivp

B_CRIT = 1.5 * math.sqrt(3.0)
PHI_MAX = 200.0  # rays still bound after this much angle count as captured
U_HORIZON = 1.0  # u >= 1 means r <= r_s

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IMAGES = os.path.join(ROOT, "images")
RESULTS = os.path.join(ROOT, "analysis", "results")


# ---------------------------------------------------------------------------
# The ODE, with an evaluation counter. State is (u, w) with w = u'.
# ---------------------------------------------------------------------------

class Rhs:
    def __init__(self):
        self.n = 0

    def __call__(self, u, w):
        self.n += 1
        return w, -u + 1.5 * u * u


# ---------------------------------------------------------------------------
# Single steps. Each returns the new (u, w).
# ---------------------------------------------------------------------------

def euler_step(f, u, w, h):
    du, dw = f(u, w)
    return u + h * du, w + h * dw


def rk4_step(f, u, w, h):
    a1, b1 = f(u, w)
    a2, b2 = f(u + 0.5 * h * a1, w + 0.5 * h * b1)
    a3, b3 = f(u + 0.5 * h * a2, w + 0.5 * h * b2)
    a4, b4 = f(u + h * a3, w + h * b3)
    return (u + h / 6.0 * (a1 + 2 * a2 + 2 * a3 + a4),
            w + h / 6.0 * (b1 + 2 * b2 + 2 * b3 + b4))


# Dormand-Prince 5(4) coefficients.
DP_A = [
    [],
    [1 / 5],
    [3 / 40, 9 / 40],
    [44 / 45, -56 / 15, 32 / 9],
    [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
    [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
    [35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84],
]
DP_B5 = [35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0.0]
DP_B4 = [5179 / 57600, 0.0, 7571 / 16695, 393 / 640, -92097 / 339200,
         187 / 2100, 1 / 40]
DP_E = [b5 - b4 for b5, b4 in zip(DP_B5, DP_B4)]


def dp45_step(f, u, w, h, k1=None):
    """One Dormand-Prince step. Returns (u5, w5, err_u, err_w, k7).

    k1 may be passed in from the previous accepted step (first-same-as-last),
    which saves one evaluation. k7 is f at the new point.
    """
    ks = [k1 if k1 is not None else f(u, w)]
    for i in range(1, 7):
        uu = u + h * sum(a * k[0] for a, k in zip(DP_A[i], ks))
        ww = w + h * sum(a * k[1] for a, k in zip(DP_A[i], ks))
        ks.append(f(uu, ww))
    # The 7th stage is evaluated at the 5th-order solution itself.
    u5 = u + h * sum(b * k[0] for b, k in zip(DP_B5, ks))
    w5 = w + h * sum(b * k[1] for b, k in zip(DP_B5, ks))
    eu = h * sum(e * k[0] for e, k in zip(DP_E, ks))
    ew = h * sum(e * k[1] for e, k in zip(DP_E, ks))
    return u5, w5, eu, ew, ks[6]


def dp45_plain(f, u, w, h):
    """A DP45 step used without error control (for the exit-crossing search)."""
    u5, w5, _, _, _ = dp45_step(f, u, w, h)
    return u5, w5


# ---------------------------------------------------------------------------
# Tracing a ray: integrate until it escapes (u back to 0) or is captured.
# ---------------------------------------------------------------------------

def locate_exit(step, f, u, w, phi, h, u_next):
    """The step from (phi, u, w) of size h ends with u < 0. Find the exact exit
    angle by Newton's method on the step size, using the method's own step,
    so the exit angle has the same order of accuracy as the method.
    """
    s = h * u / (u - u_next)
    for _ in range(3):
        us, ws = step(f, u, w, s)
        if ws == 0.0:
            break
        s -= us / ws
    return phi + s


def trace_fixed(step, b, h):
    """Fixed-step trace. Returns (escaped, phi_exit, n_evals)."""
    f = Rhs()
    u, w, phi = 0.0, 1.0 / b, 0.0
    while phi < PHI_MAX:
        u_next, w_next = step(f, u, w, h)
        if u_next < 0.0:
            return True, locate_exit(step, f, u, w, phi, h, u_next), f.n
        if u_next >= U_HORIZON:
            return False, None, f.n
        u, w, phi = u_next, w_next, phi + h
    return False, None, f.n


def trace_rk45(b, tol, h0=0.05, h_max=0.5):
    """Adaptive Dormand-Prince trace with mixed absolute/relative tolerance."""
    f = Rhs()
    u, w, phi, h = 0.0, 1.0 / b, 0.0, h0
    k1 = None
    while phi < PHI_MAX:
        u5, w5, eu, ew, k7 = dp45_step(f, u, w, h, k1)
        scale_u = tol * (1.0 + max(abs(u), abs(u5)))
        scale_w = tol * (1.0 + max(abs(w), abs(w5)))
        err = max(abs(eu) / scale_u, abs(ew) / scale_w)
        if err <= 1.0:
            if u5 < 0.0:
                return True, locate_exit(dp45_plain, f, u, w, phi, h, u5), f.n
            if u5 >= U_HORIZON:
                return False, None, f.n
            u, w, phi, k1 = u5, w5, phi + h, k7
        else:
            k1 = k1  # rejected: k1 at the old point is still valid
        factor = 5.0 if err == 0.0 else min(5.0, max(0.2, 0.9 * err ** -0.2))
        h = min(h * factor, h_max)
    return False, None, f.n


def trace(method, b, setting):
    if method == "Euler":
        return trace_fixed(euler_step, b, setting)
    if method == "RK4":
        return trace_fixed(rk4_step, b, setting)
    if method == "RK45":
        return trace_rk45(b, setting)
    raise ValueError(method)


# ---------------------------------------------------------------------------
# Ground truth.
# ---------------------------------------------------------------------------

def deflection_exact(b):
    """Exact deflection from the integral phi_total = 2 * int_0^u0 du / sqrt(F),
    F(u) = 1/b^2 - u^2 + u^3 = (u - u0)(u - u1)(u - u2).

    Substituting u = u0 - s^2 removes the square-root singularity at the
    turning point u0, leaving a smooth integrand.
    """
    roots = np.roots([1.0, -1.0, 0.0, 1.0 / b ** 2])
    roots = np.sort(roots.real)
    u2, u0, u1 = roots  # u2 < 0 < u0 < u1
    # Polish the roots with Newton so the near-critical case stays accurate.
    def polish(r):
        for _ in range(3):
            r -= (r ** 3 - r ** 2 + 1 / b ** 2) / (3 * r ** 2 - 2 * r)
        return r
    u0, u1 = polish(u0), polish(u1)

    def integrand(s):
        return 2.0 / math.sqrt((u1 - u0 + s * s) * (u0 - u2 - s * s))

    width = math.sqrt(max(u1 - u0, 1e-30))
    pts = [p for p in (width, 3 * width, 10 * width) if p < math.sqrt(u0)]
    # Asking for 1e-15 makes quad warn about roundoff near the critical b; the
    # answer still agrees with DOP853 to ~4e-12 there (printed in main).
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", IntegrationWarning)
        val, _ = quad(integrand, 0.0, math.sqrt(u0), epsabs=1e-15,
                      epsrel=1e-15, limit=500, points=pts or None)
    return 2.0 * val - math.pi


def deflection_scipy(b):
    """Independent check: scipy DOP853 with very tight tolerances."""
    def rhs(phi, y):
        return [y[1], -y[0] + 1.5 * y[0] ** 2]

    def exit_event(phi, y):
        return y[0]
    exit_event.terminal = True
    exit_event.direction = -1

    sol = solve_ivp(rhs, (0.0, PHI_MAX), [0.0, 1.0 / b], method="DOP853",
                    rtol=1e-13, atol=1e-15, events=exit_event,
                    first_step=1e-3)
    # The first event is phi = 0 itself (u starts at 0, rising), which
    # direction=-1 already excludes.
    return sol.t_events[0][0] - math.pi


# ---------------------------------------------------------------------------
# Experiments.
# ---------------------------------------------------------------------------

B_VALUES = [20.0, 6.0, 3.5, 2.8, 2.65, 2.6]
SWEEPS = {
    "Euler": np.logspace(-1.0, -4.5, 8),   # step size h
    "RK4": np.logspace(-0.3, -2.7, 9),     # step size h
    "RK45": np.logspace(-3.0, -13.0, 11),  # tolerance
}
COLORS = {"Euler": "#d1495b", "RK4": "#2e86ab", "RK45": "#3c9d5d"}


def convergence_study():
    rows = []
    truth = {}
    print("Ground truth deflection (exact integral vs. scipy DOP853):")
    for b in B_VALUES:
        ex = deflection_exact(b)
        sp = deflection_scipy(b)
        truth[b] = ex
        print(f"  b = {b:6.2f}: {ex:.15f} rad   |diff| = {abs(ex - sp):.1e}")

    for b in B_VALUES:
        for method, settings in SWEEPS.items():
            for s in settings:
                escaped, phi_exit, n = trace(method, b, float(s))
                err = abs(phi_exit - math.pi - truth[b]) if escaped else float("nan")
                rows.append({"b": b, "method": method, "setting": float(s),
                             "evals": n, "escaped": escaped, "abs_error": err})
    return rows, truth


def plot_convergence(rows, truth):
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True, sharey=True)
    floor = 1e-16
    for ax, b in zip(axes.flat, B_VALUES):
        for method in SWEEPS:
            pts = [(r["evals"], max(r["abs_error"], floor)) for r in rows
                   if r["b"] == b and r["method"] == method and r["escaped"]]
            if pts:
                x, y = zip(*sorted(pts))
                ax.loglog(x, y, "o-", color=COLORS[method], label=method, ms=4)
            # Settings where the method wrongly captured the ray.
            bad = [r["evals"] for r in rows
                   if r["b"] == b and r["method"] == method and not r["escaped"]]
            for xb in bad:
                ax.loglog(xb, 1.0, "x", color=COLORS[method], ms=8, mew=2)
        ax.set_title(f"b = {b:g}   (deflection {truth[b]:.3f} rad)", fontsize=10)
        ax.grid(True, which="major", alpha=0.3)
        ax.set_ylim(1e-16, 10)
    for ax in axes[1]:
        ax.set_xlabel("function evaluations")
    for ax in axes[:, 0]:
        ax.set_ylabel("|deflection error| (rad)")
    axes[0, 0].legend(loc="lower left")
    fig.suptitle("Deflection error vs. cost  "
                 "(x at error 1 = ray wrongly captured)", fontsize=12)
    fig.tight_layout()
    path = os.path.join(IMAGES, "integrator_error.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


CRIT_SETTINGS = {
    "Euler": [0.1, 0.01, 0.001],
    "RK4": [0.5, 0.1, 0.01],
    "RK45": [1e-3, 1e-6, 1e-9],
}


def critical_b(method, setting, lo=2.0, hi=3.5, iters=36):
    """Bisect for the smallest b that escapes (bracket width ~1e-11)."""
    assert not trace(method, lo, setting)[0] and trace(method, hi, setting)[0]
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if trace(method, mid, setting)[0]:
            hi = mid
        else:
            lo = mid
    return hi


def critical_study():
    rows = []
    for method, settings in CRIT_SETTINGS.items():
        for s in settings:
            bc = critical_b(method, s)
            rows.append({"method": method, "setting": s, "b_crit": bc,
                         "error": bc - B_CRIT})
            print(f"  {method:5s} {s:<8g} b_c = {bc:.10f}   error = {bc - B_CRIT:+.2e}")
    return rows


def write_outputs(conv_rows, crit_rows):
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "convergence.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(conv_rows[0]))
        wr.writeheader()
        wr.writerows(conv_rows)

    lines = [
        f"Exact critical impact parameter: 3*sqrt(3)/2 = {B_CRIT:.10f}",
        "",
        "| Method | Setting | Smallest escaping b | Error |",
        "|---|---|---|---|",
    ]
    for r in crit_rows:
        label = f"tol {r['setting']:g}" if r["method"] == "RK45" else f"h = {r['setting']:g}"
        lines.append(f"| {r['method']} | {label} | {r['b_crit']:.10f} | {r['error']:+.2e} |")
    with open(os.path.join(RESULTS, "critical_b.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    t0 = time.time()
    conv_rows, truth = convergence_study()
    path = plot_convergence(conv_rows, truth)
    print(f"Wrote {os.path.relpath(path, ROOT)}")
    print("Critical impact parameter (exact 2.5980762114):")
    crit_rows = critical_study()
    write_outputs(conv_rows, crit_rows)
    print(f"Wrote analysis/results/  ({time.time() - t0:.0f} s total)")


if __name__ == "__main__":
    main()
