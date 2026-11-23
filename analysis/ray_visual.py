"""Picture of light rays bending around the black hole: exact vs. Euler.

The same set of incoming rays is traced twice, with a near-exact integrator
(RK4, tiny step) and with Euler at a coarse step. Euler adds energy to the orbit,
so rays that should skim past the hole get pulled in, and the shadow comes out
too big.

Outputs:
  images/ray_paths.png        still image (both panels)
  analysis/results/ray_paths.mp4   ~6 s animation of the rays flying in, for the video

Run from the repo root:  analysis/.venv/bin/python analysis/ray_visual.py
"""

import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from integrators import B_CRIT, Rhs, euler_step, rk4_step  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VIEW = 8.0         # half-width of the plotted window, in Schwarzschild radii
R_START = 11.0     # only keep the part of each ray inside this radius
EULER_H = 0.1
EXACT_H = 0.002

B_LIST = [1.0, 2.0, 2.45, 2.63, 2.67, 3.0, 4.0, 5.5, 7.0]

BG = "#05060a"
COL_OK = "#7fb3ff"
COL_CAPTURED = "#ff6b6b"
COL_WRONG = "#ffd166"


def trace_path(step, b, h):
    """Return (x, y, captured) for the part of the ray inside R_START.
    The ray comes in from the left (x -> -inf) at height b."""
    f = Rhs()
    u, w, phi = 0.0, 1.0 / b, 0.0
    us, phis = [], []
    captured = False
    while phi < 60.0:
        u_next, w_next = step(f, u, w, h)
        phi += h
        if u_next < 0.0:
            break
        if u_next >= 1.0:
            captured = True
            us.append(1.0)
            phis.append(phi)
            break
        u, w = u_next, w_next
        if u > 1.0 / R_START:
            us.append(u)
            phis.append(phi)
    r = 1.0 / np.array(us)
    phis = np.array(phis)
    return -r * np.cos(phis), r * np.sin(phis), captured


def all_paths(step, h):
    paths = []
    for b in B_LIST:
        x, y, cap = trace_path(step, b, h)
        paths.append({"b": b, "x": x, "y": y, "captured": cap})
    return paths


def arclength(p):
    d = np.hypot(np.diff(p["x"]), np.diff(p["y"]))
    return np.concatenate([[0.0], np.cumsum(d)])


def setup_axes(ax, title):
    ax.set_facecolor(BG)
    ax.set_xlim(-VIEW, VIEW)
    ax.set_ylim(-VIEW * 0.8, VIEW * 0.8)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.add_patch(plt.Circle((0, 0), 1.5, fill=False, ls="--", lw=1,
                            color="#888888", alpha=0.6))
    ax.add_patch(plt.Circle((0, 0), 1.0, color="black", zorder=5))
    ax.add_patch(plt.Circle((0, 0), 1.0, fill=False, lw=1.2,
                            color="#444444", zorder=6))
    ax.set_title(title, color="white", fontsize=20, pad=12)


def colors(exact, euler):
    """Rays where Euler and the exact answer disagree are yellow in both panels."""
    out_exact, out_euler = [], []
    for pe, pu in zip(exact, euler):
        if pu["captured"] != pe["captured"]:
            out_exact.append(COL_WRONG)
            out_euler.append(COL_WRONG)
        else:
            out_exact.append(COL_CAPTURED if pe["captured"] else COL_OK)
            out_euler.append(COL_CAPTURED if pu["captured"] else COL_OK)
    return out_exact, out_euler


def make_figure():
    fig, axes = plt.subplots(1, 2, figsize=(19.2, 10.8), facecolor=BG)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.85, bottom=0.11, wspace=0.03)
    setup_axes(axes[0], "Accurate (RK4, tiny steps)")
    setup_axes(axes[1], f"Euler (step {EULER_H})")
    fig.suptitle("Light rays passing a black hole", color="white",
                 fontsize=28, y=0.965)
    fig.text(0.5, 0.045,
             "blue = escapes     red = falls in     "
             "yellow = the same rays: they really escape (left), but Euler pulls them in (right)",
             color="#cccccc", ha="center", fontsize=17)
    return fig, axes


def legend_extra(ax):
    pass


def still(exact, euler):
    cex, ceu = colors(exact, euler)
    fig, axes = make_figure()
    for ax, paths, cols in ((axes[0], exact, cex), (axes[1], euler, ceu)):
        for p, c in zip(paths, cols):
            lw = 4.0 if c == COL_WRONG else 2.2
            ax.plot(p["x"], p["y"], color=c, lw=lw, zorder=4 if c == COL_WRONG else 3)
        legend_extra(ax)
    path = os.path.join(ROOT, "images", "ray_paths.png")
    fig.savefig(path, dpi=100, facecolor=BG)
    plt.close(fig)
    return path


def animation(exact, euler, fps=30, grow_s=5.0, hold_s=2.0):
    cex, ceu = colors(exact, euler)
    fig, axes = make_figure()
    lines = []
    for ax, paths, cols in ((axes[0], exact, cex), (axes[1], euler, ceu)):
        legend_extra(ax)
        for p, c in zip(paths, cols):
            lw = 4.0 if c == COL_WRONG else 2.2
            (ln,) = ax.plot([], [], color=c, lw=lw,
                            zorder=4 if c == COL_WRONG else 3)
            lines.append((ln, p, arclength(p)))
    # Every ray moves at the same speed along its own path.
    max_len = max(s[-1] for _, _, s in lines)
    n_grow, n_hold = int(fps * grow_s), int(fps * hold_s)

    def update(i):
        L = max_len * min(1.0, (i + 1) / n_grow)
        for ln, p, s in lines:
            k = np.searchsorted(s, L, side="right")
            ln.set_data(p["x"][:k], p["y"][:k])
        return [ln for ln, _, _ in lines]

    anim = FuncAnimation(fig, update, frames=n_grow + n_hold, blit=True)
    out = os.path.join(ROOT, "analysis", "results", "ray_paths.mp4")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    anim.save(out, writer=FFMpegWriter(fps=fps, bitrate=8000,
                                       extra_args=["-pix_fmt", "yuv420p"]),
              dpi=100, savefig_kwargs={"facecolor": BG})
    plt.close(fig)
    return out


def main():
    exact = all_paths(rk4_step, EXACT_H)
    euler = all_paths(euler_step, EULER_H)
    for pe, pu in zip(exact, euler):
        flag = "  <-- Euler wrong" if pe["captured"] != pu["captured"] else ""
        print(f"b = {pe['b']:5.2f}  exact: {'in ' if pe['captured'] else 'out'}"
              f"  euler: {'in ' if pu['captured'] else 'out'}{flag}")
    print(f"(true cutoff b = {B_CRIT:.3f})")
    print("Wrote", os.path.relpath(still(exact, euler), ROOT))
    print("Wrote", os.path.relpath(animation(exact, euler), ROOT))


if __name__ == "__main__":
    main()
