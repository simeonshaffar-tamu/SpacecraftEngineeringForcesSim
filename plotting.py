"""Diagrams of the orbit state and Sun/Moon geometry.

Each view is a class that draws onto a supplied Axes, so it can be saved standalone or shown
together with the other view (plus a true-anomaly slider) on one interactive figure.
"""

import os
from typing import Callable, Dict, List

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from matplotlib.widgets import Slider

import constants as C
from config import SimConfig
from orbital_elements import perifocal_to_eci
from perturbations import PerturbationResult
from simulation import compute_state
from typing_aliases import Vec2, Vec3

ACC_LEGEND_TITLE = ("Acceleration arrows show direction only (Sun/Moon are differential, i.e. tidal, "
                    "so they do not point at the body)")

# Acceleration arrow colours (prominent); r and v are deliberately muted
ACC_COLORS: Dict[str, str] = {"Drag": "#d62728", "J2": "#2ca02c", "Moon": "#17becf", "Sun": "#ff9900"}


class AccelArrows:
    """The four acceleration arrows, anchored at the spacecraft, direction only.

    Each arrow has fixed nominal length `length` times the fraction of its 3-D magnitude that lies in
    the view plane (given by a projection function), so mostly out-of-plane vectors look short.
    """

    def __init__(self, ax: Axes, length: float) -> None:
        self.length = length
        arrow = dict(angles="xy", scale_units="xy", scale=1, width=0.004)
        self.quivers = {n: ax.quiver(0, 0, 0, 0, color=c, zorder=4, **arrow) for n, c in ACC_COLORS.items()}

    def update(self, anchor: Vec2, res: PerturbationResult,
               proj: Callable[[Vec3], Vec2], alpha: float = 1.0) -> List[str]:
        """Redraw arrows at `anchor` (alpha < 1 to fade them); return legend label strings in Drag, J2, Moon, Sun order."""
        labels = []
        for name, acc in (("Drag", res.drag), ("J2", res.j2), ("Moon", res.moon), ("Sun", res.sun)):
            mag = np.linalg.norm(acc)
            q = self.quivers[name]
            q.set_offsets([anchor])
            q.set_alpha(alpha)
            if mag > 0.0:
                u = proj(acc) / mag * self.length
                q.set_UVC([u[0]], [u[1]])
                labels.append(f"{name} accel (|a| = {mag:.2e} m/s²)")
            else:
                q.set_UVC([0.0], [0.0])
                labels.append(f"{name} accel (off)")
        return labels


class OrbitPlaneView:
    """View down the orbit's angular momentum vector (h out of page).

    x = perigee direction, y = 90 deg ahead of perigee. Vectors are projected onto this plane.
    Acceleration arrows show direction only: each has fixed nominal length times the fraction of its
    3-D magnitude lying in this plane (magnitudes are listed in the legend). r and v are muted.
    """

    def __init__(self, ax: Axes, cfg: SimConfig) -> None:
        self.ax, self.cfg = ax, cfg
        o = cfg.orbit
        a, e = o.a_m, o.e
        R = perifocal_to_eci(np.radians(o.i_deg), np.radians(o.raan_deg), np.radians(o.argp_deg))
        self.P, self.Q = R[:, 0], R[:, 1]
        self.v_scale_len = 0.25 * a      # nominal velocity arrow length
        self.acc_len = 0.20 * a          # nominal acceleration arrow length

        nu = np.linspace(0.0, 2.0 * np.pi, 720)
        rad = a * (1.0 - e**2) / (1.0 + e * np.cos(nu))

        ax.add_patch(Circle((0, 0), C.R_EARTH, color="#7fb3d5", alpha=0.6, zorder=1))
        h_earth = Circle((0, 0), 0, color="#7fb3d5", alpha=0.6, label=f"Earth (R = {C.R_EARTH/1e3:.1f} km)")
        h_orbit, = ax.plot(rad * np.cos(nu), rad * np.sin(nu), color="tab:blue", lw=1.5, zorder=2, label="Orbit")
        h_peri, = ax.plot([0, a * (1.0 - e)], [0, 0], color="gray", ls=":", lw=1.8, zorder=2,
                          label=f"Perigee position (r_p = {a*(1.0-e)/1e3:.1f} km)")
        ax.plot(0, 0, "+", color="k", ms=8, zorder=2)

        arrow = dict(angles="xy", scale_units="xy", scale=1, width=0.004)
        self.q_r = ax.quiver(0, 0, 0, 0, color="#bdbdbd", zorder=2, **arrow)
        self.q_v = ax.quiver(0, 0, 0, 0, color="#cdbfe3", zorder=2, **arrow)
        self.sc, = ax.plot([], [], "o", color="k", ms=7, zorder=6, label="Spacecraft")
        self.acc = AccelArrows(ax, self.acc_len)

        self.handles = [h_earth, h_orbit, h_peri, self.q_r, self.q_v, self.sc, *self.acc.quivers.values()]
        self.legend = ax.legend(handles=self.handles, labels=[h.get_label() for h in self.handles[:3]] + [""] * 7,
                                loc="upper center", bbox_to_anchor=(0.5, -0.1), fontsize=8, ncol=2,
                                title=ACC_LEGEND_TITLE, title_fontsize=8)
        self.labels = self.legend.get_texts()
        self.labels[5].set_text("Spacecraft")

        ax.set_aspect("equal")
        lim = 1.2 * a * (1.0 + e)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("Perigee direction [m]")
        ax.set_ylabel("In-plane, 90° ahead of perigee [m]")
        ax.set_title("Orbit plane — view along angular momentum vector (h out of page)", fontsize=11)
        ax.grid(alpha=0.3)
        self.update(o.nu_deg)

    def _proj(self, vec: Vec3) -> Vec2:
        return np.array([vec @ self.P, vec @ self.Q])

    def update(self, nu_deg: float) -> None:
        """Move the spacecraft to true anomaly nu_deg [deg] and redraw state and acceleration vectors."""
        r, v, res = compute_state(self.cfg, nu_deg)
        r2, v2 = self._proj(r), self._proj(v)
        self.sc.set_data([r2[0]], [r2[1]])
        self.q_r.set_UVC([r2[0]], [r2[1]])
        self.q_v.set_offsets([r2])
        self.q_v.set_UVC([v2[0] * self.v_scale_len / np.linalg.norm(v)], [v2[1] * self.v_scale_len / np.linalg.norm(v)])
        self.labels[3].set_text(f"Position r (|r| = {np.linalg.norm(r)/1e3:.1f} km)")
        self.labels[4].set_text(f"Velocity v (|v| = {np.linalg.norm(v)/1e3:.3f} km/s; arrow length fixed)")

        for k, text in enumerate(self.acc.update(r2, res, self._proj)):
            self.labels[6 + k].set_text(text)


class SunMoonView:
    """View down the angular momentum vector of Earth's orbit about the Sun (ecliptic north).

    Ecliptic normal n = [0, -sin(eps), cos(eps)] in the idealized ECI frame; view x-axis is +x_ECI
    (vernal equinox) and y = n x x. Everything is projected onto the ecliptic plane.

    Occlusion: the viewer is at ecliptic north. A point of the orbit is hidden only if it is behind Earth
    (farther from the viewer than Earth's centre) AND its projection lies inside Earth's disk. Hidden parts are
    thin/dashed and drawn under the translucent disk; everything else is solid, bold and drawn on top. The
    spacecraft is hollow (with faded arrows) when hidden. The equator is a single dotted ellipse; the spin axis
    shows N (toward viewer) / S (behind Earth).

    Earth (translucent), the orbit, the nodes and the spacecraft are drawn TO SCALE relative to each other.
    The Moon (60 Earth radii away) and Sun (~23,000 Earth radii away) cannot fit at that scale, so they are
    shown as markers along their true projected directions, pinned just outside the orbit (marker sizes and
    radial positions are NOT to scale; true distances are in the legend).
    Acceleration arrows are direction-only, projected onto the ecliptic plane, as in OrbitPlaneView.
    Filled circle = ascending node (RAAN), open circle = descending node (LDDN).
    """

    def __init__(self, ax: Axes, cfg: SimConfig) -> None:
        o, env = cfg.orbit, cfg.environment
        a, e = o.a_m, o.e
        eps = np.radians(C.OBLIQUITY_DEG)
        n = np.array([0.0, -np.sin(eps), np.cos(eps)])
        ex = np.array([1.0, 0.0, 0.0])
        ey = np.cross(n, ex)

        def proj(vec: Vec3) -> Vec2:
            return np.array([vec @ ex, vec @ ey])

        def unit(vec: Vec2) -> Vec2:
            return vec / np.linalg.norm(vec)

        self.cfg, self.proj, self.n = cfg, proj, n
        r_max = a * (1.0 + e)
        lim = 1.7 * r_max
        r_sun_marker, r_moon_marker = 1.5 * r_max, 1.25 * r_max     # pinned marker radii (not to scale)
        sun_dir = unit(proj(np.array(env.r_sun_m)))
        moon_dir = unit(proj(np.array(env.r_moon_m)))
        sun_km = np.linalg.norm(env.r_sun_m) / 1e3
        moon_km = np.linalg.norm(env.r_moon_m) / 1e3

        Rm = perifocal_to_eci(np.radians(o.i_deg), np.radians(o.raan_deg), np.radians(o.argp_deg))
        argp = np.radians(o.argp_deg)
        nu = np.linspace(0.0, 2.0 * np.pi, 1440)
        rad = a * (1.0 - e**2) / (1.0 + e * np.cos(nu))
        orb3 = np.array([Rm @ np.array([r_ * np.cos(t), r_ * np.sin(t), 0.0]) for r_, t in zip(rad, nu)])
        orb = np.array([proj(p3) for p3 in orb3])
        orb_depth = orb3 @ n                                 # >0: toward viewer (in front of Earth)

        # Earth's equator (circle in the z_ECI-normal plane) and spin axis
        th = np.linspace(0.0, 2.0 * np.pi, 720)
        eq3 = C.R_EARTH * np.stack([np.cos(th), np.sin(th), np.zeros_like(th)], axis=1)
        eq = np.array([proj(p3) for p3 in eq3])
        z_hat = np.array([0.0, 0.0, 1.0])
        pole_n = C.R_EARTH * proj(z_hat)                    # north pole (depth = R cos(eps) > 0: in front)

        nodes = {}
        for label, nu_n in (("asc", -argp), ("desc", np.pi - argp)):
            r_n = a * (1.0 - e**2) / (1.0 + e * np.cos(nu_n))
            nodes[label] = proj(Rm @ (r_n * np.array([np.cos(nu_n), np.sin(nu_n), 0.0])))

        def split(xy: np.ndarray, hidden: np.ndarray):
            """Return (visible, hidden) copies of xy with NaN elsewhere; segments overlap by one sample."""
            v = ~hidden
            ve = v | np.roll(v, 1) | np.roll(v, -1)
            he = hidden | np.roll(hidden, 1) | np.roll(hidden, -1)
            return np.where(ve[:, None], xy, np.nan), np.where(he[:, None], xy, np.nan)

        earth_col, orb_col, eq_col = "#4a90c2", "#08306b", "#b03a2e"
        # z-order plan: back things 1 < Earth disk 2 < front things 4+
        ax.add_patch(Circle((0, 0), C.R_EARTH, color=earth_col, alpha=0.45, zorder=2))
        ax.add_patch(Circle((0, 0), 0, color=earth_col, alpha=0.45,
                            label=f"Earth (to scale, R = {C.R_EARTH/1e3:.1f} km)"))

        orb_hidden = (orb_depth < 0.0) & (np.linalg.norm(orb, axis=1) < C.R_EARTH)
        orb_f, orb_b = split(orb, orb_hidden)
        ax.plot(*orb_f.T, color=orb_col, lw=2.2, zorder=4, label="Orbit, visible")
        ax.plot(*orb_b.T, color=orb_col, lw=1.0, ls=(0, (4, 3)), alpha=0.7, zorder=1,
                label="Orbit, hidden behind Earth")

        ax.plot(*eq.T, color=eq_col, lw=1.8, ls=":", zorder=3, label="Earth's equator")

        # Spin axis: north half toward the viewer (solid, over Earth), south half behind (dashed, under Earth)
        ax.plot([0, 1.35 * pole_n[0]], [0, 1.35 * pole_n[1]], color="k", lw=1.6, zorder=4)
        ax.plot([0, -1.35 * pole_n[0]], [0, -1.35 * pole_n[1]], color="k", lw=1.0, ls=(0, (4, 3)), zorder=1)
        ax.plot(*pole_n, "o", color="k", ms=5, zorder=5)
        ax.plot(*(-pole_n), "o", mfc="none", mec="k", ms=5, zorder=1)
        ax.text(*(1.55 * pole_n), "N", ha="center", va="center", fontsize=12, fontweight="bold", zorder=6)
        ax.text(*(-1.55 * pole_n), "S", ha="center", va="center", fontsize=12, color="0.35", zorder=6)
        ax.plot([], [], color="k", lw=1.6, label="Spin axis: N (solid) toward viewer, S (dashed) behind")

        ax.plot(*nodes["asc"], "o", color="k", ms=9, zorder=5, label="RAAN (ascending node)")
        ax.plot(*nodes["desc"], "o", mfc="white", mec="k", mew=1.8, ms=9, zorder=5, label="LDDN (descending node)")

        for direction, radius, color, ms, name, km in (
                (sun_dir, r_sun_marker, "orange", 24, "Sun", sun_km),
                (moon_dir, r_moon_marker, "gray", 12, "Moon", moon_km)):
            ax.plot([0, direction[0] * radius], [0, direction[1] * radius], color=color, ls="--", lw=0.8, alpha=0.6, zorder=2)
            ax.plot(*(direction * radius), "o", color=color, ms=ms, zorder=3,
                    label=f"{name} direction ({km:.3g} km away; not to scale)")

        self.sc, = ax.plot([], [], "o", mfc="k", mec="k", ms=6, zorder=6,
                           label="Spacecraft (hollow = hidden by Earth)")
        self.acc = AccelArrows(ax, 0.20 * a)
        ax.annotate("", xy=(0.6, 0.05), xytext=(0.4, 0.05), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="->"))
        ax.annotate("♈ vernal equinox direction (+x ECI)", (0.5, 0.05), xycoords="axes fraction",
                    xytext=(0, -14), textcoords="offset points", ha="center", fontsize=8)

        ax.set_aspect("equal")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("Ecliptic x, vernal equinox direction [m]")
        ax.set_ylabel("Ecliptic y [m]")
        ax.set_title("Orbit relative to Earth, Sun and Moon — view along Earth's orbital angular momentum\n"
                     "(ecliptic north out of page; Earth/orbit to scale, Sun/Moon directions only)", fontsize=10)
        ax.grid(alpha=0.3)
        handles, labels = ax.get_legend_handles_labels()
        self.legend = ax.legend(handles=handles + list(self.acc.quivers.values()),
                                labels=labels + [""] * 4, loc="upper center", bbox_to_anchor=(0.5, -0.1),
                                fontsize=8, ncol=2, markerscale=0.5, labelspacing=0.8,
                                title=ACC_LEGEND_TITLE, title_fontsize=8)
        self.labels = self.legend.get_texts()
        self.n_static = len(handles)
        self.update(o.nu_deg)

    def update(self, nu_deg: float) -> None:
        """Move the spacecraft to true anomaly nu_deg [deg] and redraw the acceleration arrows."""
        r, _, res = compute_state(self.cfg, nu_deg)
        p = self.proj(r)
        behind = bool((r @ self.n) < 0.0 and np.linalg.norm(p) < C.R_EARTH)   # hidden by Earth
        self.sc.set_data([p[0]], [p[1]])
        self.sc.set_markerfacecolor("white" if behind else "k")
        self.sc.set_zorder(1.5 if behind else 6)             # behind: drawn under the translucent Earth
        for k, text in enumerate(self.acc.update(p, res, self.proj, alpha=0.4 if behind else 1.0)):
            self.labels[self.n_static + k].set_text(text)


def _save(fig: Figure, path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")


def save_separate_figures(cfg: SimConfig, out_dir: str = "outputs") -> None:
    """Save each view as its own PNG at the epoch true anomaly in cfg (no slider)."""
    for name, View in (("orbit_plane", OrbitPlaneView), ("sun_moon_view", SunMoonView)):
        fig, ax = plt.subplots(figsize=(8, 8))
        View(ax, cfg)
        _save(fig, os.path.join(out_dir, f"{name}.png"))
        plt.close(fig)


def make_interactive_figure(cfg: SimConfig) -> Figure:
    """Both views side by side with a true-anomaly slider (0-360 deg) driving the orbit-plane view."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 9))
    fig.subplots_adjust(bottom=0.34, wspace=0.15, left=0.06, right=0.97, top=0.93)
    orbit_view = OrbitPlaneView(ax1, cfg)
    sun_moon_view = SunMoonView(ax2, cfg)

    ax_slider = fig.add_axes([0.25, 0.04, 0.5, 0.03])
    slider = Slider(ax_slider, "True anomaly [deg]", 0.0, 360.0, valinit=cfg.orbit.nu_deg % 360.0)

    def on_change(val: float) -> None:
        orbit_view.update(val)
        sun_moon_view.update(val)
        fig.canvas.draw_idle()

    slider.on_changed(on_change)
    fig._slider = slider          # keep a reference so the widget isn't garbage-collected
    fig._orbit_view = orbit_view
    return fig
