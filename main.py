"""Compute instantaneous perturbing accelerations at the epoch defined in config.py.

Usage: python main.py            # print results, save PNGs, open interactive figure
       python main.py --no-show  # print results and save PNGs only
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np

from config import CONFIG
from plotting import make_interactive_figure, save_separate_figures
from simulation import accelerations_at, initial_state, with_derived_orbit
from typing_aliases import Vec3


def fmt(vec: Vec3) -> str:
    return "[" + ", ".join(f"{x: .6e}" for x in vec) + "]"


def main(show: bool = True) -> None:
    cfg = CONFIG
    r, v = initial_state(cfg)
    res = accelerations_at(cfg, r, v)

    print(f"=== Initial state (ECI), input mode: {cfg.input_mode} ===")
    print(f"r [m]   = {fmt(r)}   |r| = {np.linalg.norm(r):.3f}")
    print(f"v [m/s] = {fmt(v)}   |v| = {np.linalg.norm(v):.3f}")
    print(f"\n=== Perturbing accelerations [m/s^2]  (switches: {cfg.switches}) ===")
    for name, a in (("Drag", res.drag), ("J2", res.j2), ("Moon", res.moon), ("Sun", res.sun)):
        print(f"{name:<6}= {fmt(a)}   |a| = {np.linalg.norm(a):.6e}")
    print(f"{'Total':<6}= {fmt(res.total)}   |a| = {np.linalg.norm(res.total):.6e}")

    try:
        plot_cfg = with_derived_orbit(cfg)          # in cartesian mode, orbit elements come from r, v
    except ValueError as err:
        print(f"\nSkipping plots: {err}")
        return
    if cfg.input_mode == "cartesian":
        o = plot_cfg.orbit
        print(f"\nDerived elements: a = {o.a_m:.3f} m, e = {o.e:.6f}, i = {o.i_deg:.4f} deg, "
              f"RAAN = {o.raan_deg:.4f} deg, argp = {o.argp_deg:.4f} deg, nu = {o.nu_deg:.4f} deg")

    save_separate_figures(plot_cfg)
    if show:
        fig = make_interactive_figure(plot_cfg)
        fig.savefig("outputs/combined_view.png", dpi=110)
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-show", action="store_true", help="save PNGs without opening the interactive figure")
    main(show=not parser.parse_args().no_show)
