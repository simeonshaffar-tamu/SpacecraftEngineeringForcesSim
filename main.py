"""Compute instantaneous perturbing accelerations at the epoch defined in config.py."""

import numpy as np

import constants as C
from config import CONFIG
from orbital_elements import keplerian_to_cartesian
from perturbations import compute_perturbations


def fmt(vec):
    return "[" + ", ".join(f"{x: .6e}" for x in vec) + "]"


def main():
    cfg = CONFIG
    o, sw, sc, env = cfg.orbit, cfg.switches, cfg.spacecraft, cfg.environment

    r, v = keplerian_to_cartesian(
        o.a_m, o.e, np.radians(o.i_deg), np.radians(o.raan_deg),
        np.radians(o.argp_deg), np.radians(o.nu_deg), C.MU_EARTH)

    res = compute_perturbations(
        r, v, sc.mass_kg, sc.drag_area_m2, sc.cd, env.rho_kg_m3,
        env.r_sun_m, env.r_moon_m,
        enable_drag=sw.drag, enable_j2=sw.j2, enable_moon=sw.moon, enable_sun=sw.sun)

    print("=== Initial state (ECI) ===")
    print(f"r [m]   = {fmt(r)}   |r| = {np.linalg.norm(r):.3f}")
    print(f"v [m/s] = {fmt(v)}   |v| = {np.linalg.norm(v):.3f}")
    print(f"\n=== Perturbing accelerations [m/s^2]  (switches: {sw}) ===")
    for name, a in (("Drag", res.drag), ("J2", res.j2), ("Moon", res.moon), ("Sun", res.sun)):
        print(f"{name:<6}= {fmt(a)}   |a| = {np.linalg.norm(a):.6e}")
    print(f"{'Total':<6}= {fmt(res.total)}   |a| = {np.linalg.norm(res.total):.6e}")


if __name__ == "__main__":
    main()
