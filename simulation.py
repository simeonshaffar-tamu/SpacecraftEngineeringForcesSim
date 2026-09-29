"""Glue: config -> spacecraft state and perturbing accelerations at a given true anomaly."""

from dataclasses import replace
from typing import Tuple

import numpy as np

import constants as C
from config import OrbitElements, SimConfig
from orbital_elements import cartesian_to_keplerian, keplerian_to_cartesian
from perturbations import PerturbationResult, compute_perturbations
from typing_aliases import Vec3


def accelerations_at(cfg: SimConfig, r: Vec3, v: Vec3) -> PerturbationResult:
    """Perturbing accelerations [m/s^2] for the ECI state (r [m], v [m/s]) using cfg's switches/parameters."""
    sw, sc, env = cfg.switches, cfg.spacecraft, cfg.environment
    return compute_perturbations(
        r, v, sc.mass_kg, sc.drag_area_m2, sc.cd, env.rho_kg_m3,
        np.array(env.r_sun_m), np.array(env.r_moon_m),
        enable_drag=sw.drag, enable_j2=sw.j2, enable_moon=sw.moon, enable_sun=sw.sun)


def initial_state(cfg: SimConfig) -> Tuple[Vec3, Vec3]:
    """ECI (r [m], v [m/s]) at the epoch, from Keplerian elements or the Cartesian input per cfg.input_mode."""
    if cfg.input_mode == "cartesian":
        return np.array(cfg.cartesian.r_m, dtype=float), np.array(cfg.cartesian.v_mps, dtype=float)
    if cfg.input_mode == "keplerian":
        return compute_state(cfg, cfg.orbit.nu_deg)[:2]
    raise ValueError(f"input_mode must be 'keplerian' or 'cartesian', got {cfg.input_mode!r}")


def with_derived_orbit(cfg: SimConfig) -> SimConfig:
    """Return a config whose `orbit` matches the initial state (derived from r, v in cartesian mode).

    Raises ValueError if the Cartesian state is not on a bound, non-degenerate orbit.
    """
    if cfg.input_mode != "cartesian":
        return cfg
    r, v = initial_state(cfg)
    a, e, i, raan, argp, nu = cartesian_to_keplerian(r, v, C.MU_EARTH)
    d = np.degrees
    return replace(cfg, orbit=OrbitElements(a_m=a, e=e, i_deg=d(i), raan_deg=d(raan),
                                            argp_deg=d(argp), nu_deg=d(nu)))


def compute_state(cfg: SimConfig, nu_deg: float) -> Tuple[Vec3, Vec3, PerturbationResult]:
    """Return (r [m], v [m/s], accelerations [m/s^2]) on the orbit in cfg.orbit at true anomaly nu_deg [deg]."""
    o = cfg.orbit
    r, v = keplerian_to_cartesian(
        o.a_m, o.e, np.radians(o.i_deg), np.radians(o.raan_deg),
        np.radians(o.argp_deg), np.radians(nu_deg), C.MU_EARTH)
    return r, v, accelerations_at(cfg, r, v)
