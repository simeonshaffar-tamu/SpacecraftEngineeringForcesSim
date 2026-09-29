"""Instantaneous perturbing accelerations in an idealized ECI frame (z = spin axis).
All inputs/outputs SI. Earth's central point-mass term is intentionally excluded."""

from dataclasses import dataclass

import numpy as np

import constants as C
from typing_aliases import Vec3


@dataclass
class PerturbationResult:
    drag: Vec3   # [m/s^2]
    j2: Vec3     # [m/s^2]
    moon: Vec3   # [m/s^2]
    sun: Vec3    # [m/s^2]

    @property
    def total(self) -> Vec3:
        return self.drag + self.j2 + self.moon + self.sun


def _vec3(x: Vec3, name: str) -> Vec3:
    v = np.asarray(x, dtype=float)
    if v.shape != (3,) or not np.all(np.isfinite(v)):
        raise ValueError(f"{name} must be a finite 3-vector")
    return v


def j2_acceleration(r: Vec3, mu: float = C.MU_EARTH, j2: float = C.J2_EARTH,
                    r_eq: float = C.R_EARTH) -> Vec3:
    """Earth J2 acceleration [m/s^2] at ECI position r [m]."""
    x, y, z = r
    r2 = r @ r
    rn = np.sqrt(r2)
    k = -1.5 * j2 * mu * r_eq**2 / rn**5
    zr2 = 5.0 * z**2 / r2
    return k * np.array([x * (1.0 - zr2), y * (1.0 - zr2), z * (3.0 - zr2)])


def drag_acceleration(r: Vec3, v: Vec3, rho: float, cd: float, area: float, mass: float,
                      omega: Vec3 = C.OMEGA_EARTH) -> Vec3:
    """Drag acceleration [m/s^2] against a rigidly co-rotating atmosphere."""
    v_rel = v - np.cross(np.asarray(omega, dtype=float), r)
    return -0.5 * rho * cd * (area / mass) * np.linalg.norm(v_rel) * v_rel


def third_body_acceleration(r: Vec3, r_body: Vec3, mu_body: float) -> Vec3:
    """Differential point-mass acceleration [m/s^2] of the spacecraft relative to Earth."""
    d = r_body - r
    return mu_body * (d / np.linalg.norm(d) ** 3 - r_body / np.linalg.norm(r_body) ** 3)


def compute_perturbations(r: Vec3, v: Vec3, mass: float, drag_area: float, cd: float, rho: float,
                          r_sun: Vec3, r_moon: Vec3,
                          enable_drag: bool = True, enable_j2: bool = True,
                          enable_moon: bool = True, enable_sun: bool = True) -> PerturbationResult:
    """Compute each enabled perturbing acceleration (disabled ones are zero vectors).

    Vec3 args: r [m], v [m/s], r_sun/r_moon (Earth->body) [m]. Scalars: mass [kg], drag_area [m^2],
    cd [-], rho [kg/m^3].
    """
    r, v = _vec3(r, "r"), _vec3(v, "v")
    r_sun, r_moon = _vec3(r_sun, "r_sun"), _vec3(r_moon, "r_moon")
    if mass <= 0.0:
        raise ValueError("mass must be positive")
    if drag_area < 0.0 or cd < 0.0 or rho < 0.0:
        raise ValueError("drag_area, cd and rho must be non-negative")
    if np.linalg.norm(r) == 0.0:
        raise ValueError("spacecraft position must be nonzero")

    zero = np.zeros(3)
    return PerturbationResult(
        drag=drag_acceleration(r, v, rho, cd, drag_area, mass) if enable_drag else zero,
        j2=j2_acceleration(r) if enable_j2 else zero,
        moon=third_body_acceleration(r, r_moon, C.MU_MOON) if enable_moon else zero,
        sun=third_body_acceleration(r, r_sun, C.MU_SUN) if enable_sun else zero,
    )
