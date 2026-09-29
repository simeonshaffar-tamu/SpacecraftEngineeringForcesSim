"""Keplerian elements -> ECI Cartesian state (two-body, no propagation)."""

import numpy as np
from typing import Tuple

from typing_aliases import Mat3, Vec3


def _rot_z(t: float) -> Mat3:
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _rot_x(t: float) -> Mat3:
    c, s = np.cos(t), np.sin(t)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def perifocal_to_eci(i: float, raan: float, argp: float) -> Mat3:
    """Rotation matrix PQW -> ECI: Rz(raan) Rx(i) Rz(argp). Columns are P (perigee), Q, W (h-hat)."""
    return _rot_z(raan) @ _rot_x(i) @ _rot_z(argp)


def keplerian_to_cartesian(a: float, e: float, i: float, raan: float, argp: float,
                           nu: float, mu: float) -> Tuple[Vec3, Vec3]:
    """Return (r, v) in ECI [m, m/s]. Angles i, raan, argp, nu in radians, a in metres, mu in m^3/s^2.
    Valid for elliptical orbits (0 <= e < 1, a > 0)."""
    if a <= 0.0:
        raise ValueError("semi-major axis must be positive")
    if not 0.0 <= e < 1.0:
        raise ValueError("eccentricity must satisfy 0 <= e < 1")

    p = a * (1.0 - e**2)
    r_mag = p / (1.0 + e * np.cos(nu))

    # Perifocal frame
    r_pqw = r_mag * np.array([np.cos(nu), np.sin(nu), 0.0])
    v_pqw = np.sqrt(mu / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])

    # PQW -> ECI: Rz(raan) Rx(i) Rz(argp)
    Q = perifocal_to_eci(i, raan, argp)
    return Q @ r_pqw, Q @ v_pqw


def _signed_angle(a: Vec3, b: Vec3, axis: Vec3) -> float:
    """Angle from a to b in [-pi, pi], positive counter-clockwise about `axis` (right-handed)."""
    return float(np.arctan2(np.cross(a, b) @ axis, a @ b))


def cartesian_to_keplerian(r: Vec3, v: Vec3, mu: float, tol: float = 1e-10
                           ) -> Tuple[float, float, float, float, float, float]:
    """Return (a [m], e, i, raan, argp, nu) with angles in radians in [0, 2pi) (i in [0, pi]).

    Only bound (elliptical) orbits are supported; raises ValueError otherwise (including h = 0,
    i.e. r parallel to v or v = 0). Singular cases: circular -> argp = 0 and nu is measured from the
    ascending node; equatorial -> raan = 0 (node line taken along +x).
    """
    r_mag, v_mag = np.linalg.norm(r), np.linalg.norm(v)
    h = np.cross(r, v)
    h_mag = np.linalg.norm(h)
    if r_mag == 0.0 or h_mag <= tol * r_mag * max(v_mag, 1.0):
        raise ValueError("degenerate orbit: angular momentum is (near) zero (v = 0 or v parallel to r)")
    energy = 0.5 * v_mag**2 - mu / r_mag
    if energy >= 0.0:
        raise ValueError("state is not on a bound (elliptical) orbit: specific energy >= 0")

    a = -mu / (2.0 * energy)
    e_vec = ((v_mag**2 - mu / r_mag) * r - (r @ v) * v) / mu
    e = float(np.linalg.norm(e_vec))
    h_hat = h / h_mag
    inc = float(np.arccos(np.clip(h_hat[2], -1.0, 1.0)))

    node = np.array([-h[1], h[0], 0.0])                  # z x h, points to the ascending node
    if np.linalg.norm(node) > tol * h_mag:
        node_hat = node / np.linalg.norm(node)
        raan = float(np.arctan2(node_hat[1], node_hat[0]))
    else:
        node_hat, raan = np.array([1.0, 0.0, 0.0]), 0.0

    p_hat = e_vec / e if e > tol else node_hat            # perigee direction (node if circular)
    argp = _signed_angle(node_hat, p_hat, h_hat)
    nu = _signed_angle(p_hat, r / r_mag, h_hat)
    two_pi = 2.0 * np.pi
    return a, e, inc, raan % two_pi, argp % two_pi, nu % two_pi
