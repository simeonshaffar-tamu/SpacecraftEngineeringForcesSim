"""Keplerian elements -> ECI Cartesian state (two-body, no propagation)."""

import numpy as np


def _rot_z(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _rot_x(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def keplerian_to_cartesian(a, e, i, raan, argp, nu, mu):
    """Return (r, v) in ECI [m, m/s]. Angles in radians, a in metres.
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
    Q = _rot_z(raan) @ _rot_x(i) @ _rot_z(argp)
    return Q @ r_pqw, Q @ v_pqw
