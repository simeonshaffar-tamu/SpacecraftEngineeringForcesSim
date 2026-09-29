"""User-editable simulation inputs. Edit values here, then run `python main.py`."""

from dataclasses import dataclass, field
from typing import Literal, Tuple

Vec3Tuple = Tuple[float, float, float]


@dataclass
class OrbitElements:
    """Classical Keplerian elements at the start of the epoch."""
    a_m: float = 7.0e6           # semi-major axis [m]
    e: float = 0.001             # eccentricity [-]
    i_deg: float = 51.6          # inclination [deg]
    raan_deg: float = 40.0       # right ascension of ascending node [deg]
    argp_deg: float = 30.0       # argument of perigee [deg]
    nu_deg: float = 30.0         # true anomaly at epoch [deg]


@dataclass
class CartesianState:
    """Initial state in the idealized ECI frame (used when SimConfig.input_mode == "cartesian")."""
    r_m: Vec3Tuple = (-8e6, 0.0, 0.0)           # position [m]
    v_mps: Vec3Tuple = (0.0, 7546.05, 0.0)       # velocity [m/s] (default ~ circular speed at 7000 km)


@dataclass
class Switches:
    """Independent on/off flags for each perturbation."""
    drag: bool = False
    j2: bool = False
    moon: bool = True
    sun: bool = False


@dataclass
class Spacecraft:
    mass_kg: float = 250.0
    drag_area_m2: float = 2.0
    cd: float = 2.2


@dataclass
class Environment:
    """Environment state at the epoch, in the same ECI frame as the spacecraft.
    Density and body positions are inputs (evaluated at the same instant)."""
    rho_kg_m3: float = 1.0e-12                                   # atmospheric density at s/c
    r_sun_m: Vec3Tuple = (1.49598e11, 0.0, 0.0)                  # Earth->Sun [m] (placeholder)
    r_moon_m: Vec3Tuple = (-2.71811e8, 0.0, 0.0)                 # Earth->Moon [m] (placeholder, in ecliptic plane, 180 deg from Sun)


@dataclass
class SimConfig:
    # "keplerian": start from `orbit` (6 elements, converted internally to r, v).
    # "cartesian": start from `cartesian` (r, v); orbit elements are derived for the plots only.
    input_mode: Literal["keplerian", "cartesian"] = "cartesian"
    orbit: OrbitElements = field(default_factory=OrbitElements)
    cartesian: CartesianState = field(default_factory=CartesianState)
    switches: Switches = field(default_factory=Switches)
    spacecraft: Spacecraft = field(default_factory=Spacecraft)
    environment: Environment = field(default_factory=Environment)


CONFIG = SimConfig()
