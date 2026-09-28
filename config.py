"""User-editable simulation inputs. Edit values here, then run `python main.py`."""

from dataclasses import dataclass, field
from typing import Tuple

Vec3 = Tuple[float, float, float]


@dataclass
class OrbitElements:
    """Classical Keplerian elements at the start of the epoch."""
    a_m: float = 7.0e6           # semi-major axis [m]
    e: float = 0.001             # eccentricity [-]
    i_deg: float = 51.6          # inclination [deg]
    raan_deg: float = 40.0       # right ascension of ascending node [deg]
    argp_deg: float = 30.0       # argument of perigee [deg]
    nu_deg: float = 0.0          # true anomaly at epoch [deg]


@dataclass
class Switches:
    """Independent on/off flags for each perturbation."""
    drag: bool = True
    j2: bool = True
    moon: bool = True
    sun: bool = True


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
    r_sun_m: Vec3 = (1.49598e11, 0.0, 0.0)                       # Earth->Sun [m] (placeholder)
    r_moon_m: Vec3 = (3.844e8, 0.0, 0.0)                         # Earth->Moon [m] (placeholder)


@dataclass
class SimConfig:
    orbit: OrbitElements = field(default_factory=OrbitElements)
    switches: Switches = field(default_factory=Switches)
    spacecraft: Spacecraft = field(default_factory=Spacecraft)
    environment: Environment = field(default_factory=Environment)


CONFIG = SimConfig()
