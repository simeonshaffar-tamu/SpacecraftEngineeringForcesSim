"""Physical constants (SI units). Earth values are the standard textbook set;
confirm against the values given in the assignment handout before final use."""

import numpy as np

# Earth
MU_EARTH = 3.986004418e14        # m^3/s^2
R_EARTH = 6.378137e6             # m (equatorial radius)
J2_EARTH = 1.08262668e-3         # dimensionless
OMEGA_EARTH = np.array([0.0, 0.0, 7.2921150e-5])  # rad/s, co-rotating atmosphere angular velocity (ECI)

# Third bodies (from assignment)
MU_SUN = 1.32712440018e20        # m^3/s^2
MU_MOON = 4.902800066e12         # m^3/s^2

# Idealized frame: ecliptic tilted about +x (vernal equinox direction) by the mean obliquity
OBLIQUITY_DEG = 23.4393
