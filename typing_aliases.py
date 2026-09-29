"""Type aliases so signatures make scalar vs. vector arguments explicit."""

import numpy as np
from numpy.typing import NDArray

Vec3 = NDArray[np.float64]   # shape (3,)  — 3-component vector (ECI, SI units unless noted)
Vec2 = NDArray[np.float64]   # shape (2,)  — 2-component in-plane vector
Mat3 = NDArray[np.float64]   # shape (3,3) — rotation matrix
