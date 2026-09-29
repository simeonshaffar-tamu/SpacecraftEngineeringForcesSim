# Part (a) — Specification and GenAI Record

## Purpose
Compute instantaneous perturbing accelerations on a spacecraft (Earth J2, atmospheric drag,
Sun and Moon differential point-mass gravity) at a single instant. No orbit propagation.

## Inputs (all defined in `config.py`)
| Group | Quantity | Units |
|---|---|---|
| Orbit (Keplerian, epoch) | a, e, i, RAAN, argument of perigee, true anomaly | m, –, deg, deg, deg, deg |
| Switches | `drag`, `j2`, `moon`, `sun` (independent booleans) | – |
| Spacecraft | mass m, drag area A_D, drag coefficient C_D | kg, m², – |
| Environment | density ρ at spacecraft; Earth→Sun and Earth→Moon position vectors (same ECI frame, same instant as ρ) | kg/m³, m |

**Alternative input:** set `input_mode = "cartesian"` in `config.py` and give ECI position `r_m` [m] and velocity
`v_mps` [m/s] directly (as in the assignment's test cases). Accelerations are computed from that state as-is. For the
plots, orbit elements are derived from (r, v); if the state is not on a bound, non-degenerate orbit (e.g. v = 0 or
v parallel to r), the accelerations are still printed and the plots are skipped.

Keplerian elements are converted internally to Cartesian ECI position **r** [m] and velocity **v** [m/s]
(`orbital_elements.py`); the acceleration function takes Cartesian **r**, **v**.

## Outputs (printed to terminal, m/s²)
Drag, J2, Moon, and Sun acceleration vectors (ECI, 3 components each), plus their vector sum.
A disabled contribution is a zero vector and is excluded from the sum.

## Models
- **J2:** a = −(3/2) J2 μ R_E² / r⁵ · [x(1 − 5z²/r²), y(1 − 5z²/r²), z(3 − 5z²/r²)]
- **Drag:** v_rel = v − ω_E × r, a = −½ ρ C_D (A_D/m) |v_rel| v_rel, ω_E = [0, 0, 7.2921150e-5] rad/s
- **Third body (Sun/Moon):** a = μ_b [ (r_b − r)/|r_b − r|³ − r_b/|r_b|³ ]
- Earth's central point-mass term is **excluded** (a propagator adds it separately).

## Constants
μ_Sun = 1.32712440018e20 m³/s², μ_Moon = 4.902800066e12 m³/s² (from assignment).
Earth μ, R_E, J2 (`constants.py`) are standard values — **verify against the assignment handout**.

## Assumptions
- Idealized ECI frame, z along Earth's spin axis; precession/nutation neglected (not an exact J2000 model).
- Atmosphere rigidly co-rotates with Earth; density is an input (no atmosphere model).
- Constant C_D, A_D (no attitude dependence); J3 and higher harmonics ignored.
- Sun/Moon are point masses at the supplied positions; supplied positions are placeholders until real ephemerides are used.
- Elliptical orbits only (0 ≤ e < 1). SI units internally.
- Invalid inputs (mass ≤ 0, negative area/C_D/density, zero position, non-finite vectors) raise `ValueError`.

## GenAI record
- **Tool/model:** Claude Code (Claude Sonnet 5.5, `claude-sonnet-5-5`)
- **Access date:** 2026-09-28
- **Initial prompt (summary):** Create a new repo for a spacecraft simulation for a GenAI-assisted assignment
  (Problem 3, perturbing accelerations). Inputs: config struct in a Python file with 6 Keplerian elements
  (converted internally to Cartesian) and on/off switches per acceleration. Outputs: drag, J2, Moon and Sun
  acceleration vectors (m/s²). Do not propagate; compute the state at epoch from true anomaly. Only step (a) for now.
  (Full assignment text pasted alongside.)
- **Key follow-up:** Use Python with a virtual environment; inputs live in a `.py` file, outputs are terminal printouts.
- **Note:** Spacecraft parameters, density and Sun/Moon vectors were added to the config by the assistant because the
  function needs them (per the assignment); Earth constants are assistant-supplied standard values.

## Diagrams (`plotting.py`, saved to `outputs/`)
1. **Orbit plane** — viewed along the orbit's angular momentum vector.
2. **Sun/Moon/orbit** — viewed along Earth's orbital angular momentum about the Sun (ecliptic north).
   The ecliptic normal is fixed in the idealized ECI frame as [0, −sin ε, cos ε], ε = 23.4393°, with +x ECI as the
   vernal equinox direction. Everything is projected onto the ecliptic plane. Earth (translucent), the orbit, the
   nodes and the spacecraft are to scale relative to each other; the Sun and Moon are far too distant to fit, so they
   are markers along their true projected directions (true distances in the legend). Filled circle = ascending node
   (RAAN), open circle = descending node (LDDN). Both diagrams show the spacecraft and direction-only acceleration
   arrows, and a slider sets the true anomaly in the interactive figure.
   The Sun/Moon plot also shows Earth's equator (tilted ellipse) and spin axis (N toward viewer, S behind). Depth cues:
   near-side orbit/equator are solid and drawn over the translucent Earth; far-side parts are dashed and drawn beneath it;
   the spacecraft is hollow with faded arrows when behind.
