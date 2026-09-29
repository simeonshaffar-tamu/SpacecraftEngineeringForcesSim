# Spacecraft Perturbation Simulation (Spacecraft Engineering HW2)

Computes instantaneous J2, drag, Sun, and Moon perturbing accelerations (m/s²) in an idealized ECI frame.
See [docs/SPEC.md](docs/SPEC.md) for inputs, outputs, models, assumptions, and the GenAI record.

## Setup & run
```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python main.py
```
Edit inputs in `config.py`. `input_mode` selects Keplerian elements (`"keplerian"`) or a direct Cartesian ECI state (`"cartesian"`, e.g. r = [7e6, 0, 0] m).

## Layout
- `config.py` — inputs (orbit elements, switches, spacecraft, environment)
- `constants.py` — physical constants
- `typing_aliases.py` — `Vec3` (shape-3 array), `Vec2`, `Mat3` aliases used in signatures
- `orbital_elements.py` — Keplerian ↔ Cartesian
- `perturbations.py` — acceleration models
- `plotting.py` — orbit-plane and Sun/Moon views; interactive combined figure with true-anomaly slider
- `simulation.py` — config + true anomaly -> state and accelerations
- `main.py` — run and print (`--no-show` skips the interactive window)
