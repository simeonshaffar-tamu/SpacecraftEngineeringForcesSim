# Spacecraft Perturbation Simulation (Spacecraft Engineering HW2)

Computes instantaneous J2, drag, Sun, and Moon perturbing accelerations (m/s²) in an idealized ECI frame.
See [docs/SPEC.md](docs/SPEC.md) for inputs, outputs, models, assumptions, and the GenAI record.

## Setup & run
```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python main.py
```
Edit inputs in `config.py`.

## Layout
- `config.py` — inputs (orbit elements, switches, spacecraft, environment)
- `constants.py` — physical constants
- `orbital_elements.py` — Keplerian → Cartesian
- `perturbations.py` — acceleration models
- `main.py` — run and print
