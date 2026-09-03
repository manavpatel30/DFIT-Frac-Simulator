# Sleeve-Based DFIT Simulator

A compact, physics-based Python model of a sleeve-operated Diagnostic Fracture Injection Test (DFIT). The project follows the pressure response from sealed-wellbore pressurization through fracture propagation and post-shut-in decline.

The repository provides two ways to explore the model:

- `dfit_v4.py` runs the example case and produces Matplotlib figures.
- `dfit_app.py` starts a local interactive dashboard for changing parameters and comparing results.

> This is a simplified educational and research prototype. It is intended to make the governing physics transparent; it is not a substitute for a calibrated commercial simulator or field-design workflow.

## Model workflow

| Stage | Physical process | Main representation |
| --- | --- | --- |
| 0 | Sleeve closed and wellbore pressurizing | Compressible wellbore storage with a ramped rate |
| 1 | Sleeve open, fracture initiation, and pumping | Fixed-height PKN fracture, Carter leakoff, and empirical friction loss |
| 2 | Injection stopped and pressure declining | Fixed-length fracture, pressure-dependent leakoff, elastic aperture, and 1D slot flow |

For equations and assumptions, see [Model overview](docs/MODEL_OVERVIEW.md). For all dashboard inputs, see the [Parameter guide](docs/PARAMETER_GUIDE.md).

## Repository structure

```text
dfit-fracture-simulator/
├── dfit_v4.py                 # Simulation engine and plotting example
├── dfit_app.py                # Local browser dashboard
├── requirements.txt           # Python dependencies
├── README.md                  # Project introduction and quick start
├── SOURCE_INTEGRITY.md        # Verification that source code is unchanged
└── docs/
    ├── GITHUB_UPLOAD.md       # Publishing instructions
    ├── MODEL_OVERVIEW.md      # Physics, equations, and limitations
    └── PARAMETER_GUIDE.md     # Inputs, units, defaults, and effects
```

## Requirements

- Python 3.9 or newer is recommended.
- NumPy
- SciPy
- Matplotlib

The dashboard uses only Python's standard-library web server on the backend and browser-native HTML, CSS, and JavaScript on the frontend. Matplotlib is still required because the dashboard loads the simulation engine from `dfit_v4.py`.

## Installation

Clone or download the repository, then work from its root directory:

```bash
python -m venv .venv
```

Activate the environment:

```bash
# macOS or Linux
source .venv/bin/activate

# Windows PowerShell
.venv\Scripts\Activate.ps1
```

Install the dependencies:

```bash
python -m pip install -r requirements.txt
```

## Run the interactive dashboard

From the repository root, run:

```bash
python dfit_app.py
```

Open [http://localhost:5050](http://localhost:5050) in a browser. Change values in the parameter panel and select **Re-run Simulation**. Stop the server with `Ctrl+C`.

The two Python files must remain together, and the command must be run from the repository root, because `dfit_app.py` loads `dfit_v4.py` by its relative filename.

## Run the plotting example

```bash
python dfit_v4.py
```

The script runs the default three-stage case, prints summary values, and opens the generated Matplotlib figures. The example inputs are defined in the `if __name__ == "__main__"` block near the end of the file.

## Outputs

The model reports and visualizes quantities including:

- wellbore and fracture pressure;
- breakdown pressure and ISIP;
- injection rate and friction-pressure loss;
- fracture half-length and average width;
- injected, leaked-off, and fracture-fluid volumes;
- shut-in pressure, aperture, and leakoff-rate histories; and
- animated fracture geometry in the dashboard.

## Important scope notes

- The fracture is represented as a symmetric, planar, fixed-height PKN fracture.
- Leakoff is based on a Carter-style relationship rather than a full reservoir-flow solution.
- The shut-in calculation holds fracture length fixed.
- Closure uses a linear-elastic aperture relation with a prescribed residual aperture.
- The `initiated` value is a diagnostic flag. In the current source, Stage 1 arrays are still evaluated when the breakdown criterion is not met, so those cases should be interpreted cautiously.
- Thermal effects, poroelastic stress changes, proppant, multiphase flow, natural fractures, layering, and complex fracture networks are outside the model.

See [Model overview](docs/MODEL_OVERVIEW.md) for the complete set of documented assumptions.

## Reference used in the source

Valkó and Economides, *Hydraulic Fracture Mechanics* (1995).

## Publishing on GitHub

The folder is ready to become a repository. See [GitHub upload instructions](docs/GITHUB_UPLOAD.md) for browser and command-line options.

## License

No license has been added because a license is a project-owner decision. Before publishing the repository for reuse, add the license that matches the permissions you want to grant. Without a license, normal copyright restrictions apply even if the repository is public.

