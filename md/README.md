# Molecular Dynamics Analysis of KCNA5 T1 Domain (WT and D166N)

## Overview

This directory contains GROMACS input files and Python analysis scripts for 100 ns
molecular dynamics simulations of the KCNA5 T1 domain tetramer, comparing wildtype
(WT) and the D166N cancer-associated mutant. Analysis focuses on inter-subunit
distances at the buried polar interface formed by D166 and its hydrogen-bonding
partners (T133, Q134, T137, Q180), across three independent repeats per system.


## Directory structure

```
md/
├── KCNA5_WT_100ns/
│   └── gromacs/
│       ├── toppar/                     # Force field parameter files
│       ├── index.ndx                   # GROMACS index file
│       ├── step6.6_equilibration.gro   # Starting structure after equilibration
│       ├── step7_production.mdp        # Production run parameters
│       └── topol.top                   # Topology file
│
├── KCNA5_D166N_100ns/
│   └── gromacs/                        # Same structure as WT above
│       ├── toppar/                     
│       ├── index.ndx                   
│       ├── step6.6_equilibration.gro   
│       ├── step7_production.mdp        
│       └── topol.top                   
├── simple_mda/                         # MDAnalysis utility package (see below)
│   └── src/
│       ├── __init__.py
│       ├── load_md_npz.py
│       ├── radial_orientation.py
│       ├── simple_angles.py
│       ├── simple_calc.py
│       ├── simple_dist.py
│       ├── simple_helix.py
│       ├── simple_ks.py
│       ├── simple_pca.py
│       ├── simple_rmsf.py
│       └── sliding_window.py
│   └── pyproject.toml
│
├── dist.ipynb                          # Computes inter-residue distances from trajectories
└── dist_analysis.ipynb                 # Analyses and visualises distance distributions
```

## Simulation details

| Parameter        | Value                          |
|------------------|-------------------------------|
| Software         | GROMACS                        |
| System           | KCNA5 T1 domain tetramer       |
| Simulation time  | 100 ns production              |
| Repeats          | 3 per system (r1, r2, r3)      |
| Systems          | WT, D166N                      |

Equilibration was performed prior to production; the starting structure for production
is `step6.6_equilibration.gro`. Production run parameters are in `step7_production.mdp`.

## Running the simulation

The equilibration steps (steps 1–6.6) are pre-run; the provided `.gro` file is the
starting point for production. To run the production simulation from the supplied files:

```bash
# 1. Generate the run input file (.tpr) from the equilibrated structure
gmx grompp \
  -f step7_production.mdp \
  -c step6.6_equilibration.gro \
  -p topol.top \
  -n index.ndx \
  -o production.tpr

# 2. Run production MD
gmx mdrun \
  -v \
  -deffnm production \
  -ntmpi 1 \
  -ntomp 8    # adjust to your available CPU cores
```


Repeat for each replicate (r1, r2, r3) and each system (WT, D166N). Rename output
`.tpr` and `.xtc` files to match the naming convention expected by the analysis
notebooks (e.g. `WT_r1.tpr`, `WT_r1.xtc`).

## Data

Trajectory files (`.xtc`) and compiled topology files (`.tpr`) are not included
due to file size. Place them in a `data/` directory before running the notebooks:

```
md/data/
├── WT_r1.tpr / WT_r1.xtc
├── WT_r2.tpr / WT_r2.xtc
├── WT_r3.tpr / WT_r3.xtc
├── D166N_r1.tpr / D166N_r1.xtc
├── D166N_r2.tpr / D166N_r2.xtc
└── D166N_r3.tpr / D166N_r3.xtc
```

Large trajectory files may be deposited on Zenodo or similar if needed for
full reproducibility.

## `simple_mda` package

`simple_mda` is a local MDAnalysis utility package providing helper functions
for distance, angle, RMSF, PCA, and helix analysis. Install it before running
the notebooks:

```bash
pip install -e simple_mda/
```

## Reproduction

### Step 1 — Compute inter-residue distances (`dist.ipynb`)

Calculates distances between D166 (CG atom) and four hydrogen-bonding partner
residues across chains A and D (subunit pair AD) for all repeats:

| Distance measured       | Atoms        |
|-------------------------|--------------|
| D166 – T133             | CG – OG1     |
| D166 – Q134             | CG – NE2     |
| D166 – T137             | CG – OG1     |
| D166 – Q180             | CG – NE2     |

The tetramer has a chain length of 408 residues; residue indices are adjusted
accordingly for each subunit pair. Results are saved as pickle files in `results/`:

```
results/
├── WT_T133_dist.pkl
├── WT_Q134_dist.pkl
├── WT_T137_dist.pkl
├── WT_Q180_dist.pkl
├── D166N_T133_dist.pkl
├── D166N_Q134_dist.pkl
├── D166N_T137_dist.pkl
└── D166N_Q180_dist.pkl
```

### Step 2 — Analyse and visualise distances (`dist_analysis.ipynb`)

Loads the pickle files and produces:

- **Correlation matrix**: Pearson r between the four distance measurements
  within WT and D166N, to assess whether perturbation of D166 propagates
  through the interface network


## Dependencies

### External
- GROMACS

### Python
- MDAnalysis
- simple_mda (local, see above)
- numpy, pandas, matplotlib, seaborn, scipy
- nglview (trajectory visualisation in Jupyter)
- pickle (standard library)