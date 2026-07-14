# Binding ΔΔG Analysis of KCNA5 T1 Domain Mutations

## Overview

This directory calculates the effect of all possible missense mutations at the
KCNA5 T1 domain inter-subunit interface on binding energy (ΔΔG), using FoldX.
Mutations are run across both orientations of the inter-subunit interface (chain A→D
and chain D→A), and results are cross-referenced against COSMIC COAD mutations to
identify cancer-associated mutations with significant effects on interface stability.

## Directory structure

```
binding_ddg/
├── data/
│   └── SWISSMODEL_WT_for_ddg.pdb      # Prepared WT structure used as FoldX input
│
├── results/
│   ├── Binding_DDG_AD.csv             # ΔΔG results: chain A mutated, chain D as partner
│   └── Binding_DDG_DA.csv             # ΔΔG results: chain D mutated, chain A as partner
│
├── binding_ddg.py                     # Runs FoldX BuildModel + AnalyseComplex for all
│                                      # interface residue mutations
├── summary_binding_ddg.py             # Collects FoldX .fxout files into a summary CSV
└── ddg_analysis.ipynb                 # Analyses and visualises ΔΔG results
```

## Input structure preparation

The wildtype structure (`data/SWISSMODEL_WT_for_ddg.pdb`) was obtained from the
SWISS-MODEL repository:

> https://swissmodel.expasy.org/repository/uniprot/P22460?range=120-527&template=9eef.1.D

This model covers the KCNA5 T1 domain (residues 120–527). The structure contains
multiple chains representing the tetrameric T1 assembly. Prior to running `binding_ddg.py`,
the chain identifiers were manually reassigned:

- Chain A (T1 subunit) → **Chain E** (protein)
- Chain D (T1 subunit) → **Chain F** (ligand/partner)

This renaming is required because `binding_ddg.py` treats one chain as the protein
and the other as the binding partner (ligand), and FoldX requires distinct chain labels.

## Scripts

### `binding_ddg.py`

Adapted from: https://github.com/shorthouse-lab/binding_ddg

For each residue on chain E (protein) within 10 Å of chain F (ligand/partner),
this script:
- Generates all 19 possible amino acid substitutions using FoldX `BuildModel`
- Runs FoldX `AnalyseComplex` on both wildtype and mutant structures to obtain
  interaction energies
- Outputs per-residue folders (`residue_XXX/`) containing FoldX `.fxout` files

Two runs were performed:
- **AD**: chain E = protein, chain F = partner → `Binding_DDG_AD.csv`
- **DA**: chain F = protein, chain E = partner → `Binding_DDG_DA.csv`

This captures mutations on both sides of the inter-subunit interface.

FoldX settings used:
- pH 7.0, ionic strength 0.05 M, `water=CRYSTAL`, `vdwDesign=2`
- FoldX version: 5.1

To run:
```bash
python binding_ddg.py data/SWISSMODEL_WT_for_ddg.pdb
```

### `summary_binding_ddg.py`

Run in the directory containing the `residue_*/` folders output by `binding_ddg.py`:

```bash
python summary_binding_ddg.py
```

Reads `Summary_*_AC.fxout` files from each residue folder, extracts interaction
energies for WT and mutant, computes ΔΔG = MUT − WT, and writes
`Binding_DDG_summary.csv`.

Output columns:
- `Mutation`: FoldX mutation code (e.g. `FE196D` = Phe196→Asp on chain E)
- `WT_energy`: wildtype interaction energy (kcal/mol)
- `MUT_energy`: mutant interaction energy (kcal/mol)
- `DDG`: ΔΔG (kcal/mol); positive = destabilising

### `ddg_analysis.ipynb`

This notebook:
- Loads `Binding_DDG_AD.csv` and `Binding_DDG_DA.csv`, assigns chain labels (A and D
  respectively), and concatenates into a single dataframe
- Reformats mutation codes from FoldX notation (e.g. `FE196D`) to standard notation
  (e.g. `F196D`)
- Filters for mutations with |ΔΔG| ≥ 0.8 kcal/mol as functionally significant
- Cross-references significant mutations against COSMIC COAD mutations
  (`KCNA5_T1_mutations.csv`)
- Plots the full ΔΔG distribution (all interface mutations) overlaid with the
  distribution for COSMIC-observed mutations, with 90th and 95th percentile thresholds

## Dependencies

### Python
- biopython
- pandas, numpy, matplotlib, seaborn

### External
- FoldX 5.1 (https://foldxsuite.crg.eu/) — must be installed and path set in
  `binding_ddg.py` (`foldxdir` variable, line 14)