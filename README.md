# E-bike Frame Improvement Project

This folder contains the modular MATLAB code for the CAMD e-bike frame FEM and topology optimization project.

## Required run order

1. `baseline/main.m`
2. `topology_optimization/main.m`
3. `redesign_and_validation/main.m`

Each module writes figures, CSV files and MAT files into output folders inside the same module folder.

## Module summary

- `baseline/`: reconstructs the 2D reference frame skeleton and topology design domain, generates the structured Q4 mesh, defines material/loads/supports and runs baseline FEM.
- `topology_optimization/`: runs SIMP minimum-compliance topology optimization from the baseline FEM data.
- `redesign_and_validation/`: interprets the topology result into the First Clean Redesign, validates it using FEM, and reports it as the final conceptual engineering design.

## MATLAB version and toolboxes

The scripts use base MATLAB functionality for plotting, sparse matrix assembly, linear solves and CSV/MAT export. No SolidWorks, ANSYS, PDE Toolbox or manual CAD import is required.

## Important modelling notes

- Units: mm, N, MPa.
- 2D plane-stress Q4 elements are used.
- The material is a simplified isotropic equivalent based on woven prepreg material values from Yasin et al. (2023), not a full composite laminate model.
- Void regions are retained with very small stiffness (`EminRatio`) to keep the global stiffness matrix solvable.
- The final result is the First Clean Redesign: a conceptual 2D MATLAB geometry with clear frame members and preserved battery cavity, not a certified 3D manufacturable e-bike frame.
