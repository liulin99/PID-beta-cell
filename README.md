# Reproducibility package for “PID-controller enhanced artificial β-cells”

This archive contains all MATLAB scripts and embedded initial conditions needed to reproduce the figures in the paper. **No external datasets are required.**

## Contents
- `paper_code.m` — generates all figures except the cell-density figure.
- `cell_density.m` — generates the cell-density figure.

## Requirements
- MATLAB R2022b or later (uses `dictionary` and `tail`); tested on R2023b
- Toolboxes: base MATLAB only (`ode23`)
- OS tested: Windows 11 / macOS 14

## How to run
In MATLAB:
```matlab
% All figures except the cell-density figure
run('paper_code.m')

% Cell-density figure
run('cell_density.m')

## Runtime notes
These simulations are compute-intensive and may take a long time to run (potentially hours), depending on hardware and MATLAB version. The code is deterministic; long runtime is expected.
