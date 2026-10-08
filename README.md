# Open-source code for the IOIM inflow

Code accompanying the inner-outer interaction model (IOIM-S) inlet for a
zero-pressure-gradient turbulent boundary layer. The three parts have
different licenses. Use the `LICENSE` in each directory; this file only
says where things are.

```
DNS_implementation/Incompact3d-IOIM/          Fortran solver, derivative of Xcompact3d
DNS_implementation/mean_flow_generation/      MATLAB mean field written for that solver
Validation_of_mean_flow_prediction/           MATLAB comparison with the Torroja DNS
```

## What to run

1. March the mean field in `mean_flow_generation` and copy `profiles/` to
   `Incompact3d-IOIM/TBL_case/profiles/`. `Nx`, `Ny`, `Lx`, `Ly`, `delta_in`
   and `dt` must match `TBL_case/input.i3d`.
2. Build and run `Incompact3d-IOIM`. `iin = 4` is the recycling inlet
   (`method_rr = 1` is IOIM-S). Any other `iin` keeps the original Blasius inlet.
3. `Validation of mean flow prediction` only checks the march against the
   Sillero, Jiménez and Moser DNS. It does not write the solver profiles.
   Download the database separately; it is not included here.

## Copyright

- `Incompact3d-IOIM` is BSD 3-Clause. The Xcompact3d copyright (2012-2022)
  is retained. Modifications are Copyright (c) 2026, Anjia Ying and Lin Fu,
  under the same license. This is not an official Xcompact3d release.
- Both MATLAB directories are MIT, Copyright (c) 2026, Anjia Ying and Lin Fu.
- The Torroja DNS files are not covered by either license. Cite
  Sillero, Jiménez and Moser, Phys. Fluids 25, 105102 (2013), if they are used.
