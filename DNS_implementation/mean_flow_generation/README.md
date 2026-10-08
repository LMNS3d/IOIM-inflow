# Mean-flow generation for the IOIM inlet

MATLAB preprocessor for the Xcompact3d TBL case. It marches an
incompressible zero-pressure-gradient mean boundary layer on the DNS grid
and writes the inlet profiles read by `iin = 4`.

At each station the composite profile (van Driest inner layer plus a
Coles-type wake) is closed by a skin-friction correlation. Stations are
coupled by the von Kármán momentum integral. The wall-normal velocity is
recovered from continuity. Viscosity is chosen so that the inlet 99%
thickness equals `delta_in` in `TBL_case/input.i3d`.

This directory does not compare against the Torroja DNS. That check is in
`Validation of mean flow prediction`.

## Layout

```
main.m                 march on the DNS grid and write profiles/
src/f_velocity.m       profile closed on theta and the Cf law
src/f_Round1.m         one profile evaluation and the Cf residual
src/f_Scheme_GV.m      composite profile at a fixed wake factor
src/f_ITBL_mean.m      inner law and wake shape function
src/f_generate_grid.m  wall-normal grid used by the Fortran stretching
src/first_derivative.m non-uniform first derivative
profiles/              written on a successful run; not part of the source
```

`kf` and `C` printed by `main.m` are documentation. The values used in the
march are those inside `f_Round1`.

## Grid

`Nx`, `Ny`, `Lx`, `Ly`, `delta_in` and `dt` must match
`DNS implementation/Incompact3d-IOIM/TBL_case/input.i3d`. The wall-normal
nodes come from the same stretching as the Fortran solver.

## IOIM bounds

After the march, `main.m` prints bounds for `input.i3d`, evaluated at the
station nearest `40*delta_in`:

- `n_sp` no larger than the value that keeps `(n_sp*dt)+ <= 1`
- `N_I` at least the value that keeps the block length at 10 outer time units
- `N_Gamma` greater than `N_I`, and at least `N_I/2` plus that same minimum
- `N_averaging` at least `10*N_I`

`width_interp` is not estimated here.

## Run

From this directory, in MATLAB:

```matlab
main
```

Copy the generated `profiles/` directory to `TBL_case/profiles/`. Do not
replace `input.i3d` with that directory. `iin = 4` stops if `profiles/U_X1.txt`
is missing. A sponge without those profiles warns and is disabled.

## License

MIT. Copyright (c) 2026, Anjia Ying and Lin Fu. See `LICENSE`.
