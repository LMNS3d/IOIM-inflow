# Validation of the mean-flow prediction

MATLAB comparison of the composite boundary-layer march with the
UPM/Torroja zero-pressure-gradient DNS of Sillero, Jiménez and Moser. The
march itself is the same model used by `mean_flow_generation`: a van Driest
inner layer plus a Coles-type wake, closed by a skin-friction correlation
and advanced with the von Kármán momentum integral. The wall-normal velocity
is recovered from continuity.

This directory does not write the Xcompact3d inlet profiles.

## Layout

```
main.m                     load DNS, march, write figures
src/f_velocity.m           profile closed on theta and the Cf law
src/f_Round1.m             one profile evaluation and the Cf residual
src/f_Scheme_GV.m          composite profile at a fixed wake factor
src/f_ITBL_mean.m          inner law and wake shape function
src/first_derivative.m     non-uniform first derivative
src/tight_subplot.m        axes layout (local reimplementation)
data/                      user-supplied MAT files (not distributed)
figures/                   written on a successful run
```

`f_velocity` returns the 99% thickness in the sixth position. That output
is a thickness, not a momentum thickness.

## Method

Inner contribution, integrated in wall units with van Driest damping:

```
d+      = (1 - exp(-y+/A+))^2,   A+ = 17
nu_t+   = kappa * y+ * d+,       kappa = 0.41
dU+/dy+ = 1 / (1 + nu_t+)
```

Wake contribution, integrated in outer units and scaled by a factor Pi:

```
dU_wake+/dy = (pi/kappa) / Re_tau * sin(pi * y/delta)
```

Pi is the root of `Cf_correlation/Cf_profile - 1 = 0`, with

```
Cf = 2 / (log(Re_theta)/kf + C)^2,   kf = 0.344, C = 1.770
```

The march is forward Euler on `d(theta)/dx = Cf/2`, with Cf evaluated at
the midpoint Reynolds number. Viscosity is fixed by matching `delta_99` at
the inlet. The inlet is the first DNS station with `x >= 10 * delta_99(1)`.

## Data

Download the database from

https://torroja.dmt.upm.es/turbdata/blayers/high_re/

and place two MAT files in `data/`. Both must store the record in a variable
named `val`. These files are derived from the DNS and are not included.
Cite the DNS if the figures are published; the database terms are those of
the Torroja Laboratory, not this MIT license.

`tbl_xdata_structure.mat` — streamwise records, fields

| field   | meaning                                      |
|---------|----------------------------------------------|
| xgrid   | streamwise grid                              |
| cf      | skin friction                                |
| Reth    | momentum-thickness Reynolds number           |
| Retau   | friction Reynolds number                     |
| d99     | 99% thickness                                |
| dels    | displacement thickness                       |
| utau    | friction velocity                            |
| theta   | momentum thickness                           |

`structure_all_profiles.mat` — struct array, one entry per profile, fields

| field | meaning                                              |
|-------|------------------------------------------------------|
| d99, Retau, Reth, dels, utau | station identifiers used to match xgrid |
| yd99  | wall-normal grid normalized by the local d99         |
| ua    | mean streamwise velocity in wall units               |
| va    | mean wall-normal velocity in wall units              |

## Run

From this directory, in MATLAB R2020a or later (`exportgraphics`, `subtitle`):

```matlab
main
```

Figures are written to `figures/`:

- `PredictionM_Cf_Re.pdf` — Cf and Re_theta against the DNS
- `PredictionM_MeanUV.pdf` — mean profiles, staggered by station
- `test_vk.pdf` — momentum-integral residual of the DNS reference
- `stations_cf.pdf` — stations selected for the profile comparison

## Dependency note

Figure axes use a local `tight_subplot` with the same call signature as
Pekka Kumpulainen's File Exchange utility (id 27991, BSD). That utility is
not redistributed. `src/tight_subplot.m` is an independent reimplementation
used only for layout. Replacing the two call sites with `tiledlayout` is
straightforward on R2019b or later if preferred.

## Citation

If you use the reference data, cite

J. A. Sillero, J. Jiménez and R. D. Moser, "One-point statistics for
turbulent wall-bounded flows at Reynolds numbers up to δ+ ≈ 2000",
Phys. Fluids 25, 105102 (2013).

## License

MIT. Copyright (c) 2026, Anjia Ying and Lin Fu. See `LICENSE`.
