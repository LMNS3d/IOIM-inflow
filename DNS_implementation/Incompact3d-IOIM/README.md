# Incompact3d-IOIM

Derivative of [Xcompact3d](https://github.com/xcompact3d/Incompact3d) with an
inner-outer interaction model inlet for a zero-pressure-gradient turbulent
boundary layer. Modifications are by Anjia Ying and Lin Fu (2026). This is
not an official Xcompact3d release, and the Xcompact3d name is not used to
endorse it.

The upstream solver remains under the BSD 3-Clause License. Its copyright
notice is retained in `LICENSE` and in each modified source file. IOIM-S
additions are distributed under the same license.

Upstream documentation: https://xcompact3d.readthedocs.io/en/latest/

## What changed

- `iin = 4` is a recycling inlet. `method_rr = 1` is IOIM-S, `2` is LWS-S
  and `3` is LWS-R. Any other `iin` keeps the original Blasius inlet.
- A cold start (`irestart = 0`, `iin = 4`) uses the MATLAB mean field and a
  synthetic fluctuation field. `ioim_checkpoint.bin` is read only on restart.
- The example case is `TBL_case/input.i3d`. Mean profiles belong in
  `TBL_case/profiles/`, one `U_X*.txt` and `V_X*.txt` pair per streamwise
  station. They are produced by `mean_flow_generation`.
- An optional sponge can be placed at the outlet or the top boundary. It
  needs the same mean field; if the profiles are missing, the run warns and
  continues without the sponge. `iin = 4` without the profiles stops.

## Build

Install 2DECOMP&FFT first and point CMake at it. Do not rely on the automatic
download if that step fails.

```bash
export decomp2d_DIR=/path/to/2decomp/install/lib/cmake/decomp2d
cmake -S . -B build
cmake --build build
```

IOIM-S also needs FFTW and LAPACK.

## License

BSD 3-Clause. See `LICENSE`.
