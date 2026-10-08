Incompact3d-IOIM installation
=============================

This solver is a derivative of Xcompact3d. The build follows the upstream
procedure, with two extra libraries required by the IOIM-S inlet.

## Source

Clone this repository, not the upstream Xcompact3d tree:

```bash
git clone https://github.com/xcompact3d/Incompact3d
```

Replace the address above with the URL of this derivative. Upstream
documentation remains at https://xcompact3d.readthedocs.io/en/latest/.

If you are behind a firewall, force Git to use HTTPS:

```bash
git config --global url."https://".insteadOf git@
```

Set `https_proxy` and `http_proxy` if a proxy is required.

## Minimum requirements

- CMake 3.20 or newer
- a recent Fortran compiler, such as gfortran 9 or newer
- an MPI library with a Fortran wrapper
- 2DECOMP&FFT 2.0.3, the version used by the upstream tests
- FFTW 3 and LAPACK, required by the IOIM-S inlet

On Debian or Ubuntu the two inlet libraries are:

```bash
sudo apt install libfftw3-dev liblapack-dev pkg-config
```

CMake finds FFTW through pkg-config. If FFTW is installed elsewhere:

```bash
export PKG_CONFIG_PATH=/path/to/fftw/lib/pkgconfig:${PKG_CONFIG_PATH}
```

This FFTW dependency is separate from the FFT backend inside 2DECOMP&FFT.
The Generic backend is enough for 2DECOMP. The inlet still needs `libfftw3`.

## Compilation

Point `FC` at the MPI Fortran wrapper so the compiler and MPI agree:

```bash
export FC=my_mpif90
```

Configure and build from the repository root:

```bash
cmake -S . -B build
cmake --build build -j
```

`cmake --build build -j -v` prints the compiler flags. Install with:

```bash
cmake --install build
```

The default prefix is `build/opt`. It contains `bin/xcompact3d` and example
input files. Change it at configure time:

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/path/to/install
```

The configuration can also be edited with `ccmake build`.

## Build with an existing 2DECOMP&FFT

If `decomp2d` is not already visible to CMake, the build tries to download
2DECOMP&FFT 2.0.3 and build it with the Generic FFT backend. To use an
existing installation instead, set

```bash
export decomp2d_DIR=/path/to/2decomp/install/opt/lib/cmake/decomp2d
```

On some systems `lib` is `lib64`. `decomp2d_DIR` must contain
`decomp2dConfig.cmake`. Installation of 2DECOMP&FFT is described at
https://github.com/2decomp-fft/2decomp-fft/blob/dev/INSTALL.md.

Other 2DECOMP backends are optional and are not required for IOIM-S:

- FFTW or MKL as the 2DECOMP FFT engine
- ADIOS2 instead of MPI-IO
- a single-precision build

For an MKL FFT backend:

```bash
export MKL_DIR=${MKLROOT}/lib/cmake/mkl
```

For ADIOS2:

```bash
cmake -S . -B build -DIO_BACKEND=adios2 -Dadios2_DIR=/path/to/adios2/install/lib/cmake/adios2
```

## Running the boundary-layer case

The IOIM example is `TBL_case/input.i3d`. Generate the mean field with
`mean_flow_generation` and copy its `profiles/` directory to
`TBL_case/profiles/` before a run with `iin = 4`.

## Examples and tests

Other input files are under `examples/`. The upstream test suite has 13
tests; see `tests/README.md`. Those tests do not exercise the IOIM inlet.

## Known issues

CTest looks for `mpirun`. If CMake picks the wrong executable:

```bash
cmake -S . -B build -DMPIEXEC_EXECUTABLE=/correct/path/to/mpirun
```

If configuration fails while downloading 2DECOMP&FFT, remove
`build/decomp2d-build` and point `decomp2d_DIR` at a preinstalled library
rather than retrying the download.
