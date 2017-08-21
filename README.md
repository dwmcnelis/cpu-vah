# cpu-vah

A CPU based (3+1) Dimensional Viscous Anisotropic Hydrodynamics Simulation.


## Prerequisites

| Target        | Purpose         |
| ------------- |-------------|
| **c++**       | C++ compiler (nvcc, gcc, or clang) |
| **libconfig**       | Configuration files |
| **gtest**        | C++ unit testing |



## Makefile targets

| Target        | Purpose         |
| ------------- |-------------|
| all      | Buld cpu-va executable |
| test      | Run test suite |
| hydro      | Run simulation |
| clean      | Remove intermediate build files |

## Usage

```
Usage: cpu-vh [OPTION...] ARG1 ARG2
Run -- A program to run a single viscous hydrodynamic simulation of a
relativistic heavy ion collision

  -c, --config=CONFIG_DIRECTORY   Path to configuration directory
  -h, --hydro[=RUN_HYDRO]    Run hydrodynamic simulation
  -o, --output=OUTPUT_DIRECTORY   Path to output directory
  -t, --test[=RUN_TEST]      Run software tests
  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version
```