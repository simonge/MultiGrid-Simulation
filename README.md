# MG Geant4 Simulation

MG is a Geant4 simulation of the ESS Multigrid detector. It uses Geant4 for particle transport, ROOT for output, and NCrystalGeant4 for crystal-material support.

## Dependencies

Install or load a mutually compatible C++ toolchain and the following dependencies:

- CMake **3.30 or newer** and a C++17 compiler. Use the same compiler/toolchain that was used to build Geant4 and ROOT.
- Geant4 with GDML, UI, and visualization support, including the Geant4 data libraries.
- ROOT with the Core, RIO, Tree, and Hist components.
- Boost headers (used by the analysis code).
- Python 3 with `venv` support and pip. The setup script uses a repository-local virtual environment for the `ncrystal-geant4` Python package; it does not require Conda.
- NCrystal and the NCrystalGeant4 CMake package. `setup.sh` installs the Python bridge and passes its CMake package location to CMake; the required C++ packages Geant4 and ROOT must already be available in the active environment.

On the ESS/CERN LCG installation used for development, load the LCG view before configuring or running MG. For example:

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_110/x86_64-ubuntu2604-gcc15-opt/setup.sh
```

Use the LCG view appropriate to your machine and Geant4/ROOT installation. Avoid mixing libraries from unrelated system, Conda, or LCG installations.

## Geant4 hadronic physics data

The selected HP hadronic physics list requires proton cross-section data. Set `G4PARTICLEHPDATA` in the shell that will run MG:

```bash
export G4PARTICLEHPDATA=/cvmfs/geant4.cern.ch/share/data/G4TENDL1.4
```

This is the data path used on the development system. Ensure that it exists and contains the Geant4 TENDL data compatible with the Geant4 release in use. The environment variable must be set at runtime; setting it only while building is not sufficient. Other Geant4 dataset variables (such as `G4NEUTRONHPDATA`) are normally set by the Geant4/LCG environment setup.

## Configure and build

From the repository root, with Geant4 and ROOT available in the environment:

```bash
./setup.sh
```

The script creates `.venv/`, installs the pinned `ncrystal-geant4` package from `requirements.txt`, configures the CMake build in `build/`, and builds the `MG` executable. Generated build and environment directories are ignored by Git.

If `python3 -m venv` reports that the venv module is missing, install your distribution's Python venv package and retry. To configure/build manually, first install `ncrystal-geant4` in an environment where its `ncrystal-geant4` command is available, then pass the directory reported by `ncrystal-geant4 --cmakedir` as `NCrystalGeant4_DIR` to CMake.

## Run

Source the Geant4/ROOT environment and export the HP data variable in each new shell before running. Examples from the repository root:

```bash
./build/MG --help
./build/MG --preset EMMA
```

A preset named `NAME` selects `macros/NAME.mac` and, by default, `macros/NAME.dat`. The macro controls the run (including its event count) and can set the ROOT output filename with `/MG/analysis/setOutputFile`.

### Command-line options

```text
-m, --macro FILE       Execute a Geant4 macro (batch mode)
-o, --output FILE      Set the ROOT output path; overrides a macro's output setting
-f, --geometry FILE    Read detector geometry parameters from FILE
-r, --seed INTEGER     Set the random seed (0 chooses a time-based seed)
-n, --events INTEGER   Run this number of events without a macro
-u, --session TYPE     Select the Geant4 interactive UI session
    --preset NAME      Use macros/NAME.mac and macros/NAME.dat by default
-h, --help             Show help
```

For example, run the EMMA macro while choosing a new output file and seed:

```bash
./build/MG --preset EMMA --output output/emma.root --seed 12345
```

`--output` takes precedence over output-file commands in the macro. `--events` runs directly without a macro and therefore cannot be combined with `--macro` or a preset that selects a macro. Use the macro when you need its source, detector, or physics configuration.

Output paths are interpreted relative to the current working directory. Ensure the destination directory exists before starting a run.

## GitHub Actions

The workflow in `.github/workflows/build.yml` builds MG and runs a 100-event `EMMAtest0` smoke test. It uses the LCG view and CVMFS data paths above, so it requires a self-hosted GitHub Actions runner with CVMFS mounted and the `self-hosted`, `linux`, `x64`, and `lcg-110` labels. It is not configured for GitHub-hosted runners.
