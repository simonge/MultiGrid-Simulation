# MG Geant4 Simulation

MG is a Geant4 simulation of the European Spallation Source (ESS) Multigrid neutron detector. Its purpose is to model particle transport through the detector and surrounding geometry, and to study the detector response to neutron beams and other particles, including scattering and energy deposition. Configurable geometries and source/run settings are supplied through parameter files and Geant4 macros. NCrystalGeant4 provides crystal-material support for simulations that use crystalline materials.

## What the simulation tracks

MG transports particles using Geant4 and records accepted steps in its sensitive detector volumes. These include Multigrid Ar/CO$_2$ detector gas, optional $^3$He detector volumes, sample/scattering volumes, the beam monitor, and selected test or surrounding-gas regions. The sensitive-detector code applies volume-specific filters: for example, gas hits with zero energy deposition are generally discarded, while sample hits are selected using neutron and volume-boundary conditions. A recorded row is therefore a detector hit/step, not a complete record of every Geant4 step in the simulation.

For each event with at least one accepted hit, MG writes one entry to a ROOT `TTree` named `PD`. The entry has event-level primary-beam fields plus variable-length hit arrays:

| Branch | Contents |
| --- | --- |
| `PD_Nhits` | Number of accepted hits recorded for the event; also the length of the hit arrays. |
| `PD_T0n`, `PD_E0n` | Source neutron time-of-flight and energy values supplied by the primary generator. Their interpretation depends on the configured source mode. |
| `PD_pdg` | Hit particle PDG code. |
| `PD_vid` | Detector voxel, cell, sample, or monitor identifier, depending on the sensitive volume. |
| `PD_m` | Particle rest mass. |
| `PD_t` | Global hit time. |
| `PD_x`, `PD_y`, `PD_z` | Hit position. |
| `PD_vx`, `PD_vy`, `PD_vz` | Primary vertex position associated with the event. |
| `PD_px`, `PD_py`, `PD_pz` | Particle momentum components at the hit. |
| `PD_th`, `PD_ph`, `PD_p` | Momentum direction angles and magnitude. |
| `PD_E` | Particle kinetic energy at the hit, as recorded by the sensitive detector. |
| `PD_Ed` | Energy deposited at the hit. |

The position, time, and energy quantities are stored in Geant4's internal units (typically mm, ns, and MeV). Momentum components and mass are explicitly converted to MeV-based units in the hit code. Events with zero accepted hits do not create `PD` tree entries, so the tree contains hit-bearing events rather than a complete event-by-event count. ROOT files are intended for later analysis and plotting; the `output/` directory contains analysis programs and examples.

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
    --export-gdml FILE Export geometry to GDML only (no physics or events)
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

### Export detector geometry to ROOT

With MG built and the Geant4/ROOT environment loaded, export a detector geometry to a ROOT `TGeoManager` file with:

```bash
./scripts/export_geometry_root.sh macros/EMMAtest0.dat output/EMMAtest0-geometry.root
```

The script exports Geant4 geometry to a temporary GDML file, imports it with ROOT, and writes the geometry to the requested `.root` file. The `--export-gdml` mode constructs the detector geometry and exits before creating the physics list, initializing physics, or running events. It can also be used directly with `./build/MG --geometry macros/EMMAtest0.dat --export-gdml output/EMMAtest0.gdml`. This is geometry-only; simulation hits remain in the regular MG output file and are not copied or overlaid.

## GitHub Actions

The workflow in `.github/workflows/build.yml` builds MG, runs a 100-event `EMMAtest0` smoke test, and renders the detector headlessly to `detector-visualization.pdf` using Geant4's ToolsSG off-screen driver. The PDF is uploaded as a workflow artifact. It uses a GitHub-hosted Ubuntu 24.04 runner and installs CVMFS to load `LCG_110/x86_64-ubuntu2404-gcc13-opt`, with the Geant4 data path above. It also installs the system X11/Xmu and OpenGL/GLU development packages required by Geant4's enabled visualization drivers, even for headless runs. The runner OS must match the LCG view's target OS; loading the Ubuntu 26.04 development view on Ubuntu 24.04 can fail because its binaries require newer system libraries such as OpenSSL.
