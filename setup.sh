#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
venv_dir="${repo_dir}/.venv"
build_dir="${repo_dir}/build"

# LCG may export these for its own Python; do not let them override the
# interpreter or site configuration used by the repository-local venv.
unset PYTHONHOME PYTHONPATH

if [[ -n "${CONDA_PREFIX:-}" ]]; then
  echo "Deactivate the Conda environment first; this setup uses system Geant4/ROOT and a repo-local Python venv." >&2
  exit 1
fi

python_executable="${MG_PYTHON:-$(command -v python3 || true)}"
if [[ -z "${python_executable}" || ! -x "${python_executable}" ]]; then
  echo "A usable Python interpreter is required to create the local NCrystalGeant4 environment." >&2
  exit 1
fi
if ! "${python_executable}" -c 'import ssl' >/dev/null 2>&1; then
  echo "Python interpreter '${python_executable}' has no working ssl module; set MG_PYTHON to an SSL-enabled Python." >&2
  exit 1
fi

if [[ ! -x "${venv_dir}/bin/python" ]]; then
  "${python_executable}" -m venv "${venv_dir}"
fi
if ! "${venv_dir}/bin/python" -c 'import ssl' >/dev/null 2>&1; then
  echo "The existing virtual environment '${venv_dir}' does not have a working ssl module; remove it and rerun setup with SSL-enabled MG_PYTHON." >&2
  exit 1
fi

# Keep the pip-installed bridge and its helper command local to this checkout.
"${venv_dir}/bin/python" -m pip install -r "${repo_dir}/requirements.txt"
export PATH="${venv_dir}/bin:${PATH}"
ncrystal_geant4_cmake_dir="$(ncrystal-geant4 --cmakedir)"
if [[ -z "${ncrystal_geant4_cmake_dir}" ]]; then
  echo "ncrystal-geant4 did not report its CMake package directory." >&2
  exit 1
fi

cmake -S "${repo_dir}" -B "${build_dir}" \
  -DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}" \
  -DNCrystalGeant4_DIR="${ncrystal_geant4_cmake_dir}"
cmake --build "${build_dir}" --parallel

echo "Build complete: ${build_dir}/MG"
