#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
venv_dir="${repo_dir}/.venv"
build_dir="${repo_dir}/build"

if [[ -n "${CONDA_PREFIX:-}" ]]; then
  echo "Deactivate the Conda environment first; this setup uses system Geant4/ROOT and a repo-local Python venv." >&2
  exit 1
fi

if ! command -v python3 >/dev/null 2>&1; then
  echo "python3 is required to create the local NCrystalGeant4 environment." >&2
  exit 1
fi

if [[ ! -x "${venv_dir}/bin/python" ]]; then
  python3 -m venv "${venv_dir}"
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
