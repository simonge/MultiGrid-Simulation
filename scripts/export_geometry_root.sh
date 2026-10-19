#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"

usage() {
  echo "Usage: $0 GEOMETRY.dat OUTPUT.root" >&2
  exit 2
}

[[ $# -eq 2 ]] || usage

geometry_file="$(realpath -- "$1")"
output_file="$(realpath -m -- "$2")"

if [[ ! -f "$geometry_file" ]]; then
  echo "Geometry input not found: $geometry_file" >&2
  exit 1
fi
if [[ ! -x "$repo_dir/build/MG" ]]; then
  echo "MG executable not found; build the project with ./setup.sh first." >&2
  exit 1
fi
if ! command -v root >/dev/null 2>&1; then
  echo "ROOT is not on PATH; source the Geant4/ROOT environment first." >&2
  exit 1
fi
if [[ -e "$output_file" ]]; then
  echo "Output already exists: $output_file" >&2
  exit 1
fi

output_dir="$(dirname -- "$output_file")"
mkdir -p -- "$output_dir"

work_dir="$(mktemp -d)"
trap 'rm -rf -- "$work_dir"' EXIT
gdml_file="$work_dir/geometry.gdml"
root_temp="$work_dir/geometry.root"

(
  cd "$work_dir"
  "$repo_dir/build/MG" --geometry "$geometry_file" --export-gdml "$gdml_file"
)

MG_GDML_INPUT="$gdml_file" MG_ROOT_OUTPUT="$root_temp" \
  root -l -b -q "$repo_dir/scripts/gdml_to_root.C"

if [[ ! -s "$root_temp" ]]; then
  echo "ROOT geometry export did not produce a non-empty file." >&2
  exit 1
fi
if [[ -e "$output_file" ]]; then
  echo "Output already exists: $output_file" >&2
  exit 1
fi
mv -- "$root_temp" "$output_file"

echo "Wrote ROOT geometry to $output_file"
