#!/usr/bin/env bash
#
# render.sh - Render PlantUML .puml diagrams to .svg on this machine.
#
# No installation required. This uses:
#   1. java            - already present at /usr/bin/java
#   2. plantuml jar    - the git-ignored jar in the PLUM repo root
#                        (plantuml-1.2026.6.jar); you run it, you don't install it
#   3. Smetana layout  - PlantUML's built-in pure-Java layout engine, so we do
#                        NOT need Graphviz/dot (which is not installed here).
#
# The Smetana engine is selected per-diagram by the line
#     !pragma layout smetana
# near the top of each .puml file. If a .puml lacks that pragma, PlantUML will
# try to call Graphviz and fail on this machine.
#
# Usage:
#   ./render.sh                       # render every *.puml in this directory
#   ./render.sh foo.puml bar.puml     # render only the named file(s)
#
# Output SVGs are written to the svg/ subdirectory next to this script.

set -euo pipefail

# Directory this script lives in (doc/uml), regardless of where it's called from.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Repo root is two levels up: doc/uml -> doc -> <repo root>.
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Locate the PlantUML jar. Prefer the one committed-convention location in the
# repo root; fall back to a shared copy under ~/apps if present.
JAR=""
for candidate in \
    "${REPO_ROOT}"/plantuml-*.jar \
    "${HOME}"/apps/plant_uml/plantuml-*.jar
do
    if [[ -f "${candidate}" ]]; then
        JAR="${candidate}"
        break
    fi
done

if [[ -z "${JAR}" ]]; then
    echo "ERROR: No PlantUML jar found." >&2
    echo "  Looked in: ${REPO_ROOT}/plantuml-*.jar and ${HOME}/apps/plant_uml/plantuml-*.jar" >&2
    echo "  Place a plantuml-<version>.jar in the repo root and re-run." >&2
    exit 1
fi

if ! command -v java >/dev/null 2>&1; then
    echo "ERROR: 'java' not found on PATH." >&2
    exit 1
fi

OUT_DIR="${SCRIPT_DIR}/svg"
mkdir -p "${OUT_DIR}"

# Which files to render: named args, or all *.puml in SCRIPT_DIR.
if [[ $# -gt 0 ]]; then
    FILES=("$@")
else
    FILES=("${SCRIPT_DIR}"/*.puml)
fi

echo "Using jar: ${JAR}"
echo "Output to: ${OUT_DIR}"
for f in "${FILES[@]}"; do
    echo "  rendering ${f##*/} ..."
    java -jar "${JAR}" -tsvg -o "${OUT_DIR}" "${f}"
done
echo "Done."
