#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
echo "========================================"
echo " Multi-Gaussian Fit GUI"
echo "========================================"
echo ""

JULIA=""
for candidate in julia /usr/local/bin/julia; do
    if command -v "$candidate" &>/dev/null 2>&1; then
        JULIA="$candidate"; break
    fi
done
if [ -z "$JULIA" ]; then
    for dir in /opt/julia-* ~/julia-* /Applications/Julia-*; do
        if [ -x "$dir/bin/julia" ]; then JULIA="$dir/bin/julia"; break; fi
    done
fi
if [ -z "$JULIA" ]; then
    echo "ERROR: Julia not found. Install Julia 1.10+ from https://julialang.org/downloads/"
    exit 1
fi

echo "  Julia: $JULIA"
echo "  Installing dependencies..."

cd "$SCRIPT_DIR"

# Install the core package if not already available
$JULIA -e 'using Pkg; Pkg.develop(path="'$HOME'/Git/GaussianFit1D.jl"); Pkg.instantiate()' 2>/dev/null || true
$JULIA --project=. -e 'using Pkg; Pkg.develop(path="'$HOME'/Git/GaussianFit1D.jl"); Pkg.instantiate()'

echo "  Starting GUI..."
$JULIA --project=. app.jl
