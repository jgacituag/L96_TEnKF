#!/usr/bin/env bash
# =============================================================================
# compile.sh  —  Build all Fortran modules as Python extensions via f2py
#
# IMPORTANT: activate the conda environment before running:
#   conda activate tenkf_l96
#   bash compile.sh
#
# Usage:
#   bash compile.sh            # gfortran + OpenMP
#   bash compile.sh --no-omp   # disable OpenMP (macOS or systems without libomp)
# =============================================================================
set -e
export SETUPTOOLS_USE_DISTUTILS=stdlib   # fix for setuptools >= 52 + f2py

OMP_FLAGS="-fopenmp -lgomp"
for arg in "$@"; do
  case $arg in --no-omp) OMP_FLAGS="" ;; esac
done
FFLAGS="-O3 ${OMP_FLAGS}"

ROOT="$(cd "$(dirname "$0")" && pwd)"
COMMON="${ROOT}/fortran/common"
DA="${ROOT}/fortran/da"
MODEL="${ROOT}/fortran/model"

# Use "python -m numpy.f2py" — works regardless of whether the f2py
# binary exists or has a broken shebang (common in conda environments)
PYTHON=$(which python)
F2PY="${PYTHON} -m numpy.f2py"

if [ -z "$PYTHON" ]; then
    echo "ERROR: python not found. Activate the conda environment first:"
    echo "  conda activate tenkf_l96"
    exit 1
fi

echo "============================================"
echo "  TEnKF L96 — Fortran compilation"
echo "  Python  : ${PYTHON}"
echo "  OMP     : ${OMP_FLAGS:-'(disabled)'}"
echo "============================================"

echo ""
echo "[1/3] Lorenz-96 model  →  fortran/model/model"
cd "${MODEL}"
${F2PY} -c --fcompiler=gnu95 --f90flags="${FFLAGS}" \
    "${COMMON}/common_tools.f90" \
    lorenzN.f90 \
    -m model
echo "      OK"

echo ""
echo "[2/3] Observation operator  →  fortran/da/obsope"
cd "${DA}"
${F2PY} -c --fcompiler=gnu95 --f90flags="${FFLAGS}" \
    "${COMMON}/netlib.f90" \
    "${COMMON}/common_tools.f90" \
    "${COMMON}/common_mtx.f90" \
    common_obs_lorenzN.f90 \
    -m obsope
echo "      OK"

echo ""
echo "[3/3] DA tools (LETKF + TEnKF)  →  fortran/da/da"
cd "${DA}"
${F2PY} -c --fcompiler=gnu95 --f90flags="${FFLAGS}" \
    "${COMMON}/netlib.f90" \
    "${COMMON}/common_tools.f90" \
    "${COMMON}/common_mtx.f90" \
    common_letkf.f90 \
    common_da_tools_1d.f90 \
    -m da
echo "      OK"

echo ""
echo "============================================"
echo "  Done. Run from repo root:"
echo "    python run_nature.py"
echo "    python run_example.py"
echo "============================================"
