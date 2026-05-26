#!/bin/bash
# Convergence analysis: how μ* and λη stabilize with flights and reps
# Fixed: ρ=0.02, rv=0.5, λ/rv=100, destructive Viswanathan

set -e

RHO=0.02
RV=0.5
SEED=42
DX=50

OUTDIR="convergence"
mkdir -p "$OUTDIR"

echo "=== Convergence in flights (reps=50 fixed) ==="
REPS=50
for NF in 500 1000 2000 5000 10000 20000 50000 100000; do
    OUT="${OUTDIR}/flights_${NF}.csv"
    if [ ! -f "$OUT" ]; then
        echo "  flights=$NF..."
        ./levy_viswanathan \
            --rho $RHO --rv $RV --remove-prob 1.0 \
            --mu-min 1.0 --mu-max 3.0 --mu-step 0.1 \
            --flights $NF --reps $REPS --seed $SEED --dx $DX \
            --out "$OUT" 2>&1 | tail -1
    fi
done

echo ""
echo "=== Convergence in reps (flights=20000 fixed) ==="
NF=20000
for NR in 5 10 20 50 100 200; do
    OUT="${OUTDIR}/reps_${NR}.csv"
    if [ ! -f "$OUT" ]; then
        echo "  reps=$NR..."
        ./levy_viswanathan \
            --rho $RHO --rv $RV --remove-prob 1.0 \
            --mu-min 1.0 --mu-max 3.0 --mu-step 0.1 \
            --flights $NF --reps $NR --seed $SEED --dx $DX \
            --out "$OUT" 2>&1 | tail -1
    fi
done

echo ""
echo "=== Convergence in reps for NON-DESTRUCTIVE (flights=20000) ==="
for NR in 10 30 50 100 200; do
    OUT="${OUTDIR}/nd_reps_${NR}.csv"
    if [ ! -f "$OUT" ]; then
        echo "  non-destr reps=$NR..."
        ./levy_viswanathan \
            --rho $RHO --rv $RV --remove-prob 0.0 \
            --mu-min 1.0 --mu-max 3.0 --mu-step 0.1 \
            --flights $NF --reps $NR --seed $SEED --dx $DX \
            --out "$OUT" 2>&1 | tail -1
    fi
done

echo "=== Convergence complete ==="
