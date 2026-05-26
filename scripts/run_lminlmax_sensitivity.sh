#!/bin/bash
# Sensitivity to lmin/lmax truncation
# Baseline: ρ=0.02, rv=0.5, λ=50, λ/rv=100
# Vary lmax/λ and lmin/rv

set -e

RHO=0.02
RV=0.5
LAMBDA=50
FLIGHTS=15000
REPS=50
SEED=42
DX=50

OUTDIR="sensitivity"
mkdir -p "$OUTDIR"

echo "=== lmax sensitivity (lmin=rv=$RV fixed) ==="
# lmax/λ = 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0
for RATIO in 0.1 0.2 0.5 1.0 2.0 5.0 10.0; do
    LMAX=$(echo "scale=4; $RATIO * $LAMBDA" | bc -l)
    OUT="${OUTDIR}/lmax_ratio_${RATIO}.csv"
    if [ ! -f "$OUT" ]; then
        echo "  lmax/λ=$RATIO (lmax=$LMAX)..."
        # Use destructive Viswanathan for clean signal
        ./levy_viswanathan \
            --rho $RHO --rv $RV --remove-prob 1.0 \
            --lmin $RV --lmax $LMAX \
            --mu-min 1.0 --mu-max 3.0 --mu-step 0.1 \
            --flights $FLIGHTS --reps $REPS --seed $SEED --dx $DX \
            --out "$OUT" 2>&1 | tail -1
    fi
done

echo ""
echo "=== lmin sensitivity (lmax=λ=$LAMBDA fixed) ==="
# lmin/rv = 0.1, 0.25, 0.5, 1.0, 2.0, 5.0
for RATIO in 0.1 0.25 0.5 1.0 2.0 5.0; do
    LMIN=$(echo "scale=4; $RATIO * $RV" | bc -l)
    OUT="${OUTDIR}/lmin_ratio_${RATIO}.csv"
    if [ ! -f "$OUT" ]; then
        echo "  lmin/rv=$RATIO (lmin=$LMIN)..."
        ./levy_viswanathan \
            --rho $RHO --rv $RV --remove-prob 1.0 \
            --lmin $LMIN --lmax $LAMBDA \
            --mu-min 1.0 --mu-max 3.0 --mu-step 0.1 \
            --flights $FLIGHTS --reps $REPS --seed $SEED --dx $DX \
            --out "$OUT" 2>&1 | tail -1
    fi
done

echo ""
echo "=== Sensitivity runs complete ==="
ls "$OUTDIR"/*.csv | wc -l
echo "CSV files."
