#!/bin/bash
# =============================================================================
# MEGA Sensitivity: чувствительность к lmin/lmax
# Бюджет: ~92 CPU-часа
# При ρ=0.02, rv=0.5, λ=50, visw_destr
# =============================================================================
set -e

RHO=0.02
RV=0.5
LAMBDA=50
FLIGHTS=200000
REPS=1000
SEED=42
DX=50
CHUNK=500
JOBS=20

OUTDIR="sensitivity_mega"
mkdir -p "$OUTDIR"

TASKFILE=$(mktemp /tmp/mega_sens_tasks.XXXXXX)

echo "=== Generating Sensitivity MEGA tasks ==="

# ─── lmax sensitivity (lmin=rv=0.5 фикс.) ───

for RATIO in 0.1 0.2 0.5 1.0 2.0 5.0 10.0; do
    LMAX=$(echo "scale=4; $RATIO * $LAMBDA" | bc -l)
    OUT="${OUTDIR}/lmax_ratio_${RATIO}.csv"
    if [ ! -s "$OUT" ]; then
        echo "./levy_viswanathan --rho $RHO --rv $RV --chunk-size $CHUNK --remove-prob 1.0 --lmin $RV --lmax $LMAX --mu-min 1.0 --mu-max 3.0 --mu-step 0.1 --flights $FLIGHTS --reps $REPS --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"
    fi
done

# ─── lmin sensitivity (lmax=λ=50 фикс.) ───

for RATIO in 0.1 0.25 0.5 1.0 2.0 5.0; do
    LMIN=$(echo "scale=4; $RATIO * $RV" | bc -l)
    OUT="${OUTDIR}/lmin_ratio_${RATIO}.csv"
    if [ ! -s "$OUT" ]; then
        echo "./levy_viswanathan --rho $RHO --rv $RV --chunk-size $CHUNK --remove-prob 1.0 --lmin $LMIN --lmax $LAMBDA --mu-min 1.0 --mu-max 3.0 --mu-step 0.1 --flights $FLIGHTS --reps $REPS --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"
    fi
done

NTASKS=$(wc -l < "$TASKFILE")
echo "=== $NTASKS sensitivity tasks ==="

if [ "$NTASKS" -eq 0 ]; then
    echo "All sensitivity experiments already completed!"
    rm "$TASKFILE"
    exit 0
fi

if command -v parallel &> /dev/null; then
    echo "=== Running with GNU parallel (jobs=$JOBS) ==="
    parallel --jobs $JOBS --progress --joblog "${OUTDIR}/joblog.txt" < "$TASKFILE"
else
    echo "=== Running sequentially ==="
    while IFS= read -r cmd; do
        echo "Running: $cmd"
        eval "$cmd" 2>&1 | tail -1
    done < "$TASKFILE"
fi

rm "$TASKFILE"

echo ""
echo "=== Sensitivity MEGA complete ==="
ls "$OUTDIR"/*.csv 2>/dev/null | wc -l
echo "CSV files."
