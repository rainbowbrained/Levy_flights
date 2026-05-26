#!/bin/bash
# =============================================================================
# MEGA Convergence: сходимость по flights и reps
# Бюджет: ~170 CPU-часов
# Тестируем при 2 плотностях: lrv=100 (средняя) и lrv=1000 (разреженная)
# =============================================================================
set -e

RV=0.5
MU_MIN=1.0
MU_MAX=3.0
MU_STEP=0.1
DX=50
SEED=42
JOBS=20

OUTDIR="convergence_mega"
mkdir -p "$OUTDIR"

TASKFILE=$(mktemp /tmp/mega_conv_tasks.XXXXXX)

echo "=== Generating Convergence MEGA tasks ==="

# ─── A. Сходимость по flights (reps=1000 фикс., visw_destr) ───

for LRV in 100 1000; do
    if [ "$LRV" = "100" ]; then
        RHO=0.02; CHUNK=500
    else
        RHO=0.002; CHUNK=1000
    fi

    for NF in 10000 50000 100000 200000 500000; do
        OUT="${OUTDIR}/lrv_${LRV}_flights_${NF}.csv"
        if [ ! -s "$OUT" ]; then
            echo "./levy_viswanathan --rho $RHO --rv $RV --chunk-size $CHUNK --remove-prob 1.0 --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP --flights $NF --reps 1000 --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"
        fi
    done
done

# ─── B. Сходимость по reps (F=200k фикс., visw_destr) ───

for LRV in 100 1000; do
    if [ "$LRV" = "100" ]; then
        RHO=0.02; CHUNK=500
    else
        RHO=0.002; CHUNK=1000
    fi

    for NR in 50 200 500 1000 2000; do
        OUT="${OUTDIR}/lrv_${LRV}_reps_${NR}.csv"
        if [ ! -s "$OUT" ]; then
            echo "./levy_viswanathan --rho $RHO --rv $RV --chunk-size $CHUNK --remove-prob 1.0 --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP --flights 200000 --reps $NR --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"
        fi
    done
done

# ─── C. Сходимость non-destructive по reps (F=200k) ───

for LRV in 100 1000; do
    if [ "$LRV" = "100" ]; then
        RHO=0.02; CHUNK=500
    else
        RHO=0.002; CHUNK=1000
    fi

    for NR in 100 500 1000 2000; do
        OUT="${OUTDIR}/lrv_${LRV}_nd_reps_${NR}.csv"
        if [ ! -s "$OUT" ]; then
            echo "./levy_viswanathan --rho $RHO --rv $RV --chunk-size $CHUNK --remove-prob 0.0 --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP --flights 200000 --reps $NR --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"
        fi
    done
done

NTASKS=$(wc -l < "$TASKFILE")
echo "=== $NTASKS convergence tasks ==="

if [ "$NTASKS" -eq 0 ]; then
    echo "All convergence experiments already completed!"
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
echo "=== Convergence MEGA complete ==="
ls "$OUTDIR"/*.csv 2>/dev/null | wc -l
echo "CSV files."
