#!/bin/bash
# =============================================================================
# MEGA Phase Diagram: масштабные эксперименты для гладких графиков
# Бюджет: ~570 CPU-часов, ~30 часов wall-time на 24 ядрах
# =============================================================================
set -e

RV=0.5
MU_MIN=1.0
MU_MAX=3.0
MU_STEP=0.1
DX=50
SEED=42
JOBS=20  # Оставляем 4 ядра свободными для системы и других задач

OUTDIR="phase_diagram_mega"
mkdir -p "$OUTDIR"

# λ/rv -> "ρ  flights  reps_destr  reps_nondestr  chunk_size"
declare -A GRID
GRID[10]="0.2     100000  1500  2000  200"
GRID[20]="0.1     200000  1500  2000  300"
GRID[50]="0.04    300000  1500  2000  500"
GRID[100]="0.02   300000  1500  2000  500"
GRID[200]="0.01   300000  1500  2000  500"
GRID[500]="0.004  300000  1500  2000  1000"
GRID[1000]="0.002 300000  1500  2000  1000"
GRID[2000]="0.001 300000  1500  2000  2000"

LRV_VALS=(10 20 50 100 200 500 1000 2000)

# Генерация списка задач для GNU parallel
TASKFILE=$(mktemp /tmp/mega_tasks.XXXXXX)

echo "=== Generating task list for Phase Diagram MEGA ==="

for lrv in "${LRV_VALS[@]}"; do
    read -r RHO FLIGHTS REPS_D REPS_ND CHUNK <<< "${GRID[$lrv]}"

    PREFIX="${OUTDIR}/lrv_${lrv}"

    # 1) No interrupt
    OUT="${PREFIX}_no_interrupt.csv"
    if [ ! -s "$OUT" ]; then
        echo "./levy_no_interrupt --rho $RHO --rv $RV --chunk-size $CHUNK --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP --flights $FLIGHTS --reps $REPS_D --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"
    fi

    # 2) Teleport (lc = 5*rv)
    LC=$(echo "scale=4; 5*$RV" | bc -l)
    OUT="${PREFIX}_teleport.csv"
    if [ ! -s "$OUT" ]; then
        echo "./levy_teleport --rho $RHO --rv $RV --chunk-size $CHUNK --lc $LC --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP --flights $FLIGHTS --reps $REPS_D --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"
    fi

    # 3) Viswanathan destructive
    OUT="${PREFIX}_visw_destr.csv"
    if [ ! -s "$OUT" ]; then
        echo "./levy_viswanathan --rho $RHO --rv $RV --chunk-size $CHUNK --remove-prob 1.0 --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP --flights $FLIGHTS --reps $REPS_D --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"
    fi

    # 4) Viswanathan non-destructive
    OUT="${PREFIX}_visw_nondestr.csv"
    if [ ! -s "$OUT" ]; then
        echo "./levy_viswanathan --rho $RHO --rv $RV --chunk-size $CHUNK --remove-prob 0.0 --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP --flights $FLIGHTS --reps $REPS_ND --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"
    fi
done

NTASKS=$(wc -l < "$TASKFILE")
echo "=== $NTASKS tasks to run (skipped already completed) ==="

if [ "$NTASKS" -eq 0 ]; then
    echo "All experiments already completed!"
    rm "$TASKFILE"
    exit 0
fi

echo ""
echo "Estimated total CPU time: ~570 hours"
echo "Estimated wall time with $JOBS parallel jobs: ~30 hours"
echo ""
echo "Starting in 5 seconds... (Ctrl+C to abort)"
sleep 5

# Запуск через GNU parallel
if command -v parallel &> /dev/null; then
    echo "=== Running with GNU parallel (jobs=$JOBS) ==="
    parallel --jobs $JOBS --progress --joblog "${OUTDIR}/joblog.txt" < "$TASKFILE"
else
    echo "=== GNU parallel not found, running sequentially ==="
    echo "WARNING: Sequential execution will take ~570 hours!"
    echo "Install GNU parallel: sudo apt install parallel"
    echo ""
    while IFS= read -r cmd; do
        echo "Running: $cmd"
        eval "$cmd" 2>&1 | tail -1
    done < "$TASKFILE"
fi

rm "$TASKFILE"

echo ""
echo "=== Phase Diagram MEGA complete ==="
echo "Results in: $OUTDIR/"
ls "$OUTDIR"/*.csv 2>/dev/null | wc -l
echo "CSV files generated."
