#!/bin/bash
# =============================================================================
# MASTER LAUNCHER: все 3 MEGA эксперимента параллельно
# Общий бюджет: ~832 CPU-часа из 1152 доступных (48ч × 24 ядра)
# Ожидаемое wall-time: ~30-36 часов
#
# Использование:
#   nohup bash run_all_mega.sh > mega_all.log 2>&1 &
#   tail -f mega_all.log
# =============================================================================
set -e

JOBS=20  # 20 параллельных jobs, 4 ядра в резерве

echo "============================================================"
echo "  MEGA EXPERIMENTS — Full Launch"
echo "  Budget: 48h × 24 cores = 1152 CPU-hours"
echo "  Planned: ~832 CPU-hours (72% utilization)"
echo "  Start: $(date)"
echo "============================================================"
echo ""
echo "Experiments:"
echo "  1. Phase Diagram MEGA  — ~570 CPU-h, 32 jobs"
echo "  2. Convergence MEGA    — ~170 CPU-h, 28 jobs"
echo "  3. Sensitivity MEGA    —  ~92 CPU-h, 13 jobs"
echo "  Total: 73 jobs across $JOBS cores"
echo ""

# Проверяем наличие бинарников
for BIN in ./levy_no_interrupt ./levy_teleport ./levy_viswanathan; do
    if [ ! -x "$BIN" ]; then
        echo "ERROR: Binary not found or not executable: $BIN"
        echo "Compile with: gcc -O3 -ffast-math -march=native -o ${BIN#./} ${BIN#./}.c -lm"
        exit 1
    fi
done
echo "All binaries found."

# Определяем метод параллелизации
if command -v parallel &> /dev/null; then
    USE_PARALLEL="gnu"
    echo "Using GNU parallel"
else
    USE_PARALLEL="xargs"
    echo "Using xargs -P (install 'parallel' for better progress tracking)"
fi
echo ""

# Создаём директории
mkdir -p phase_diagram_mega convergence_mega sensitivity_mega

# ─── Генерация ЕДИНОГО списка задач ───
TASKFILE=$(mktemp /tmp/mega_all_tasks.XXXXXX)

# ═══ Phase Diagram MEGA ═══
RV=0.5; MU_MIN=1.0; MU_MAX=3.0; MU_STEP=0.1; DX=50; SEED=42

declare -A GRID
GRID[10]="0.2     100000  1500  2000  200"
GRID[20]="0.1     200000  1500  2000  300"
GRID[50]="0.04    300000  1500  2000  500"
GRID[100]="0.02   300000  1500  2000  500"
GRID[200]="0.01   300000  1500  2000  500"
GRID[500]="0.004  300000  1500  2000  1000"
GRID[1000]="0.002 300000  1500  2000  1000"
GRID[2000]="0.001 300000  1500  2000  2000"

for lrv in 10 20 50 100 200 500 1000 2000; do
    read -r RHO FLIGHTS REPS_D REPS_ND CHUNK <<< "${GRID[$lrv]}"
    PREFIX="phase_diagram_mega/lrv_${lrv}"
    LC=$(echo "scale=4; 5*$RV" | bc -l)

    OUT="${PREFIX}_no_interrupt.csv"
    [ ! -s "$OUT" ] && echo "./levy_no_interrupt --rho $RHO --rv $RV --chunk-size $CHUNK --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP --flights $FLIGHTS --reps $REPS_D --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"

    OUT="${PREFIX}_teleport.csv"
    [ ! -s "$OUT" ] && echo "./levy_teleport --rho $RHO --rv $RV --chunk-size $CHUNK --lc $LC --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP --flights $FLIGHTS --reps $REPS_D --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"

    OUT="${PREFIX}_visw_destr.csv"
    [ ! -s "$OUT" ] && echo "./levy_viswanathan --rho $RHO --rv $RV --chunk-size $CHUNK --remove-prob 1.0 --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP --flights $FLIGHTS --reps $REPS_D --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"

    OUT="${PREFIX}_visw_nondestr.csv"
    [ ! -s "$OUT" ] && echo "./levy_viswanathan --rho $RHO --rv $RV --chunk-size $CHUNK --remove-prob 0.0 --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP --flights $FLIGHTS --reps $REPS_ND --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"
done

# ═══ Convergence MEGA ═══
for LRV in 100 1000; do
    if [ "$LRV" = "100" ]; then RHO=0.02; CHUNK=500; else RHO=0.002; CHUNK=1000; fi

    # Flights sweep (visw_destr, R=1000)
    for NF in 10000 50000 100000 200000 500000; do
        OUT="convergence_mega/lrv_${LRV}_flights_${NF}.csv"
        [ ! -s "$OUT" ] && echo "./levy_viswanathan --rho $RHO --rv $RV --chunk-size $CHUNK --remove-prob 1.0 --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP --flights $NF --reps 1000 --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"
    done

    # Reps sweep (visw_destr, F=200k)
    for NR in 50 200 500 1000 2000; do
        OUT="convergence_mega/lrv_${LRV}_reps_${NR}.csv"
        [ ! -s "$OUT" ] && echo "./levy_viswanathan --rho $RHO --rv $RV --chunk-size $CHUNK --remove-prob 1.0 --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP --flights 200000 --reps $NR --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"
    done

    # Non-destructive reps sweep (F=200k)
    for NR in 100 500 1000 2000; do
        OUT="convergence_mega/lrv_${LRV}_nd_reps_${NR}.csv"
        [ ! -s "$OUT" ] && echo "./levy_viswanathan --rho $RHO --rv $RV --chunk-size $CHUNK --remove-prob 0.0 --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP --flights 200000 --reps $NR --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"
    done
done

# ═══ Sensitivity MEGA ═══
RHO_S=0.02; LAMBDA_S=50; CHUNK_S=500

for RATIO in 0.1 0.2 0.5 1.0 2.0 5.0 10.0; do
    LMAX=$(echo "scale=4; $RATIO * $LAMBDA_S" | bc -l)
    OUT="sensitivity_mega/lmax_ratio_${RATIO}.csv"
    [ ! -s "$OUT" ] && echo "./levy_viswanathan --rho $RHO_S --rv $RV --chunk-size $CHUNK_S --remove-prob 1.0 --lmin $RV --lmax $LMAX --mu-min 1.0 --mu-max 3.0 --mu-step 0.1 --flights 200000 --reps 1000 --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"
done

for RATIO in 0.1 0.25 0.5 1.0 2.0 5.0; do
    LMIN=$(echo "scale=4; $RATIO * $RV" | bc -l)
    OUT="sensitivity_mega/lmin_ratio_${RATIO}.csv"
    [ ! -s "$OUT" ] && echo "./levy_viswanathan --rho $RHO_S --rv $RV --chunk-size $CHUNK_S --remove-prob 1.0 --lmin $LMIN --lmax $LAMBDA_S --mu-min 1.0 --mu-max 3.0 --mu-step 0.1 --flights 200000 --reps 1000 --seed $SEED --dx $DX --out $OUT" >> "$TASKFILE"
done

# ─── Статистика и запуск ───
NTASKS=$(wc -l < "$TASKFILE")

if [ "$NTASKS" -eq 0 ]; then
    echo "All experiments already completed!"
    rm "$TASKFILE"
    exit 0
fi

N_PHASE=$(grep -c 'phase_diagram_mega' "$TASKFILE" 2>/dev/null || echo 0)
N_CONV=$(grep -c 'convergence_mega' "$TASKFILE" 2>/dev/null || echo 0)
N_SENS=$(grep -c 'sensitivity_mega' "$TASKFILE" 2>/dev/null || echo 0)

echo "Tasks to run: $NTASKS"
echo "  Phase diagram: $N_PHASE / 32"
echo "  Convergence:   $N_CONV / 28"
echo "  Sensitivity:   $N_SENS / 13"
echo ""
echo "Starting in 10 seconds... (Ctrl+C to abort)"
sleep 10

START_TIME=$(date +%s)

# Функция для отображения прогресса
show_progress() {
    while true; do
        sleep 300  # каждые 5 минут
        DONE_P=$(ls phase_diagram_mega/*.csv 2>/dev/null | wc -l)
        DONE_C=$(ls convergence_mega/*.csv 2>/dev/null | wc -l)
        DONE_S=$(ls sensitivity_mega/*.csv 2>/dev/null | wc -l)
        TOTAL_DONE=$((DONE_P + DONE_C + DONE_S))
        NOW=$(date +%s)
        ELAPSED_MIN=$(( (NOW - START_TIME) / 60 ))
        echo "[${ELAPSED_MIN}min] Progress: phase=$DONE_P/32 conv=$DONE_C/28 sens=$DONE_S/13 total=$TOTAL_DONE/73"
    done
}

# Запускаем мониторинг в фоне
show_progress &
MONITOR_PID=$!
trap "kill $MONITOR_PID 2>/dev/null" EXIT

if [ "$USE_PARALLEL" = "gnu" ]; then
    parallel --jobs $JOBS \
        --progress \
        --joblog mega_joblog.txt \
        --timeout 180000 \
        < "$TASKFILE"
else
    # Fallback: xargs -P
    echo "=== Running with xargs -P $JOBS ==="
    cat "$TASKFILE" | xargs -I {} -P $JOBS bash -c '{} 2>/dev/null'
fi

kill $MONITOR_PID 2>/dev/null || true

END_TIME=$(date +%s)
ELAPSED_H=$(( (END_TIME - START_TIME) / 3600 ))
ELAPSED_M=$(( ((END_TIME - START_TIME) % 3600) / 60 ))

rm -f "$TASKFILE"

echo ""
echo "============================================================"
echo "  ALL MEGA EXPERIMENTS COMPLETE"
echo "  Wall time: ${ELAPSED_H}h ${ELAPSED_M}m"
echo "  Finished: $(date)"
echo "============================================================"
echo ""
echo "Results:"
echo "  Phase diagram: $(ls phase_diagram_mega/*.csv 2>/dev/null | wc -l) / 32 files"
echo "  Convergence:   $(ls convergence_mega/*.csv 2>/dev/null | wc -l) / 28 files"
echo "  Sensitivity:   $(ls sensitivity_mega/*.csv 2>/dev/null | wc -l) / 13 files"
echo ""
echo "Next steps:"
echo "  python3 plot_full_study.py --input-dir phase_diagram_mega"
