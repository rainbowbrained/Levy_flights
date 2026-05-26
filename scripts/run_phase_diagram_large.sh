#!/bin/bash
# Phase diagram: LARGE runs — 3-5× more flights, 2-3× more reps
set -e

RV=0.5
MU_MIN=1.0
MU_MAX=3.0
MU_STEP=0.1
DX=50
SEED=42

OUTDIR="phase_diagram_large"
mkdir -p "$OUTDIR"

# λ/rv -> "ρ  flights  reps  chunk_size"
declare -A GRID
GRID[10]="0.2     20000  100  200"
GRID[20]="0.1     20000  100  300"
GRID[50]="0.04    30000  100  500"
GRID[100]="0.02   50000  120  500"
GRID[200]="0.01   60000  120  500"
GRID[500]="0.004  80000  120  1000"
GRID[1000]="0.002 100000 120  1000"
GRID[2000]="0.001 150000 120  2000"

LRV_VALS=(10 20 50 100 200 500 1000 2000)

echo "=== Phase diagram LARGE: ${#LRV_VALS[@]} densities × 4 strategies ==="

for lrv in "${LRV_VALS[@]}"; do
    read -r RHO FLIGHTS REPS CHUNK <<< "${GRID[$lrv]}"
    LAMBDA=$(echo "scale=6; 1/(2*$RHO*$RV)" | bc -l)

    echo ""
    echo "--- λ/rv=$lrv  ρ=$RHO  λ=$LAMBDA  flights=$FLIGHTS reps=$REPS ---"

    PREFIX="${OUTDIR}/lrv_${lrv}"

    # 1) No interrupt
    OUT="${PREFIX}_no_interrupt.csv"
    if [ ! -f "$OUT" ]; then
        echo "  [1/4] No interrupt..."
        ./levy_no_interrupt \
            --rho $RHO --rv $RV --chunk-size $CHUNK \
            --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP \
            --flights $FLIGHTS --reps $REPS --seed $SEED --dx $DX \
            --out "$OUT" 2>&1 | tail -1
    else
        echo "  [1/4] skip (exists)"
    fi

    # 2) Teleport
    LC=$(echo "scale=4; 5*$RV" | bc -l)
    OUT="${PREFIX}_teleport.csv"
    if [ ! -f "$OUT" ]; then
        echo "  [2/4] Teleport lc=$LC..."
        ./levy_teleport \
            --rho $RHO --rv $RV --chunk-size $CHUNK --lc $LC \
            --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP \
            --flights $FLIGHTS --reps $REPS --seed $SEED --dx $DX \
            --out "$OUT" 2>&1 | tail -1
    else
        echo "  [2/4] skip (exists)"
    fi

    # 3) Viswanathan destructive
    OUT="${PREFIX}_visw_destr.csv"
    if [ ! -f "$OUT" ]; then
        echo "  [3/4] Viswanathan destructive..."
        ./levy_viswanathan \
            --rho $RHO --rv $RV --chunk-size $CHUNK --remove-prob 1.0 \
            --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP \
            --flights $FLIGHTS --reps $REPS --seed $SEED --dx $DX \
            --out "$OUT" 2>&1 | tail -1
    else
        echo "  [3/4] skip (exists)"
    fi

    # 4) Viswanathan non-destructive (3× reps)
    REPS_ND=$((REPS * 3))
    OUT="${PREFIX}_visw_nondestr.csv"
    if [ ! -f "$OUT" ]; then
        echo "  [4/4] Viswanathan non-destructive (reps=$REPS_ND)..."
        ./levy_viswanathan \
            --rho $RHO --rv $RV --chunk-size $CHUNK --remove-prob 0.0 \
            --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP \
            --flights $FLIGHTS --reps $REPS_ND --seed $SEED --dx $DX \
            --out "$OUT" 2>&1 | tail -1
    else
        echo "  [4/4] skip (exists)"
    fi
done

echo ""
echo "=== Phase diagram LARGE complete ==="
ls "$OUTDIR"/*.csv 2>/dev/null | wc -l
echo "CSV files."
