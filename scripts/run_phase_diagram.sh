#!/bin/bash
# Phase diagram: systematic grid of experiments
# λ/rv = 2/ρ (for rv=0.5), so ρ = 2/(λ/rv)
#
# λ/rv values: 10, 20, 50, 100, 200, 500, 1000, 2000
# Strategies: no_interrupt, teleport, viswanathan_destructive, viswanathan_nondestructive

set -e

RV=0.5
MU_MIN=1.0
MU_MAX=3.0
MU_STEP=0.1
DX=50
SEED=42

OUTDIR="phase_diagram"
mkdir -p "$OUTDIR"

# λ/rv -> ρ mapping (rv=0.5, λ=1/(2ρrv), so λ/rv = 1/(ρ rv²) => but actually
# λ = 1/(2ρrv) => λ/rv = 1/(2ρrv²) = 1/(2ρ×0.25) = 2/ρ => ρ = 2/(λ/rv))
declare -A GRID
# λ/rv  ->  "ρ  flights  reps  chunk_size"
GRID[10]="0.2     5000  40  200"
GRID[20]="0.1     5000  40  300"
GRID[50]="0.04   10000  40  500"
GRID[100]="0.02  15000  50  500"
GRID[200]="0.01  20000  50  500"
GRID[500]="0.004 30000  50  1000"
GRID[1000]="0.002 40000 50  1000"
GRID[2000]="0.001 50000 50  2000"

# Sorted λ/rv values
LRV_VALS=(10 20 50 100 200 500 1000 2000)

echo "=== Phase diagram: ${#LRV_VALS[@]} density points × 4 strategies ==="

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
        echo "  [1/4] No interrupt — already exists, skip"
    fi

    # 2) Teleport (lc = 5*rv)
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
        echo "  [2/4] Teleport — already exists, skip"
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
        echo "  [3/4] Visw. destructive — already exists, skip"
    fi

    # 4) Viswanathan non-destructive (more reps)
    REPS_ND=$((REPS * 2))
    OUT="${PREFIX}_visw_nondestr.csv"
    if [ ! -f "$OUT" ]; then
        echo "  [4/4] Viswanathan non-destructive (reps=$REPS_ND)..."
        ./levy_viswanathan \
            --rho $RHO --rv $RV --chunk-size $CHUNK --remove-prob 0.0 \
            --mu-min $MU_MIN --mu-max $MU_MAX --mu-step $MU_STEP \
            --flights $FLIGHTS --reps $REPS_ND --seed $SEED --dx $DX \
            --out "$OUT" 2>&1 | tail -1
    else
        echo "  [4/4] Visw. non-destructive — already exists, skip"
    fi
done

echo ""
echo "=== Phase diagram complete ==="
echo "Results in: $OUTDIR/"
ls -la "$OUTDIR"/*.csv 2>/dev/null | wc -l
echo "CSV files generated."
