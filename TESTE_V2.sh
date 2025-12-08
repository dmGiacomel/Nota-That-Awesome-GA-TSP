#!/bin/bash

EXEC=./exec_parte_2.x
DIR_TSP=./arquivos_tsp

# ------------------------------------------
# Lista de arquivos para testar (*ordem importa*)
# ------------------------------------------
FILES=(
    # "eil51.tsp"
    # "kroA100.tsp"
    # "a280.tsp"
    # "rd400.tsp"
    # "d657.tsp"
    "d1291.tsp"
    "pcb3038.tsp"
)

# ------------------------------------------
# Mapa de grids por índice da lista FILES
# (0-based index)
# 0,1,2 → 2x2
# 3,4   → 3x3
# 5,6   → 4x4
# ------------------------------------------
GRID_ROWS=(4 4)
GRID_COLS=(4 4)

# ------------------------------------------
# Parâmetros do GA local
# ------------------------------------------
POP_L=350
MUT_L=0.20
CROSS_L=0.85
HEUR_L=0.05
ELIT_L=0.01
MAX_ITER_L=1700
TOUR_L=3

# ------------------------------------------
# Parâmetros do GA global
# ------------------------------------------
POP_G=700
MUT_G=0.15
CROSS_G=0.80
HEUR_G=0.05
ELIT_G=0.005
MAX_ITER_G=5000
TOUR_G=3

mkdir -p logs_2

# ------------------------------------------
# Execução
# ------------------------------------------

for IDX in "${!FILES[@]}"; do

    FILE="${FILES[$IDX]}"
    INPUT="$DIR_TSP/$FILE"
    LOG="logs_2/${FILE%.tsp}.log"

    # seleciona grid correto conforme índice
    GRID_R=${GRID_ROWS[$IDX]}
    GRID_C=${GRID_COLS[$IDX]}

    echo "== Testando $FILE =="
    echo "Grid usado: ${GRID_R}x${GRID_C}"
    echo "Gerando log em $LOG"

    {
        echo "Arquivo: $FILE"
        echo "Executado em: $(date)"
        echo "LOCAL:  pop=$POP_L mut=$MUT_L cross=$CROSS_L heur=$HEUR_L elit=$ELIT_L max_iter=$MAX_ITER_L tour=$TOUR_L"
        echo "GLOBAL: pop=$POP_G mut=$MUT_G cross=$CROSS_G heur=$HEUR_G elit=$ELIT_G max_iter=$MAX_ITER_G tour=$TOUR_G"
        echo "GRID:   rows=$GRID_R cols=$GRID_C"
        echo ""
        echo "seed,best_fitness,global_seconds,local_seconds,best_iter"
    } > "$LOG"

    for SEED in {1..5}; do
        echo -n "Rodando seed $SEED... "

        OUTPUT="$(
            $EXEC \
                "$INPUT" \
                "$POP_L" "$MUT_L" "$CROSS_L" "$HEUR_L" "$ELIT_L" "$MAX_ITER_L" "$TOUR_L" \
                "$POP_G" "$MUT_G" "$CROSS_G" "$HEUR_G" "$ELIT_G" "$MAX_ITER_G" "$TOUR_G" \
                "$GRID_R" "$GRID_C" \
                "$SEED"
        )"

        if [[ -z "$OUTPUT" ]]; then
            echo "$SEED,ERROR,ERROR,ERROR,ERROR" >> "$LOG"
            echo "ERRO"
            continue
        fi

        BF=$(echo "$OUTPUT" | cut -d',' -f1)
        TG=$(echo "$OUTPUT" | cut -d',' -f2)
        TL=$(echo "$OUTPUT" | cut -d',' -f3)
        BI=$(echo "$OUTPUT" | cut -d',' -f4)

        echo "$SEED,$BF,$TG,$TL,$BI" >> "$LOG"
        echo "ok"
    done

    echo ""
done

echo "Finalizado!"
