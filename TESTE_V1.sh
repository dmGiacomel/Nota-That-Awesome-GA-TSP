#!/bin/bash

EXEC=./exec_parte_1.x
DIR_TSP=./arquivos_tsp

# ------------------------------------------
# Lista de arquivos para testar (EDITE AQUI)
# ------------------------------------------
FILES=(
    "eil51.tsp"
    "kroA100.tsp"
    "a280.tsp"
    "rd400.tsp"
    "d657.tsp"
    "d1291.tsp"
    "pcb3038.tsp"
    "rl5915.tsp"
    # adicione/remova aqui
)

# ------------------------------------------
# Parâmetros fixos do GA (EDITE SE QUISER)
# ------------------------------------------
POP=700
MUT=0.15
CROSS=0.80
HEUR=0.05
ELIT=0.005
MAX_ITER=5000
TOUR=3
# semente será sobrescrita no loop

# ------------------------------------------
# Execução
# ------------------------------------------

mkdir -p logs

for FILE in "${FILES[@]}"; do
    INPUT="$DIR_TSP/$FILE"
    LOG="logs/${FILE%.tsp}.log"

    echo "== Testando $FILE ==" 
    echo "Gerando log em $LOG"
    echo "Arquivo: $FILE" > "$LOG"
    echo "Executado em: $(date)" >> "$LOG"
    echo "pop=$POP mut=$MUT cross=$CROSS heur=$HEUR elit=$ELIT max_iter=$MAX_ITER tour=$TOUR" >> "$LOG"
    echo "" >> "$LOG"
    echo "seed,fitness,seconds,best_iter" >> "$LOG"

    for SEED in {1..10}; do
        echo -n "Rodando seed $SEED... "

        OUTPUT=$($EXEC "$INPUT" "$POP" "$MUT" "$CROSS" "$HEUR" "$ELIT" "$MAX_ITER" "$SEED" "$TOUR")

        # OUTPUT = "fitness,seconds"
        FITNESS=$(echo "$OUTPUT" | cut -d',' -f1)
        EXEC_SECONDS=$(echo "$OUTPUT" | cut -d',' -f2)
        BEST_ITER=$(echo "$OUTPUT" | cut -d',' -f3)

        echo "$SEED,$FITNESS,$EXEC_SECONDS,$BEST_ITER" >> "$LOG"

        echo "ok"
    done

    echo ""
done

echo "Finalizado!"