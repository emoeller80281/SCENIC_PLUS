#!/bin/bash -l
#SBATCH --job-name="submit_multiple_scenic_plus_jobs"
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH -c 1
#SBATCH --mem=4G

SCRIPT_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS"

MAX_JOBS_IN_QUEUE=25

submit_run_scenic_plus_job() {
    local SAMPLE_NAME=$1
    local CELL_TYPE=$2
    local SPECIES=$3
    local INPUT_DIR=$4
    local RNA_FILE_NAME=$5
    local ATAC_FILE_NAME=$6

    # Ensure the log directory exists
    mkdir -p "LOGS/${CELL_TYPE}_logs/${SAMPLE_NAME}_logs"

    # Submit the job
    sbatch \
        --export=ALL,SAMPLE_NAME="$SAMPLE_NAME",CELL_TYPE="$CELL_TYPE",SPECIES="$SPECIES",INPUT_DIR="$INPUT_DIR",RNA_FILE_NAME="$RNA_FILE_NAME",ATAC_FILE_NAME="$ATAC_FILE_NAME",SCRIPT_DIR="$SCRIPT_DIR" \
        --output="LOGS/${CELL_TYPE}_logs/${SAMPLE_NAME}_logs/scenic_plus_${CELL_TYPE}_${SAMPLE_NAME}.out" \
        --error="LOGS/${CELL_TYPE}_logs/${SAMPLE_NAME}_logs/scenic_plus_${CELL_TYPE}_${SAMPLE_NAME}.err" \
        --job-name="SCENIC+_${CELL_TYPE}_${SAMPLE_NAME}" \
        "${SCRIPT_DIR}/run_scenic_plus.sh"
}

run_macrophage() {
    local CELL_TYPE="macrophage"
    local SAMPLE_NAMES=(
        # "muon_buffer_1"
        "muon_buffer_2"
        "muon_buffer_3"
        # "muon_buffer_4"
        )
    local SPECIES="human"

    for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"; do

        local INPUT_DIR="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/LINGER/LINGER_MACROPHAGE/muon_${SAMPLE_NAME}"
        # Check how many jobs are currently queued/running
        while [ "$(squeue -u $USER | grep SCENIC+ | wc -l)" -ge "$MAX_JOBS_IN_QUEUE" ]; do
            echo "[INFO] Maximum jobs ($MAX_JOBS_IN_QUEUE) in queue. Waiting 60 seconds..."
            sleep 60
        done

        local RNA_FILE_NAME="${SAMPLE_NAME}_RNA.csv"
        local ATAC_FILE_NAME="${SAMPLE_NAME}_ATAC.csv"

        # Submit the job for each sample
        submit_run_scenic_plus_job \
            "$SAMPLE_NAME" \
            "$CELL_TYPE" \
            "$SPECIES" \
            "$INPUT_DIR" \
            "$RNA_FILE_NAME" \
            "$ATAC_FILE_NAME"
    done
}

run_mESC(){
    local CELL_TYPE="mESC"

    local SAMPLE_NAMES=(
        # "muon_E7.5_rep1"
        # "muon_E7.5_rep2"
        # "muon_E8.5_rep1"
        "muon_E8.5_rep2"
    )
    local SPECIES="mouse"

    

    # Submit each SAMPLE_NAME as a separate job
    for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"; do
        # Check how many jobs are currently queued/running
        while [ "$(squeue -u $USER | grep SCENIC+ | wc -l)" -ge "$MAX_JOBS_IN_QUEUE" ]; do
            echo "[INFO] Maximum jobs ($MAX_JOBS_IN_QUEUE) in queue. Waiting 60 seconds..."
            sleep 60
        done

        local INPUT_DIR="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/LINGER/LINGER_MESC_SC_DATA/FULL_MESC_SAMPLES/muon_${SAMPLE_NAME}"

        local RNA_FILE_NAME="${SAMPLE_NAME}_RNA.csv"
        local ATAC_FILE_NAME="${SAMPLE_NAME}_ATAC.csv"

        # Submit the job for each sample
        submit_run_scenic_plus_job \
            "$SAMPLE_NAME" \
            "$CELL_TYPE" \
            "$SPECIES" \
            "$INPUT_DIR" \
            "$RNA_FILE_NAME" \
            "$ATAC_FILE_NAME"

    done
}

run_K562(){
    local CELL_TYPE="K562"
    local SAMPLE_NAMES=(
        # "muon_sample_1"
    )
    local SPECIES="human"

    

    echo "Submitting SCENIC+ jobs"
    # Submit each SAMPLE_NAME as a separate job
    for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"; do
        # Check how many jobs are currently queued/running
        while [ "$(squeue -u $USER | grep SCENIC+ | wc -l)" -ge "$MAX_JOBS_IN_QUEUE" ]; do
            echo "Maximum jobs ($MAX_JOBS_IN_QUEUE) in queue. Checking again in 5 minutes..."
            sleep 300
        done

        local INPUT_DIR="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/LINGER/LINGER_K562/muon_${SAMPLE_NAME}"

        local RNA_FILE_NAME="${SAMPLE_NAME}_RNA.csv"
        local ATAC_FILE_NAME="${SAMPLE_NAME}_ATAC.csv"

        # Submit the job for each sample
        submit_run_scenic_plus_job \
            "$SAMPLE_NAME" \
            "$CELL_TYPE" \
            "$SPECIES" \
            "$INPUT_DIR" \
            "$RNA_FILE_NAME" \
            "$ATAC_FILE_NAME"

        echo "  - Running SCENIC+ for ${SAMPLE_NAME}"
    done
}

run_iPSC(){
    local CELL_TYPE="iPSC"
    local SAMPLE_NAMES=(
        "muon_WT_D13_rep1"
    )
    local SPECIES="human"

    

    echo "Submitting SCENIC+ jobs"
    # Submit each SAMPLE_NAME as a separate job
    for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"; do
        # Check how many jobs are currently queued/running
        while [ "$(squeue -u $USER | grep SCENIC+ | wc -l)" -ge "$MAX_JOBS_IN_QUEUE" ]; do
            echo "Maximum jobs ($MAX_JOBS_IN_QUEUE) in queue. Checking again in 5 minutes..."
            sleep 300
        done

        local INPUT_DIR="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/LINGER/IPS_CELL_DATA/muon_${SAMPLE_NAME}"
        

        local RNA_FILE_NAME="${SAMPLE_NAME}_RNA.csv"
        local ATAC_FILE_NAME="${SAMPLE_NAME}_ATAC.csv"

        # Submit the job for each sample
        submit_run_scenic_plus_job \
            "$SAMPLE_NAME" \
            "$CELL_TYPE" \
            "$SPECIES" \
            "$INPUT_DIR" \
            "$RNA_FILE_NAME" \
            "$ATAC_FILE_NAME"

        echo "  - Running SCENIC+ for ${SAMPLE_NAME}"
    done
}

# run_K562
run_macrophage
run_mESC
# run_iPSC