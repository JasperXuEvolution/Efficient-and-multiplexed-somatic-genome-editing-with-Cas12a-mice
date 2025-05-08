#!/bin/bash
# 01-info_extraction_cas12a_quadruple_guide.bash
# This shell script is part of the UltraSeq pipeline for processing Cas12a quadruple guide data.
# It performs sequence merging, gRNA/barcode extraction, clustering, and data aggregation.
# Additionally, it cleans up intermediate files and outputs a summary log for quality control.

# ---------------------------------------------------------------------
# SBATCH directives: Job submission options for the SLURM scheduler.
# ---------------------------------------------------------------------
# The interpreter used to execute the script

# SBATCH directives that convey submission options:
#SBATCH --job-name=01-info_extraction_cas12a_quadruple_guide
#SBATCH --mail-user=xxx
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50g 
#SBATCH --time=24:00:00
#SBATCH --account=xxx
#SBATCH --partition=batch

# ---------------------------------------------------------------------
# Environment Setup: Load configuration, modules, and activate Conda environment.
# ---------------------------------------------------------------------
# Load required modules and activate the environment
source ../config.sh
module load adapterremoval/2.3.1
source ~/miniconda3/etc/profile.d/conda.sh
conda activate UltraSeq || { echo "Error: Failed to activate Conda environment"; exit 1; }

# ---------------------------------------------------------------------
# Directories and Input Files Setup
# ---------------------------------------------------------------------
# Define directories and input paths
working_dir="$PROJECT_DIR/01_data_collection"
Python_script_address="${working_dir}/main_code"
guide_ref="$working_dir/data/guide_reference-cas12a_efficiency.csv"
Input_data_info_address="${working_dir}/data/NGS_address"
Step1_address="${working_dir}/data/Bartender"
Step3_address="${working_dir}/data/Processed_data"

# Ensure the input data information file exists
if [[ ! -f "${Input_data_info_address}" ]]; then
    echo "Error: Input data information file not found: ${Input_data_info_address}"
    exit 1
fi

# Create main directories
mkdir -p "${Step1_address}" || { echo "Error: Failed to create directory ${Step1_address}"; exit 1; }
mkdir -p "${Step3_address}" || { echo "Error: Failed to create directory ${Step3_address}"; exit 1; }

# ---------------------------------------------------------------------
# Main Pipeline: Process each sample listed in the input data info file.
# ---------------------------------------------------------------------
while read -r line; do
    r1=$(echo "$line" | cut -d',' -f1)
    r2=$(echo "$line" | cut -d',' -f2)
    sampleID=$(echo "$line" | cut -d',' -f3)

    # Check if R1 and R2 files exist
    if [[ ! -f "${r1}" || ! -f "${r2}" ]]; then
        echo "Warning: Missing input files for sample ${sampleID} (R1: ${r1}, R2: ${r2}). Skipping..."
        continue
    fi

    # Create sample-specific directories
    sample_folder="${Step1_address}/${sampleID}"
    mkdir -p "${sample_folder}" || { echo "Error: Failed to create directory ${sample_folder}"; exit 1; }

    clonal_barcode_folder="${sample_folder}/Clonal_barcode"
    mkdir -p "${clonal_barcode_folder}" || { echo "Error: Failed to create directory ${clonal_barcode_folder}"; exit 1; }

    # Step 1: Extract gRNA, barcode info
    python3 "${Python_script_address}/cas12a_quadruple_guide_parsing.py" \
        --a1 "${r1}" --a2 "${r2}" --b "${guide_ref}" \
        --o "${sample_folder}" || { echo "Error: Step 1 failed for sample ${sampleID}"; exit 1; }
    
    # Step 2: Cluster clonal barcodes
    bartender_input="${sample_folder}/Bartender_input_address"
    if [[ -f "${bartender_input}" ]]; then
        while read -r line2; do
            new_name="${line2/.bartender/}"
            bartender_single_com -z -1 -d 1 -l 3 -f "${line2}" -o "${new_name}" || {
                echo "Error: Bartender clustering failed for ${line2}"
                exit 1
            }
        done < "${bartender_input}"
    else
        echo "Warning: Bartender input file not found for sample ${sampleID}. Skipping clustering..."
    fi

    # Combine data from one sample
    processed_sample_folder="${Step3_address}/${sampleID}"
    mkdir -p "${processed_sample_folder}" || { echo "Error: Failed to create directory ${processed_sample_folder}"; exit 1; }
    python3 "${Python_script_address}/cas12a_aggregate_barcode.py" \
        --a "${sample_folder}" --o "${processed_sample_folder}/" || {
        echo "Error: Step 2 failed for sample ${sampleID}"; exit 1;
    }

    # Clean up intermediate files to save disk space
    rm -rf "${clonal_barcode_folder}" || { echo "Error: Failed to clean up ${clonal_barcode_folder}"; exit 1; }
done < "${Input_data_info_address}"

# ---------------------------------------------------------------------
# Final Steps: Combine and Process All Data Across Samples
# ---------------------------------------------------------------------
# Step 3: Combine all data
python3 "${Python_script_address}/cas12a_aggregate_sample.py" --o "${Step3_address}/" || {
    echo "Error: Step 3 failed"; exit 1;
}

# Step 4: Combine all processed full data
python3 "${Python_script_address}/cas12a_aggregate_sample_for_QC.py" \
    --a "${Step1_address}/" --o "${Step3_address}/" || {
    echo "Error: Step 4 failed"; exit 1;
}

# ---------------------------------------------------------------------
# Job Summary: Display SLURM job statistics.
# ---------------------------------------------------------------------
# Output job statistics
sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 \
    --units=M -j "${SLURM_JOBID}"
