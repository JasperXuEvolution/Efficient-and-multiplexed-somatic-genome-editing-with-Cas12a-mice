#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=01-Bootstrapping-TK-KCa12a
#SBATCH --mail-user=xxx
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50g 
#SBATCH --time=10-24:00:00
#SBATCH --account=xxx
#SBATCH --partition=batch

# general input and output address
source ../config.sh
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate UltraSeq
input_data_address="${PROJECT_DIR}/02_data_cleaning_and_QC/data/Cas12a_TK_final_df.parquet"
working_dir="${PROJECT_DIR}/03_bootstrapping"

# command
python3 "$working_dir/main_code/UltraSeq_Bootstrapping_Cas12a.py" \
--a0 "$input_data_address" \
--a2 100 --a3 100 --a4 1000 --a5 KTHCas12a KTCas12a --a6 KT --a7 10 \
--o1 "$working_dir/data/TripleKnockout_BT" \
--o2 "$working_dir/data/TripleKnockout_BT" \
--l1 50 60 70 80 90 95 96 97 98 99 \
--m 'N' --c 'Yes' \
--gplist gene_combination gene_combination_unordered

sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID