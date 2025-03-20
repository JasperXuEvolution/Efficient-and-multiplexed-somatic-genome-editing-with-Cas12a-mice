#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=xxxx
#SBATCH --mail-user=xxx
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=50g 
#SBATCH --time=24:00:00
#SBATCH --account=xxx
#SBATCH --partition=batch



# Overall directory
LP="Cancer_model/"
# this project name and address 
Input_experiment_ID="PE150"
Project_directory=$LP$Input_experiment_ID
# this is the address for python scripts used  
Python_script_address=$LP$Input_experiment_ID

# modules
module load adapterremoval/2.3.1
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate UltraSeq

# input and output address
Input_data_info_address="$Project_directory/NGS_address"
# Generate Spikein read dataframe
python3 $Python_script_address/UltraSeq_SpikeInCount_Cas12a.py --a $Input_data_info_address --o $Project_directory
sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID