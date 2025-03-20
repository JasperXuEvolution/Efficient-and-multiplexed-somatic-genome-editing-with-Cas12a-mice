#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=xxx
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
# NGS_address is a file contains the directory of raw reads
Input_data_info_address="$Project_directory/NGS_address"
Step1_address="$Project_directory/Bartender"
Step3_address="$Project_directory/Processed_data"

# comand
mkdir -p $Step1_address
mkdir -p $Step3_address
while read -r line;
do
   r1=$(echo "$line" | cut -d',' -f1);
   r2=$(echo "$line" | cut -d',' -f2);
   sampleID=$(echo "$line" | cut -d',' -f3);
   
   # Step1
   temp_folder1=$Step1_address/$sampleID
   mkdir -p $temp_folder1
   temp_folder2=$Step1_address/$sampleID/Clonal_barcode
   mkdir -p $temp_folder2

   # Generate bartender input
   python3 $Python_script_address/UltraSeq_Step1_Cas12a.py --a1 $r1 --a2 $r2 --b "$Project_directory/sgRNA_Cas12a_combine.csv"\
   --o "$temp_folder1"

   
   # Step2 
   # This step deals with clustering of clonal barcode
   while read -r line2;
   do 
      new_name=${line2/.bartender/}
      # echo $line2
      # echo $new_name
      bartender_single_com -z -1 -d 1 -l 3 -f "$line2" -o "$new_name"
   done < $Step1_address/$sampleID/Bartender_input_address

   # Step3
   temp_folder3=$Step3_address/$sampleID
   mkdir -p $temp_folder3
done < $Input_data_info_address

# Step4
# Combined all the data
python3 $Python_script_address/UltraSeq_Step4_Cas12a.py --a "$Step1_address" \
--o "$Step3_address/"

sacct --format=JobID,JobName,Submit,Start,End,State,Partition,ReqTRES%30,CPUTime,MaxRSS,NodeList%30 --units=M -j $SLURM_JOBID