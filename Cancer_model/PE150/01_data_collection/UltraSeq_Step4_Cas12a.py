#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import argparse
import glob
import os


def Combine_sgRNA_barcode_from_the_Same_mouse(input_folder_address):
    # find all file for each sgRNA
    temp_pattern = '/**/Clonal_barcode/*_cluster.csv' # When recursive is set, ** followed by a path separator matches 0 or more subdirectories.
    barcode_address_list = glob.glob(input_folder_address+temp_pattern, recursive=True)
    # address for sgRNA
    temp_ref_address = input_folder_address+'/Intermediate_df.csv'
    temp_bc_df_list = [] # store all df for each sgRNA into one list 
    for temp_a in barcode_address_list:
        t1 = temp_a.replace('cluster', 'barcode') # address for cluster.csv
        t2 = temp_a # address for barcode.csv
        t3 = temp_a.replace('_cluster.csv', '.bartender') # # address for bartender input
        temp_df = merge_barcode_and_sgRNA_output(t1,t2,t3)
        temp_bc_df_list.append(temp_df)
    temp_total_bc_df = pd.concat(temp_bc_df_list).reset_index(drop = True)
    temp_sgRNA_ref_df = pd.read_csv(temp_ref_address)
    temp_final = temp_total_bc_df.merge(temp_sgRNA_ref_df, on=['Read_ID','Clonal_barcode'])
    # deduplex
    temp_name_list = ['gRNA_combination', 'Clonal_barcode_center',
                  'gRNA1','gRNA3', 'Clonal_barcode', 'Sample_ID']
    temp_final_raw = temp_final.groupby(temp_name_list, as_index=False)['Read_ID'].count()
    temp_final_raw.rename(columns={'Read_ID':'Frequency'}, inplace=True)
    temp_final_true_sgRNA = temp_final # I only keep those perfect match, which I already did for cas12a
    temp_final_completely_deduplexed = temp_final_true_sgRNA.groupby(['gRNA_combination', 'Clonal_barcode_center','gRNA1','gRNA3', 'Sample_ID'], as_index=False)['Read_ID'].count()
    temp_final_completely_deduplexed.rename(columns={'Read_ID':'Frequency',
                                                    'Clonal_barcode_center':'Clonal_barcode'}, inplace=True)
    return (temp_final_raw, temp_final_completely_deduplexed)

def merge_barcode_and_sgRNA_output(input_barcode_address,input_cluster_address,bartender_input):
    temp_df1 = pd.read_csv(input_barcode_address).drop(columns = ['Frequency'])
    temp_df2 = pd.read_csv(input_cluster_address).drop(columns = ['Cluster.Score','time_point_1'])
    temp_df3 = pd.read_csv(bartender_input, sep =',', names=['Clonal_barcode','Read_ID'], header=None)
    temp_merge=pd.merge(temp_df1, temp_df2, how='inner', on=['Cluster.ID'],
         left_index=False, right_index=False, sort=True, copy=True, indicator=False,
         validate=None)
    temp_merge.rename(columns={'Unique.reads':'Clonal_barcode',
                              'Center':'Clonal_barcode_center'}, inplace=True)
    temp_merge = temp_merge.drop(columns = ['Cluster.ID']).merge(temp_df3, on='Clonal_barcode', how='right')
    return(temp_merge)
    
def main():
    parser = argparse.ArgumentParser(description='A function to combine sgRNA and clonal barcode information')
    parser.add_argument("--a", required=True, help="This is the input bartender folder")
    parser.add_argument("--o", required=True, help="This is the prefix of output file")
    args = parser.parse_args()
    temp_output_prefix = args.o
    MYDIR = args.a
    folder_paths = []
    for entry_name in os.listdir(MYDIR):
        entry_path = os.path.join(MYDIR, entry_name)
        if os.path.isdir(entry_path):
            folder_paths.append(entry_path)

    temp_final_df_list = []
    for temp_folder in folder_paths:
        temp_sampleID = temp_folder.split('/')[-1]
        temp_o1 = temp_output_prefix + temp_sampleID + '/Combined_ND_df.csv'
        temp_o2 = temp_output_prefix + temp_sampleID + '/Combined_deduplexed_df.csv'
        temp_df1, temp_df2 = Combine_sgRNA_barcode_from_the_Same_mouse(temp_folder)
        temp_df1.to_csv(temp_o1, index=False)
        temp_df2.to_csv(temp_o2, index=False)
        temp_final_df_list.append(temp_df2)
    temp_final_df = pd.concat(temp_final_df_list)
    temp_final_df.to_csv(temp_output_prefix + 'gRNA_clonalbarcode_combined.csv', index=False)
if __name__ == "__main__":
    main()  


