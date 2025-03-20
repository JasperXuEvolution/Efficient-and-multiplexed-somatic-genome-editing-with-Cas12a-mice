#!/usr/bin/env python
# coding: utf-8
#function
import gzip
import regex
import argparse
import sys
import pandas as pd

def find_reverse_complementary(input_string):
    temp_dic = {'A':'T','G':'C','T':'A','C':'G','N':'N',
                'a':'t','g':'c','t':'a','g':'c','n':'n'}
    return(''.join(tuple([temp_dic.get(x) for x in input_string[::-1]])))
SpikeInList = ['GAGTTCTGCCTCAAGCAAGT','GTTGAATCGTCCGTACATGT','TATACTTGCACCATGCCATA']

def count_spikein(fastqgz_input_address1,fastqgz_input_address2,spikein_list):
    # given data adress and spikein sequence, this function will count the reads number for each spikein
    # fastqgz_input_address1 is the address for the read1 
    # fastqgz_input_address2 is the address for the read1 
    temp_spikein_dic = {}
    temp_pattern = regex.compile('(TAGTT){e<2}' + '(.{16})' + 'TATGG'+'(.{20})' + 'GTT(TAAGA){e<2}')
    with gzip.open(fastqgz_input_address1,'rt') as handler1, gzip.open(fastqgz_input_address2,'rt') as handler2:
        temp_readID = handler1.readline().rstrip() # read ID
        handler2.readline().rstrip()
        temp_sequence1 = handler1.readline().rstrip()
        temp_sequence2 = find_reverse_complementary(handler2.readline().rstrip())
        handler1.readline() # skip two lines
        handler1.readline()
        handler2.readline() # skip two lines
        handler2.readline()
        while temp_readID:
            temp_search_result1 = temp_pattern.search(temp_sequence1)
            if temp_search_result1:
                temp_search_result2 = temp_pattern.search(temp_sequence2)
                if temp_search_result2:
                    temp_gRNA_r1 = temp_search_result1.group(3) # gRNA based on read1
                    temp_clonal_barcode_r1 = temp_search_result1.group(2) # barcode based on read1
                    temp_gRNA_r2= temp_search_result2.group(3) # gRNA based on read2
                    temp_clonal_barcode_r2 = temp_search_result2.group(2) # barcode based on read2
                    if (temp_gRNA_r1==temp_gRNA_r2)&(temp_clonal_barcode_r1==temp_clonal_barcode_r2): # both barcode and gRNA are identical:
                        if temp_gRNA_r1 in spikein_list:
                            temp_key = (temp_gRNA_r1,temp_clonal_barcode_r1)
                            if temp_key in temp_spikein_dic.keys():
                                temp_spikein_dic[temp_key]+=1
                            else:
                                temp_spikein_dic[temp_key]=1
            temp_readID = handler1.readline().rstrip() # read ID
            handler2.readline().rstrip() # read ID
            temp_sequence1 = handler1.readline().rstrip()
            temp_sequence2 = find_reverse_complementary(handler2.readline().rstrip())
            handler1.readline() # skip two lines
            handler1.readline()
            handler2.readline() # skip two lines
            handler2.readline()
    temp_gRNA_list = []
    temp_barcode_list = []
    temp_count_list = []
    for key,value in temp_spikein_dic.items():
        temp_gRNA_list.append(key[0])
        temp_barcode_list.append(key[1])
        temp_count_list.append(value)
    # generate output dataframe
    temp_df = pd.DataFrame({'gRNA':temp_gRNA_list, 'Clonal_barcode':temp_barcode_list,'Count':temp_count_list})
    return(temp_df)

def main():
    parser = argparse.ArgumentParser(description='A function to count the number of spike-in the data')
    parser.add_argument("--a", required=True, help="This is the adress of NGS adress file")
    parser.add_argument("--o", required=True, help="This is the dir of output file")
    args = parser.parse_args()
    NGS_address = args.a # data address
    output_dir = args.o # output dir
    temp_df_list = []
    with open(NGS_address, 'rt') as handler:
        line = handler.readline().rstrip()
        while line:
            r1_address,r2_address, Sample_ID = line.split(',')
            temp_df = count_spikein(r1_address,r2_address,SpikeInList)
            temp_df['Sample_ID'] = Sample_ID
            temp_df_list.append(temp_df)
            line = handler.readline().rstrip()

    Final_df = pd.concat(temp_df_list, ignore_index=True) # combine all the df from each sample
    temp_address_file = output_dir+'/SpikeInRead.csv' # output address
    Final_df.to_csv(temp_address_file, index=False)

if __name__ == "__main__":
    main()  


