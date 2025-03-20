# Efficient and multiplexed somatic genome editing with Cas12a mice
Data analysis for "Efficient and multiplexed somatic genome editing with Cas12a mice"

## 1. Barcode extraction
* Refer to the folder S1_Extracting_barcodes
* The pipe line is in Ultra_seq_pipeline.sh
* Barcodes refer to both sgRNA sequence and clonal barcode
* Ultra_seq_pipeline.sh contains the necessary steps for barcode extraction
* The output of this steps is a master table contains the sgRNA, colonal barcode, sample ID and reads count information.

## 2. QC and preprocessing
* Refer to the folder S2_QC_and_Preprocessing
* Map clustered sgRNA to the reference.
* Filter out barcode that does not match the clonal barcode pattern
* QC of the data
* An example of reference experimental setup and sgRNA information is in folder Data/Chromatin_58Q
* For detailed steps of QC, please refer to the jupyter notebooks.

## 3. Bootstrapping and analysis
* Refer to the folder S3_Bootstrapping_and_Analysis
* We used a two-step boostrapping process to estimate tumor metrics.
* Example output is in folder Data/Chromatin_58Q
* 
* 

