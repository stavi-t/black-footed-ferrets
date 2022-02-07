#!/bin/sh
#SBATCH -N 1                 # nodes=4
#SBATCH --ntasks-per-node=4  # ppn=20
#SBATCH -J MACAU1              
#SBATCH -t 144:00:00             
#SBATCH --mail-user=st5978@princeton.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mem=50GB    

# 11.19.2021
# dataset = blood8; covariates = intercept, age; predictor = avg litter size
dirr_out = /tigress/VONHOLDT/ALD/BFF_RRBS/9_BFF_MACAU_PED_df_v2

cd /tigress/VONHOLDT/ALD/BFF_RRBS/9_BFF_MACAU_PED_df_v2
/home/vonholdt/VONHOLDT/BIN/macau-1.00/bin/macau -g BFF_CG_Meth_10x_merged_df_BLOOD_modified_methylCounts_header.txt -t BFF_CG_Meth_10x_merged_df_BLOOD_modified_totalCounts_header.txt -p predictor_avglittersize_blood8.txt -k relatedness_pedigree_blood8.txt -c covariates_blood8_df_omitPC.txt -bmm -o MACAU_output_avglittersize_blood8_df -outdir dirr_out


