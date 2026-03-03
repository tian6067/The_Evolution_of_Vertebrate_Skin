#!/bin/bash
#SBATCH --job-name=scFEA
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --partition=Fnode1
#SBATCH --mem 800G
#SBATCH --mail-user=837588643@qq.com
#SBATCH --mail-type=ALL
#SBATCH -o %A

cd ~/miniconda3/bin/
source activate
conda activate scFEA
cd /public/home/b20213040320/scFEA/scFEA
python /public/home/b20213040320/scFEA/scFEA/src/scFEA.py \
    --data_dir /public/home/b20213040320/scFEA/scFEA/data \
    --input_dir /public/home/b20213040320/AAA_skin/kc/scFEA \
    --moduleGene_file module_gene_m168.csv \
    --test_file count_matrix.csv \
    --cName_file cName_c70_m168.csv \
    --sc_imputation True \
    --stoichiometry_matrix cmMat_c70_m168.csv \
    --output_flux_file output/adj_flux.csv \
    --output_balance_file output/adj_balance.csv
