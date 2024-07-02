cd /share/pub/dengcy/Cancer_Gwas/src1.0/4.0.Run_scPagwas/

#'HeadandNeck' 'Melanoma' 
items=('EsophagealCancer' 'LungCancer' 'Pancreatic' 'BreastCancer' 'ThyroidCancer' 'BladderCancer' 'ColorectalCancer' 'LiverCancer' 'OvarianCancer' 'EndometrialCarcinoma' 'KidneyCancer' 'GastricCancer'  'ProstateCancer')
#items=()
nodes=("in009" "in007" "in008")
node_count=${#nodes[@]}
counter=0

for item in "${items[@]}"
do
echo $item
    node="${nodes[$counter % $node_count]}"
    filename="scdrs_${item}.sh"
    cat <<EOF > "$filename"
#!/bin/bash
#SBATCH -e scdrs_$item.err
#SBATCH -o scdrs_$item.out
#SBATCH -J scdrs_$item
#SBATCH -w $node
#SBATCH --mem=200000
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=dengcyelena@gmail.com
#SBATCH --time=1000:00:00
source activate mypy
cd /share/pub/dengcy/Cancer_Gwas/src1.0/4.0.Run_scPagwas/
python scDRS_pipeline.py --scRNA_h5ad_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/${item}_Immune_cell.h5ad --gene_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/${item}/merge_gene_PCC.csv --top_gene_num 500 --n_ctrl 200  --score_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/${item}/merge_scDRS_score.csv --group first_celltype_annotation --celltype_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/${item}/merge_celltypeP.csv
EOF
    counter=$((counter + 1))
    sbatch $filename
done

for i in {1561880..1561888}
do
    scancel $i
done

sbatch --begin=now+30minutes scdrs_LungCancer.sh




gwass=c("ieu-b-4810","ieu-b-4969","ieu-b-4954","bbj-a-133","ieu-b-4874","ukb-b-1316","ieu-b-4963","bbj-a-119","ieu-a-822","bbj-a-158","ieu-b-4965","bbj-a-107","finngen_r7_c3_oesophagus_exallc")#,
cancers=c("BreastCancer","Melanoma","LungCancer","LungCancer","BladderCancer","KidneyCancer","OvarianCancer","GastricCancer","GastricCancer","LiverCancer","ColorectalCancer","ColorectalCancer","EsophagealCancer")

source activate mypy
cd /share/pub/dengcy/Cancer_Gwas/src1.0/4.0.Run_scPagwas/
python scDRS_pipeline.py --scRNA_h5ad_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/LiverCancer_Immune_cell.h5ad --gene_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/LiverCancer/bbj-a-158_gene_PCC.csv --top_gene_num 500 --n_ctrl 200  --score_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/LiverCancer/bbj-a-158_scDRS_score.csv --group first_celltype_annotation --celltype_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/LiverCancer/bbj-a-158_celltypeP.csv
python scDRS_pipeline.py --scRNA_h5ad_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/ColorectalCancer_Immune_cell.h5ad --gene_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/ColorectalCancer/finngen_r7_c3_colorectal_gene_PCC.csv --top_gene_num 500 --n_ctrl 200  --score_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/ColorectalCancer/finngen_r7_c3_colorectal_scDRS_score.csv --group first_celltype_annotation --celltype_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/ColorectalCancer/finngen_r7_c3_colorectal_celltypeP.csv
cancers<-c("BladderCancer","GastricCancer","LiverCancer","ColorectalCancer")
gwass<-c('finngen_r7_c3_bladder_exallc','finngen_r7_c3_stomach_exallc','bbj-a-158','finngen_r7_c3_colorectal')


python scDRS_pipeline.py --scRNA_h5ad_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/GastricCancer_Immune_cell.h5ad --gene_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/GastricCancer/finngen_r7_c3_stomach_exallc_gene_PCC.csv --top_gene_num 500 --n_ctrl 200  --score_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/GastricCancer/finngen_r7_c3_stomach_exallc_scDRS_score.csv --group first_celltype_annotation --celltype_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/GastricCancer/finngen_r7_c3_stomach_exallc_celltypeP.csv

python /share/pub/dengcy/Cancer_Gwas/src1.0/4.0.Run_scPagwas/scDRS_pipeline.py --scRNA_h5ad_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/KidneyCancer_Immune_cell.h5ad --gene_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/KidneyCancer/finngen_r7_c3_kidney_notrenalpelvis_exallc_gene_PCC.csv --top_gene_num 500 --n_ctrl 200  --score_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/KidneyCancer/finngen_r7_c3_kidney_notrenalpelvis_exallc_scDRS_score.csv --group first_celltype_annotation --celltype_file /share/pub/dengcy/Cancer_Gwas/Runtime2.0/2.1scPagwas_run_mp/KidneyCancer/finngen_r7_c3_kidney_notrenalpelvis_exallc_celltypeP.csv
