#下载参考数据集
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz


tar -xzvf cellranger-7.0.1.tar.gz
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz

cd /share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653
wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/scrBT1430m_S0_L002_R1_001.fastq.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/scrBT1430m_S0_L002_R2_001.fastq.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/scrBT1431m_S0_L003_R1_001.fastq.gz

wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/scrBT1431m_S0_L003_R2_001.fastq.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/scrBT1432m_S0_L004_R1_001.fastq.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/scrBT1432m_S0_L004_R2_001.fastq.gz

wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/BT1375_S5_L004_R1_001.fastq.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/BT1375_S5_L004_R2_001.fastq.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/BT1376_S6_L005_R1_001.fastq.gz

wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/BT1376_S6_L005_R2_001.fastq.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/BT1377_S7_L006_R1_001.fastq.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/BT1377_S7_L006_R2_001.fastq.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/scrBT1425_hg19_S11_L006_R1_001.fastq.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/scrBT1425_hg19_S11_L006_R2_001.fastq.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/scrBT1426_hg19_S12_L006_R1_001.fastq.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/scrBT1426_hg19_S12_L006_R2_001.fastq.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/scrBT1427_hg19_S13_L007_R1_001.fastq.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6653/scrBT1427_hg19_S13_L007_R2_001.fastq.gz

## cellranger count命令
## --id给你这次的运行七个名字，如sample345
## --fastqs 输入分析数据所在路径
## --transcriptome 输入参考基因组所在路径
cd /share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653
/share/pub/dengcy/software/cellranger-7.0.1/cellranger count --id=BT1376_S6_L005 \
--sample BT1376 \
--fastqs=/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653 \
--transcriptome=/share/pub/dengcy/software/refdata-gex-GRCh38-2020-A \
--nosecondary

cd /share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653
/share/pub/dengcy/software/cellranger-7.0.1/cellranger count --id=BT1377_S7_L006 \
--sample BT1377 \
--fastqs=/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653 \
--transcriptome=/share/pub/dengcy/software/refdata-gex-GRCh38-2020-A \
--nosecondary

cd /share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653
/share/pub/dengcy/software/cellranger-7.0.1/cellranger count --id=BT1375_S5_L004 \
--sample BT1375 \
--fastqs=/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653 \
--transcriptome=/share/pub/dengcy/software/refdata-gex-GRCh38-2020-A \
--nosecondary

cd /share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653
/share/pub/dengcy/software/cellranger-7.0.1/cellranger count --id=scrBT1426_hg19_S12_L006 \
--sample scrBT1426_hg19 \
--fastqs=/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653 \
--transcriptome=/share/pub/dengcy/software/refdata-gex-GRCh38-2020-A \
--nosecondary

cd /share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653
/share/pub/dengcy/software/cellranger-7.0.1/cellranger count --id=scrBT1432m_S0_L004 \
--sample scrBT1432m \
--fastqs=/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653 \
--transcriptome=/share/pub/dengcy/software/refdata-gex-GRCh38-2020-A \
--nosecondary

############################
cd /share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653
/share/pub/dengcy/software/cellranger-7.0.1/cellranger count --id=scrBT1427_hg19_S13_L007 \
--sample scrBT1427_hg19 \
--fastqs=/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653 \
--transcriptome=/share/pub/dengcy/software/refdata-gex-GRCh38-2020-A \
--nosecondary

cd /share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653
/share/pub/dengcy/software/cellranger-7.0.1/cellranger count --id=scrBT1431m_S0_L003 \
--sample scrBT1431m \
--fastqs=/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653 \
--transcriptome=/share/pub/dengcy/software/refdata-gex-GRCh38-2020-A \
--nosecondary

cd /share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653
/share/pub/dengcy/software/cellranger-7.0.1/cellranger count --id=scrBT1430m_S0_L002 \
--sample scrBT1430m \
--fastqs=/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653 \
--transcriptome=/share/pub/dengcy/software/refdata-gex-GRCh38-2020-A \
--nosecondary

cd /share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653
/share/pub/dengcy/software/cellranger-7.0.1/cellranger count --id=scrBT1425_hg19_S11_L006 \
--sample scrBT1425_hg19 \
--fastqs=/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/LungCancer/E-MTAB-6653 \
--transcriptome=/share/pub/dengcy/software/refdata-gex-GRCh38-2020-A \
--nosecondary
