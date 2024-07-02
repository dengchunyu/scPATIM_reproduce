from scipy.stats import norm

BaseDirectory = '/share/pub/dengcy/Cancer_Gwas/Runtime2.0/1.cancer_scdata/'
#
cancers=['HeadandNeck','BreastCancer','ProstateCancer','Melanoma','LungCancer','BladderCancer','KidneyCancer','ThyroidCancer','OvarianCancer','Pancreatic','GastricCancer','LiverCancer','ColorectalCancer','EsophagealCancer','EndometrialCarcinoma']
cl_annotation = {'DC':'DC', 'Bcells':'B/Plasma' , 'Treg':'Treg', 'Bgc':'B/Plasma', 'Tem/Temra.cyt':'T.cyt', 'NK.CD16neg':'NK','Bpgc':'B/Plasma', 'Tcm/N.cyt':'T.cyt', 'DC2':'DC', 'pDC':'DC', 'Mac':'Mac/Mono', 'Mono.nc':'Mac/Mono','Mono.c':'Mac/Mono', 'Tcm/N.h':'Th', 'Tem/Trm.cyt':'T.cyt', 'NK':'NK', 'Trm.cyt':'T.cyt', 'Th17':'Th', 'Tem/Eff.h':'Th', 'NK.CD16pos':'NK', 'DC.mig':'DC', 'Bnaive':'B/Plasma', 'Mast':'Mast', 'Macrophages':'Mac/Mono', 'ILC':'ILC', 'Plasma':'B/Plasma', 'Bmem':'B/Plasma', 'Mac.inter':'Mac/Mono', 'Tgd':'Tgd','MAIT': 'MAIT','Ery.Mac':'Mac/Mono','HSC/MPP':'HSC/MPP'}
os.chdir("/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/")

for i in cancers:
    print(i)
    results_file=os.path.join(BaseDirectory,i+'_Immune_cell.h5ad')
    adata=sc.read_h5ad(results_file)
    adata.obs["merge_celltype_annotation"] = adata.obs.first_celltype_annotation.map(cl_annotation)
    adata.obs['merge_celltype_annotation'] = adata.obs['merge_celltype_annotation'].astype('category')
    adata.write(results_file)