####饼图
"/share/pub/dengcy/Cancer_Gwas/Runtime2.0/3.scPagwas_analysis/"
num_plots = len(cancers)
fig, axes = plt.subplots(1, num_plots, figsize=(60, 5), gridspec_kw={'width_ratios': [1] * num_plots})
df = pd.DataFrame(list(color_annotation.items()), columns=['Cell', 'Color'])
for i in range(0,15):
    print(i)
    cancer=cancers[i]
    p_df=pd.read_csv("proportion_"+cancer+"_celltype.csv")
    sorted_df = p_df.sort_values(by='propor', ascending=False).head(5)
    labels = sorted_df['celltype']
    sizes = sorted_df['propor']
    color_df = pd.DataFrame(list(color_annotation.items()), columns=['Cell', 'Color'])
    cell_color_dict = dict(zip(color_df['Cell'], color_df['Color']))
    color_vector = labels.map(cell_color_dict).tolist()
    axes[i].pie(sizes, labels=labels, autopct=None, startangle=140, colors=color_vector)
    axes[i].set_title(cancer)
    axes[i].axis('equal')
    axes[i].tick_params(labelsize=10)
plt.tight_layout()
plt.savefig('propor_celltype_pie_charts.pdf')
plt.close()