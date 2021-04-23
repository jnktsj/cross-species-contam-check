library(pheatmap)

# create species labels
species = read.table("species_list.txt", row.names=2, as.is=TRUE)
labels = paste0(species$V1, " (", rownames(species),")")
names(labels) = rownames(species)

# create heatmap
input_data = "chr21_alt_blat_scores_v5.txt"
#input_data = "chr21_alt_blat_scores.txt"
data = read.table(input_data, header=TRUE, row.names=1, as.is=TRUE)
colnames(data) = labels[colnames(data)]
#pheatmap(t(data), show_colnames=FALSE, show_rownames=TRUE, cluster_cols=FALSE, cluster_rows=TRUE)
# plot without clustering species species 
pheatmap(t(data), show_colnames=FALSE, show_rownames=TRUE, cluster_cols=FALSE, cluster_rows=FALSE)

# plot VAF distribution between v3 and v5 tumors
data = read.table("chr21_VAF_v3_and_v5.txt", header=TRUE)
library(ggpubr)
ggscatter(data, x="vaf_v5", y="vaf_v3", add ="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="kendall", xlab="v5", ylab="v3", shape=20, color=rgb(t(col2rgb("steelblue")/255), alpha=0.7), title="Empirical VAF")