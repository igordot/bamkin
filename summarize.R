#!/usr/bin/env Rscript


##
## Plot correlations based on a table of frequencies.
##
## usage: Rscript --vanilla summarize.R snp.freq.csv
##


# increase output width
options(width = 150)

# command line arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 1) stop("usage: Rscript --vanilla summarize.R snp.freq.csv")
freq_file = args[1]

# check that input files exist
if (!file.exists(freq_file)) stop("file does not exist: ", freq_file)

# try to install from CRAN if package is not already installed and load it
package_list = c("RColorBrewer", "matrixStats", "corrplot", "ggplot2", "ggdendro", "ggrepel", "cowplot")
for (p in package_list) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) {
    message("installing package ", p)
    install.packages(p, quiet = TRUE, repos = "http://cran.rstudio.com")
  }
  library(p, character.only = TRUE, quietly = TRUE)
}

message("importing ", freq_file)
freq_table = read.csv(freq_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
freq_table = as.matrix(freq_table)
message("cols/samples: ", ncol(freq_table))
message("rows/snps input: ", nrow(freq_table))

# filter for informative SNPs (varying frequencies across samples)
freq_table = freq_table[rowSds(freq_table) > 0.1,]
message("rows/snps variable: ", nrow(freq_table))

# subset to most variable SNPs to keep table size reasonable
max_snps = 10000
if (nrow(freq_table) > max_snps) {
  var_snps = order(rowVars(freq_table), decreasing = TRUE)[1:max_snps]
  freq_table = freq_table[var_snps,]
}

message("correlation")
freq_cor = cor(freq_table, method = "pearson")

# correlation plot (clustered using hclust)
message("correlation plot")
corrplot_colors = colorRampPalette(c("#FFFFFF", "#FAFAFB", "#F6F6F7", "#F1F2F3", "#CFD2D6", "#C0C5CA", "#031B5B"))(40)
png("plot.corr.png", width = 12, height = 12, units = "in", res = 200, pointsize = 10, family = "Arial Black")
corrplot(freq_cor, method = "color", type = "upper", order = "hclust", tl.col = "black", cl.ratio = 0.1, col = corrplot_colors)
dev.off()

# hiearchical cluster analysis using the euclidan distance between variables based on the correlation
message("hierarchical cluster plot")
freq_hclust = hclust(dist(freq_cor))
hclust_plot = ggdendrogram(freq_hclust, rotate = TRUE) +
theme(axis.text.x = element_blank())
save_plot(filename = "plot.hclust.png", plot = hclust_plot, base_width = 3, base_height = 5, units = "in", dpi = 300)

# PCA
message("PCA")
pca = prcomp(t(freq_table))
pca_data = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
pca_plot = ggplot(data = pca_data, aes_string(x = "PC1", y = "PC2")) +
geom_point(size = 2) +
theme(axis.text = element_blank(), axis.ticks = element_blank()) +
geom_text_repel(aes(label = rownames(pca_data)))
save_plot(filename = "plot.pca.png", plot = pca_plot, base_width = 5, base_height = 5, units = "in", dpi = 300)


# end
