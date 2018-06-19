#!/usr/bin/env Rscript


##
## Plot correlations based on a table of frequencies.
##
## usage: Rscript --vanilla summarize.R snp.freq.csv
##


# increase output width
options(width = 120)

# command line arguments
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript --vanilla summarize.R <snp.freq.csv>")
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

# limits
min_snps = 10
max_snps = 10000

message("importing ", freq_file)
freq_table = read.csv(freq_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
freq_table = as.matrix(freq_table)

# stop if too few SNPs in input file
if (nrow(freq_table) < min_snps) stop("too few SNPs across all samples: ", nrow(freq_table))

num_input_samples = ncol(freq_table)
message("cols/samples: ", num_input_samples)
num_input_snps = nrow(freq_table)
message("rows/SNPs input: ", num_input_snps)

# stop if too few variable SNPs
num_var_snps = count(rowSds(freq_table) > 0.1)
if (num_var_snps < min_snps) stop("too few variable SNPs: ", num_var_snps)

# keep only informative SNPs with varying frequencies across samples
freq_table = freq_table[rowSds(freq_table) > 0.1, ]

# stop if too few highly variable (present at both high and low frequencies) SNPs
num_high_var_snps = count((rowMaxs(freq_table) - rowMins(freq_table)) > 0.7)
if (num_high_var_snps < min_snps) stop("too few ref/alt SNPs:", num_high_var_snps)

message("rows/SNPs variable: ", nrow(freq_table))

# subset to just most variable SNPs if table is too large
if (nrow(freq_table) > max_snps) {
  message("subsetting table to most variable SNPs")
  var_order_snps = order(rowVars(freq_table), decreasing = TRUE)[1:max_snps]
  freq_table = freq_table[var_order_snps, ]
}

# save filtered table
write.csv(cbind("#POS" = rownames(freq_table), freq_table), file = "snp.freq.filtered.csv", row.names = FALSE, quote = FALSE)

# set ggplot geom_text font size based on number of samples (measured in mm, not points)
if (num_input_samples > 50) {
  font_size = 1.5
} else if (num_input_samples > 20) {
  font_size = 2
} else if (num_input_samples > 10) {
  font_size = 3
} else {
  font_size = 4
}

# correlation
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
freq_dendro = dendro_data(freq_hclust, type = "rectangle")
hclust_plot =
  ggplot() +
  geom_segment(data = segment(freq_dendro), aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = label(freq_dendro), aes(x = x, y = y, label = label, hjust = 0), size = font_size) +
  coord_flip() +
  scale_y_reverse(expand = c(0.8, 0)) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), text = element_blank(), line = element_blank())
save_plot(filename = "plot.hclust.png", plot = hclust_plot, base_width = 3, base_height = 5, units = "in", dpi = 300)

# PCA
message("PCA")
pca = prcomp(t(freq_table))
pca_data = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
pca_plot =
  ggplot(data = pca_data, aes_string(x = "PC1", y = "PC2")) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = rownames(pca_data)),
                  size = font_size, point.padding = unit(0.5, "lines"), color = "black") +
  theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())
save_plot(filename = "plot.pca.png", plot = pca_plot, base_width = 5, base_height = 5, units = "in", dpi = 300)

# delete Rplots.pdf
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")



# end
