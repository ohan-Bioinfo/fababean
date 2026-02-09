# Load required packages
if (!require("ape", quietly = TRUE)) {
    install.packages("ape", repos="http://cran.r-project.org", quiet=TRUE)
}
if (!require("phangorn", quietly = TRUE)) {
    install.packages("phangorn", repos="http://cran.r-project.org", quiet=TRUE)
}

library(ape)
library(phangorn)

cat("Reading PHYLIP format data...\n")

# Read the PHYLIP file we already created
phy_data <- read.dna("output/faba_150.phy", format="sequential")
cat(paste("Data dimensions:", dim(phy_data)[1], "samples,", dim(phy_data)[2], "SNPs\n"))

# Convert to phyDat format for phylogenetic analysis
phy_dat <- phyDat(phy_data, type="DNA")

# Calculate distance matrix
cat("Calculating distance matrix...\n")
dist_matrix <- dist.dna(phy_data, model="raw")

# Build Neighbor-Joining tree
cat("Building Neighbor-Joining tree...\n")
nj_tree <- nj(dist_matrix)

# Root the tree (using first sample as outgroup)
rooted_tree <- root(nj_tree, out=1)

# Save the tree
write.tree(rooted_tree, "output/faba_150_nj_tree.nwk")

# Plot the tree with better formatting
pdf("plots/faba_150_phylogenetic_tree.pdf", width=12, height=8)
plot(rooted_tree, main="Phylogenetic Tree - 150 SNPs\n(Neighbor-Joining)")
add.scale.bar()
dev.off()

# Also try with hclust
cat("Building tree with hierarchical clustering...\n")
hc <- hclust(dist_matrix, method="average")
hc_tree <- as.phylo(hc)

write.tree(hc_tree, "output/faba_150_hclust_tree.nwk")

pdf("plots/faba_150_hclust_tree.pdf", width=12, height=8)
plot(hc_tree, main="Phylogenetic Tree - 150 SNPs\n(Hierarchical Clustering)")
add.scale.bar()
dev.off()

cat("R analysis completed!\n")
cat("Trees saved to: output/faba_150_nj_tree.nwk and output/faba_150_hclust_tree.nwk\n")
