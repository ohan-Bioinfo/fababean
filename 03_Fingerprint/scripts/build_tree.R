# Load required packages
if (!require("ape", quietly = TRUE)) {
    install.packages("ape", repos="http://cran.r-project.org", quiet=TRUE)
}
if (!require("phangorn", quietly = TRUE)) {
    install.packages("phangorn", repos="http://cran.r-project.org", quiet=TRUE)
}

library(ape)
library(phangorn)

cat("Reading genotype data...\n")

# Read the genotype TSV file
genotype_data <- read.table("output/faba_150_gt.tsv", header=TRUE, row.names=1, sep="\t", na.strings="NA")

# Convert genotypes to numeric (simple approach)
convert_genotypes <- function(gt_matrix) {
    num_matrix <- matrix(NA, nrow=nrow(gt_matrix), ncol=ncol(gt_matrix))
    
    for(i in 1:nrow(gt_matrix)) {
        for(j in 1:ncol(gt_matrix)) {
            gt <- as.character(gt_matrix[i,j])
            if(is.na(gt) || gt == "NA") {
                num_matrix[i,j] <- NA
            } else if(gt == "A/A" || gt == "T/T" || gt == "C/C" || gt == "G/G") {
                num_matrix[i,j] <- 0  # Homozygous
            } else if(gt == "A/T" || gt == "T/A" || gt == "A/C" || gt == "C/A" || 
                      gt == "A/G" || gt == "G/A" || gt == "T/C" || gt == "C/T" ||
                      gt == "T/G" || gt == "G/T" || gt == "C/G" || gt == "G/C") {
                num_matrix[i,j] <- 1  # Heterozygous
            } else {
                num_matrix[i,j] <- NA
            }
        }
    }
    return(num_matrix)
}

cat("Converting genotypes to numeric format...\n")
numeric_genotypes <- convert_genotypes(genotype_data)
rownames(numeric_genotypes) <- rownames(genotype_data)

# Remove samples with too much missing data
missing_per_sample <- apply(numeric_genotypes, 1, function(x) sum(is.na(x)))
cat("Missing data per sample:\n")
print(missing_per_sample)

# Use complete cases only for distance calculation
complete_cases <- numeric_genotypes[missing_per_sample == 0, ]
cat(paste("Using", nrow(complete_cases), "samples with complete data\n"))

# Calculate distance matrix
cat("Calculating distance matrix...\n")
dist_matrix <- dist(complete_cases, method="manhattan")

# Build Neighbor-Joining tree
cat("Building Neighbor-Joining tree...\n")
nj_tree <- nj(dist_matrix)

# Root the tree (using first sample as outgroup)
rooted_tree <- root(nj_tree, out=1)

# Save the tree
write.tree(rooted_tree, "output/faba_150_nj_tree.nwk")

# Plot the tree
pdf("plots/faba_150_phylogenetic_tree.pdf", width=10, height=8)
plot(rooted_tree, main="Phylogenetic Tree - 150 SNPs\n(Neighbor-Joining)")
add.scale.bar()
dev.off()

# Also try with hclust
cat("Building tree with hierarchical clustering...\n")
hc <- hclust(dist_matrix, method="average")
hc_tree <- as.phylo(hc)

write.tree(hc_tree, "output/faba_150_hclust_tree.nwk")

pdf("plots/faba_150_hclust_tree.pdf", width=10, height=8)
plot(hc_tree, main="Phylogenetic Tree - 150 SNPs\n(Hierarchical Clustering)")
add.scale.bar()
dev.off()

cat("R analysis completed!\n")
cat("Trees saved to: output/faba_150_nj_tree.nwk and output/faba_150_hclust_tree.nwk\n")
