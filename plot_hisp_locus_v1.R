# Load required libraries
library(locuszoomr)
library(bigsnpr)
library(data.table)
library(EnsDb.Hsapiens.v75)

# Set working directory
setwd("/s/nath-lab/alsamman/____MyCodes____/locusZoom")

# 1. Load GWAS data
cat("Loading GWAS data...\n")
hisp_data <- fread("Hisp.txt")
head(hisp_data)

# Convert scientific notation if needed
hisp_data$p <- as.numeric(hisp_data$p)

# Add chromosome column with "chr" prefix if needed by locuszoomr
colnames(hisp_data) <- c("chrom", "pos", "p")
hisp_data<-hisp_data[1:100,]

# 2. Load reference panel using bigsnpr
cat("\nLoading reference panel...\n")

# If .rds doesn't exist, create it from bed/bim/fam
if (!file.exists("ref_panel/g1000_amr.rds")) {
  cat("Converting BED files to bigsnpr format...\n")
  snp_readBed("ref_panel/g1000_amr.bed")
}

# Now attach the reference panel
ref_panel <- snp_attach("ref_panel/g1000_amr.rds")

# Get genotype matrix and map
G <- ref_panel$genotypes
map <- ref_panel$map
cat("Reference panel loaded:", nrow(map), "SNPs\n")

# 3. Match SNPs between GWAS data and reference panel
cat("\nMatching SNPs...\n")
# Filter reference panel to chromosome 1
map_chr1 <- map[map$chromosome == 1, ]
G_chr1 <- G[, map$chromosome == 1]

# Find overlapping positions
hisp_data$in_ref <- hisp_data$pos %in% map_chr1$physical.pos

# Create a match between GWAS and reference panel positions
matched_indices <- match(hisp_data$pos, map_chr1$physical.pos)
hisp_data$ref_idx <- matched_indices

cat("SNPs in GWAS data:", nrow(hisp_data), "\n")
cat("SNPs in reference panel (chr1):", nrow(map_chr1), "\n")
cat("Overlapping SNPs:", sum(!is.na(matched_indices)), "\n")

# 4. Calculate LD for overlapping SNPs
cat("\nCalculating LD matrix...\n")
# Get indices of SNPs present in both datasets
valid_snps <- which(!is.na(hisp_data$ref_idx))

if (length(valid_snps) > 0) {
  # Get reference panel indices
  ref_indices <- hisp_data$ref_idx[valid_snps]
  
  # Find index SNP (lowest p-value)
  index_idx <- valid_snps[which.min(hisp_data$p[valid_snps])]
  index_ref_idx <- hisp_data$ref_idx[index_idx]
  
  cat("Index SNP position:", hisp_data$pos[index_idx], 
      "with p-value:", hisp_data$p[index_idx], "\n")
  
  # Calculate LD (r2) between index SNP and all other SNPs
  # Using correlation on genotype matrix
  index_geno <- G_chr1[, which(map_chr1$physical.pos == hisp_data$pos[index_idx])]
  
  # Initialize r2 column
  hisp_data$r2 <- NA
  
  # Calculate r2 for each overlapping SNP
  for (i in seq_along(valid_snps)) {
    snp_idx <- valid_snps[i]
    ref_col <- which(map_chr1$physical.pos == hisp_data$pos[snp_idx])
    
    if (length(ref_col) == 1) {
      snp_geno <- G_chr1[, ref_col]
      # Calculate correlation and square it for r2
      r_val <- cor(index_geno[], snp_geno[], use = "complete.obs")
      hisp_data$r2[snp_idx] <- r_val^2
    }
  }
  
  cat("LD calculated for", sum(!is.na(hisp_data$r2)), "SNPs\n")
  
  # Set index SNP r2 to 1
  hisp_data$r2[index_idx] <- 1.0
  
} else {
  cat("Warning: No overlapping SNPs found!\n")
  hisp_data$r2 <- 0
}

# 5. Create rsid column if not present (use chr:pos format)
hisp_data$rsid <- paste0("chr", hisp_data$chrom, ":", hisp_data$pos)

# 6. Create locus plot
cat("\nCreating locus plot...\n")

# Define region to plot (adjust as needed)
region_start <- min(hisp_data$pos)
region_end <- max(hisp_data$pos)
region_center <- (region_start + region_end) / 2

# Create locus object
loc <- locus(
  data = hisp_data,
  seqname = 1,
  xrange = c(region_start, region_end),
  LD = "r2",
  ens_db = "EnsDb.Hsapiens.v75"
)

# Summary of locus
summary(loc)

# Create the plot
pdf("hisp_locus_plot.pdf", width = 10, height = 6)
locus_plot(loc, 
           border = TRUE,
           labels = "index")
dev.off()

cat("\nPlot saved as: hisp_locus_plot.pdf\n")

# Also create a ggplot2 version
pdf("hisp_locus_plot_gg.pdf", width = 10, height = 6)
locus_ggplot(loc, labels = "index")
dev.off()

cat("ggplot2 version saved as: hisp_locus_plot_gg.pdf\n")

# Optional: Create interactive plotly version
# Uncomment to use:
# library(plotly)
# p <- locus_plotly(loc)
# htmlwidgets::saveWidget(p, "hisp_locus_plot.html")
# cat("Interactive plot saved as: hisp_locus_plot.html\n")

cat("\nDone!\n")

