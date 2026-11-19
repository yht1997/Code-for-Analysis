# Clear workspace
rm(list = ls()) 

# Load required packages (auto-install missing ones)
required_pkgs <- c("clusterProfiler", "enrichplot", "ggplot2", "GOplot", "BiocManager")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  if ("BiocManager" %in% missing_pkgs) install.packages("BiocManager")
  BiocManager::install(setdiff(missing_pkgs, "BiocManager"), update = FALSE)
}

# Load species annotation databases (batch check + install + load)
target_dbs <- c("org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db")  # human/mouse/rat databases
for (db in target_dbs) {
  if (!require(db, character.only = TRUE, quietly = TRUE)) {
    cat("Installing database: ", db, "\n")
    BiocManager::install(db, update = FALSE, quiet = TRUE)
  }
  library(db, character.only = TRUE, quietly = TRUE)
  cat("Loaded database: ", db, "\n")
}

# Load core analysis packages
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(GOplot)

# --------------------------
# Core configuration: species selection and parameter mapping
# --------------------------
# 1. Define database and KEGG code for each species (key: add kegg_code)
species_list <- list(
  human = list(org_db = "org.Hs.eg.db", kegg_code = "hsa"),  # human KEGG code: hsa
  mouse = list(org_db = "org.Mm.eg.db", kegg_code = "mmu"),  # mouse KEGG code: mmu
  rat   = list(org_db = "org.Rn.eg.db", kegg_code = "rno")   # rat KEGG code: rno
)

# 2. Select target species (only modify here to switch)
species <- "mouse"  # options: human/mouse/rat

# 3. Validate species (keep only one necessary check)
if (!species %in% names(species_list)) {
  stop("Invalid species parameter! Options: ", paste(names(species_list), collapse = "/"))
}

# 4. Extract configuration for current species (fix parameter extraction)
species_info <- species_list[[species]]
org_db_name <- species_info$org_db  # extract database name
kegg_org <- species_info$kegg_code  # extract KEGG code

# Load database object for current species
org_db <- get(org_db_name)
cat("Current species: ", species, ", database: ", org_db_name, ", KEGG code: ", kegg_org, "\n")

# --------------------------
# Data preparation and cleaning
# --------------------------
# Set working directory and create result folders
setwd("D:/Network_Toxicology/GO_KEGG_analysis")
dir.create("GO_analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("KEGG_analysis", showWarnings = FALSE, recursive = TRUE)

# Core: only modify the file name here to switch input file
input_file <- "toxicant_disease_intersect_genes.csv"  # input file name (with extension)

# Read and clean gene data (simplified validation logic)
if (!file.exists(input_file)) {
  stop("Input file does not exist! Please check whether '", input_file, "' is in the working directory")
}
genes_df <- tryCatch(
  read.csv(input_file, stringsAsFactors = FALSE, fileEncoding = "UTF-8"),
  error = function(e) stop("File read failed: ", e$message)
)
if (!"gene" %in% colnames(genes_df)) {
  stop("Input file must contain a column named 'gene' (gene symbols)")
}

# Gene symbol cleaning
genes_clean <- unique(na.omit(genes_df$gene))
genes_clean <- gsub("^\\s+|\\s+$", "", genes_clean)  # remove leading/trailing spaces
genes_clean <- genes_clean[genes_clean != ""]
cat("Number of valid genes after cleaning: ", length(genes_clean), "\n")
if (length(genes_clean) == 0) stop("No valid genes! Please check the input file")

# --------------------------
# Gene ID conversion (SYMBOL -> ENTREZID)
# --------------------------
entrezIDs <- bitr(
  gene = genes_clean,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org_db
)
if (nrow(entrezIDs) == 0) {
  stop("ID conversion failed! Possible reasons: 1. gene symbols do not match the species; 2. spelling errors in symbols")
}
cat("Number of successfully converted ENTREZIDs: ", nrow(entrezIDs), " (", 
    round(nrow(entrezIDs)/length(genes_clean)*100, 1), "% of cleaned genes)\n", sep = "")
gene_entrez <- entrezIDs$ENTREZID

# --------------------------
# GO enrichment analysis
# --------------------------
go_enrich <- enrichGO(
  gene = gene_entrez,
  OrgDb = org_db,
  ont = "all",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

if (nrow(go_enrich) == 0) {
  warning("No significantly enriched GO terms! Suggest relaxing pvalueCutoff to 0.1 and retrying")
} else {
  write.table(go_enrich, "GO_analysis/GO_enrichment_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Bar plot
  pdf("GO_analysis/GO_barplot.pdf", width = 10, height = 8)
  print(barplot(go_enrich, drop = TRUE, showCategory = 6, label_format = 50, 
                split = "ONTOLOGY", color = "p.adjust") +
          facet_grid(ONTOLOGY ~ ., scale = "free") +
          theme(plot.title = element_text(hjust = 0.5, size = 14)) +
          ggtitle("GO Analysis (BP/CC/MF)"))
  dev.off()
  
  # Dot plot
  pdf("GO_analysis/GO_dotplot.pdf", width = 10, height = 8)
  print(dotplot(go_enrich, showCategory = 6, label_format = 50, 
                split = "ONTOLOGY", color = "p.adjust") +
          facet_grid(ONTOLOGY ~ ., scale = "free") +
          theme(plot.title = element_text(hjust = 0.5, size = 14)) +
          ggtitle("GO Analysis"))
  dev.off()
}

# --------------------------
# KEGG pathway enrichment analysis
# --------------------------
kegg_enrich <- enrichKEGG(
  gene = gene_entrez,
  organism = kegg_org,
  keyType = "kegg",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  pAdjustMethod = "fdr"
)

if (nrow(kegg_enrich) == 0) {
  warning("No significantly enriched KEGG pathways! Suggestions: 1. check whether the species is correct; 2. relax pvalueCutoff")
} else {
  write.table(kegg_enrich, "KEGG_analysis/KEGG_enrichment_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Bar plot
  pdf("KEGG_analysis/KEGG_barplot.pdf", width = 10, height = 7)
  print(barplot(kegg_enrich, drop = TRUE, showCategory = 8, label_format = 50, color = "p.adjust") +
          theme(plot.title = element_text(hjust = 0.5, size = 14)) +
          ggtitle("KEGG Analysis"))
  dev.off()
  
  # Dot plot
  pdf("KEGG_analysis/KEGG_dotplot.pdf", width = 10, height = 7)
  print(dotplot(kegg_enrich, showCategory = 8, label_format = 50, color = "p.adjust") +
          theme(plot.title = element_text(hjust = 0.5, size = 14)) +
          ggtitle("KEGG Analysis"))
  dev.off()
}

cat("\nAnalysis complete! Results have been saved to the corresponding folders by category.\n",
    "If no significant results were obtained, please adjust parameters according to the warning messages.\n")
