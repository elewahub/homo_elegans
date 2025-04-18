wbpheno <- read.csv('/Users/ahmedelewa/Library/CloudStorage/GoogleDrive-pi@elewalabs.com/My Drive/01_projects/074_augsburg/WBPhenotypes.csv',header=TRUE,row.names = 1)
wbpheno <- wbpheno[wbpheno$RNAi.Phenotype.Observed!='N.A.'|wbpheno$Allele.Phenotype.Observed!='N.A.',]
dim(wbpheno)
#[1] 10411     2

# Convert the data frame into the desired format
wpo <- do.call(rbind, lapply(rownames(wbpheno), function(id) {
  # Extract features from columns A and B for each row
  features <- unlist(strsplit(c(wbpheno[id, "RNAi.Phenotype.Observed"], wbpheno[id, "Allele.Phenotype.Observed"]), ",\\s*"))  # Split by comma and trim spaces
  features <- features[features != "N.A."]  # Remove N.A.
  # Create a new data frame for this ID
  data.frame(ID = id, Feature = features)
}))

colnames(wpo) <- c('ID','phenotype')
wpo <- unique(wpo)

ortho <- read.csv('/Users/ahmedelewa/Library/CloudStorage/GoogleDrive-pi@elewalabs.com/My Drive/01_projects/074_augsburg/ortholist_master.csv',header=TRUE)
ortho <- ortho[ortho$WormBase.ID%in%rownames(wbpheno),]
dim(ortho)
#[1] 17334    11
w2h <- ortho[,c(1,6)]
colnames(w2h) <- c('worm','human')
w2h <- unique(w2h)

hpo <- read.csv('/Users/ahmedelewa/Library/CloudStorage/GoogleDrive-pi@elewalabs.com/My Drive/01_projects/074_augsburg/phenotype_to_genes.csv',header=TRUE)
hpo <- hpo[hpo$gene_symbol%in%ortho$HGNC.Symbol,]
dim(hpo)
#[1] 627480      5
hpo <- hpo[,c(4,2)]
colnames(hpo) <- c('ID','phenotype')
hpo <- unique(hpo)



poi <- c('Autistic behavior')
# Function to perform hypergeometric test
find_phenologs <- function(hpo, w2h, wpo) {
  results <- list()
  
  # Get the first 100 unique human phenotypes
  #human_phenotypes <- unique(hpo$phenotype)[101:1000]
  human_phenotypes <- poi 
  for (h_phenotype in human_phenotypes) {
    # Get human genes associated with this phenotype
    human_genes <- hpo$ID[hpo$phenotype == h_phenotype]
    
    # Map human genes to worm orthologs
    worm_genes_from_human <- unique(w2h$worm[w2h$human %in% human_genes])
    
    # Number of orthologs linked to the human phenotype
    n <- length(worm_genes_from_human)
    
    # Perform analysis for each worm phenotype
    for (w_phenotype in unique(wpo$phenotype)) {
      # Get worm genes associated with this worm phenotype
      worm_genes_in_pheno <- unique(wpo$ID[wpo$phenotype == w_phenotype])
      
      # Number of orthologs linked to the worm phenotype
      m <- length(intersect(w2h$worm, worm_genes_in_pheno))
      
      # Number of common orthologs between human and worm phenotypes
      c <- length(intersect(worm_genes_from_human, worm_genes_in_pheno))
      
      # Total number of orthologs shared between the two species
      N <- nrow(w2h)
      
      # Perform the hypergeometric test
      p_value <- phyper(c - 1, m, N - m, n, lower.tail = FALSE)
      
      # Store significant results
      if (p_value < 0.05) {
        results <- append(results, list(data.frame(
          human_phenotype = h_phenotype,
          worm_phenotype = w_phenotype,
          N = N,
          n = n,
          m = m,
          c = c,
          p_value = p_value
        )))
      }
    }
  }
  
  # Combine results into a single data frame
  do.call(rbind, results)
}

# Run the function for the first 100 human phenotypes
phenolog_results <- find_phenologs(hpo, w2h, wpo)

write.csv(phenolog_results, file = "/Users/ahmedelewa/Library/CloudStorage/GoogleDrive-pi@elewalabs.com/My Drive/01_projects/074_augsburg/phenolog_results_Autisticbehavior.csv", row.names = FALSE)









#poi <- c("Abnormal limb bone morphology", "Abnormal digit morphology", "Abnormal finger morphology", "Abnormal finger phalanx morphology", "Abnormal distal phalanx morphology of finger", "Broad distal phalanx of finger", "Broad distal phalanges of all fingers", "Radial dysplasia", "Abnormal morphology of the radius", "Abnormality of the upper limb", "Abnormal forearm bone morphology", "Abnormal upper limb bone morphology")
poi <- c("Broad distal phalanx of finger", "Broad distal phalanges of all fingers", "Radial dysplasia")
# Function to calculate phenologs with FDR correction
find_phenologs_with_fdr <- function(hpo, w2h, wpo, num_permutations = 10, fdr_threshold = 0.05) { #1000
  results <- list()
  permuted_pvalues <- c()  # Store all permuted p-values
  
  # Get unique human phenotypes (first 100 as required)
  #human_phenotypes <- unique(hpo$phenotype)[1:10]
  human_phenotypes <- poi[1:3] 
  for (h_phenotype in human_phenotypes) {
    # Get human genes associated with this phenotype
    human_genes <- hpo$ID[hpo$phenotype == h_phenotype]
    
    # Map human genes to worm orthologs
    worm_genes_from_human <- unique(w2h$worm[w2h$human %in% human_genes])
    
    # Number of orthologs linked to the human phenotype
    n <- length(worm_genes_from_human)
    
    # Perform analysis for each worm phenotype
    for (w_phenotype in unique(wpo$phenotype)) {
      # Get worm genes associated with this worm phenotype
      worm_genes_in_pheno <- unique(wpo$ID[wpo$phenotype == w_phenotype])
      
      # Number of orthologs linked to the worm phenotype
      m <- length(intersect(w2h$worm, worm_genes_in_pheno))
      
      # Number of common orthologs between human and worm phenotypes
      c <- length(intersect(worm_genes_from_human, worm_genes_in_pheno))
      
      # Total number of orthologs shared between the two species
      N <- nrow(w2h)
      
      # Perform the hypergeometric test
      p_value <- phyper(c - 1, m, N - m, n, lower.tail = FALSE)
      
      # Store real results
      results <- append(results, list(data.frame(
        human_phenotype = h_phenotype,
        worm_phenotype = w_phenotype,
        N = N,
        n = n,
        m = m,
        c = c,
        p_value = p_value
      )))
      
      # Perform permutations to generate null p-values
      for (i in 1:num_permutations) {
        # Shuffle worm orthologs within the phenotype
        shuffled_worm_genes <- sample(w2h$worm, n)  # Shuffle worm orthologs
        
        # Calculate overlap with worm phenotype genes
        shuffled_c <- length(intersect(shuffled_worm_genes, worm_genes_in_pheno))
        
        # Compute p-value for the shuffled data
        permuted_pvalue <- phyper(shuffled_c - 1, m, N - m, n, lower.tail = FALSE)
        permuted_pvalues <- c(permuted_pvalues, permuted_pvalue)
      }
    }
  }
  
  # Combine real results into a data frame
  real_results <- do.call(rbind, results)
  
  # Rank real and permuted p-values
  all_pvalues <- c(real_results$p_value, permuted_pvalues)
  real_labels <- c(rep(TRUE, nrow(real_results)), rep(FALSE, length(permuted_pvalues)))
  ranked_data <- data.frame(p_value = all_pvalues, is_real = real_labels)
  ranked_data <- ranked_data[order(ranked_data$p_value), ]
  
  # Calculate FDR for each threshold
  ranked_data$cum_permuted <- cumsum(!ranked_data$is_real)
  ranked_data$cum_real <- cumsum(ranked_data$is_real)
  ranked_data$fdr <- ranked_data$cum_permuted / ranked_data$cum_real
  
  # Find the maximum p-value threshold where FDR <= fdr_threshold
  fdr_cutoff <- max(ranked_data$p_value[ranked_data$fdr <= fdr_threshold], na.rm = TRUE)
  
  # Filter real results based on the FDR cutoff
  significant_results <- real_results[real_results$p_value <= fdr_cutoff, ]
  
  # Return significant results and the FDR cutoff
  list(
    significant_results = significant_results,
    fdr_cutoff = fdr_cutoff
  )
}

# Run the function
fdr_results <- find_phenologs_with_fdr(hpo, w2h, wpo)

# Save the results to a CSV
#write.csv(fdr_results$significant_results, file = "/Users/ahmedelewa/Library/CloudStorage/GoogleDrive-pi@elewalabs.com/My Drive/01_projects/074_augsburg/phenolog_results_with_fdr.csv", row.names = FALSE)
write.csv(fdr_results$significant_results, file = "/Users/ahmedelewa/Library/CloudStorage/GoogleDrive-pi@elewalabs.com/My Drive/01_projects/074_augsburg/phenolog_results_with_fdr_SP.csv", row.names = FALSE)
cat("Results saved to 'phenolog_results_with_fdr.csv'.\n")
cat("FDR cutoff used: ", fdr_results$fdr_cutoff, "\n")
