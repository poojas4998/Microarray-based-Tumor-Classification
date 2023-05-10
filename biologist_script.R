# Load required packages 
library(hgu133plus2.db)
library(tidyverse)
library(AnnotationDbi)
BiocManager::install("GSEABase") # Install or update the GSEABase package
library(GSEABase)
library('stats')

#PROBLEM 1: Map probe IDs to gene symbols 

# Load data file and rename column names
diff_exp_stats <- readr::read_csv('/projectnb/bf528/users/group_5/project_1/project_1_group_5/expression_data/part_five_passed_two_filters_statistics.csv') %>%
  rename("PROBEID" = "...1")

# Match probe IDs to gene symbols, merge to dataframe containing expression values and drop repeats 
probesets <- AnnotationDbi::select(hgu133plus2.db, keys = diff_exp_stats$PROBEID,
                                   columns = "SYMBOL",
                                   keytype = "PROBEID")
merged_data <- full_join(probesets, diff_exp_stats, by="PROBEID")
gene_exp_unique <- subset(merged_data, !duplicated(SYMBOL)) 

#PROBLEM 3: Find top 1000 up and down-regulated genes 

# All genes with adjusted p-value <0.05 have significant differential expression 
degs <- gene_exp_unique[gene_exp_unique$p_adjusted < 0.05, ]
upregulated <- degs[order(degs$t_statistics, decreasing=TRUE), ]
up_1000 <- head(upregulated, 1000)
downregulated <- degs[order(degs$t_statistics, decreasing=FALSE), ]
down_1000 <- head(downregulated, 1000)

# To find top 10 upregulated and downregulated genes 
top_10_up <- head(up_1000, 10)
top_10_down <- head(down_1000, 10)

#PROBLEM 4: Load gene sets

#Use GSEABase to load gene sets 
kegg_symbols <- GSEABase::getGmt("/projectnb/bf528/users/group_5/project_1/project_1_group_5/Biologist/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
hallmark_symbols <- GSEABase::getGmt("/projectnb/bf528/users/group_5/project_1/project_1_group_5/Biologist/h.all.v2022.1.Hs.symbols.gmt")
go_symbols <- GSEABase::getGmt("/projectnb/bf528/users/group_5/project_1/project_1_group_5/Biologist/c5.go.v2022.1.Hs.symbols.gmt")

#PROBLEM 5: Use fisher test to compare overlap for each gene set and each of the top 1000 increased and 1000 decreased genes

#create function to make contingency table 
make_contingency <- function(dg_list, GeneSetobj) {
  #convert dg_list into GeneSet using GSEABase
  dg_geneset <- GeneSet(dg_list, setName='1000 DE genes')
  #total number of genes in initial data file 
  bg_len <- length(gene_exp_unique$SYMBOL)
  
  a <- length((dg_geneset & GeneSetobj)@geneIds) #no. of degs present gene set
  b <- length(setdiff(dg_geneset, GeneSetobj)@geneIds) #no. of degs not in gene set  
  c <- length(setdiff(GeneSetobj, dg_geneset)@geneIds) #no. of genes not degs but present in gene sets  
  d <- bg_len - (a + b + c) #no. of genes not degs and not in gene sets

  #return a matrix of the contingency values
  return(matrix(c(a, b, c, d), nrow=2))
}

#initialising matrix to store overlapping genes in KEGG gene set
KEGG_overlap <- data.frame(matrix(ncol = 3, nrow = 0))

colnames(KEGG_overlap) <- c('set', 'up_1000', 'down_1000')
KEGG_overlap$set <- as.character(KEGG_overlap$set)
KEGG_overlap$up_1000 <- as.numeric(KEGG_overlap$up_1000)
KEGG_overlap$down_1000 <- as.numeric(KEGG_overlap$down_1000)

#populating KEGG overlap matrix 
for (set in kegg_symbols)
{
  cont_table <- make_contingency(up_1000$SYMBOL, set)
  fisher_results <- fisher.test(cont_table)
  up <- fisher_results$p.value
  cont_table <- make_contingency(down_1000$SYMBOL, set)
  fisher_results <- fisher.test(cont_table)
  down <- fisher_results$p.value
  KEGG_overlap[nrow(KEGG_overlap) + 1,] = c(set@setName, up, down)
}

#initialising matrix to store overlapping genes in GO gene set
GO_overlap <- data.frame(matrix(ncol = 3, nrow = 0))

colnames(GO_overlap) <- c('set', 'up_1000', 'down_1000')
GO_overlap$set <- as.character(GO_overlap$set)
GO_overlap$up_1000 <- as.numeric(GO_overlap$up_1000)
GO_overlap$down_1000 <- as.numeric(GO_overlap$down_1000)

#populating GO overlap matrix
for (set in go_symbols)
{
  cont_table <- make_contingency(up_1000$SYMBOL, set)
  fisher_results <- fisher.test(cont_table)
  up <- fisher_results$p.value
  cont_table <- make_contingency(down_1000$SYMBOL, set)
  fisher_results <- fisher.test(cont_table)
  down <- fisher_results$p.value
  
  GO_overlap[nrow(GO_overlap) + 1,] = c(set@setName, up, down)
}

#initialising matrix to store overlapping genes in Hallmark gene set
Hallmark_overlap <- data.frame(matrix(ncol = 3, nrow = 0))

colnames(Hallmark_overlap) <- c('set', 'up_1000', 'down_1000')
Hallmark_overlap$set <- as.character(Hallmark_overlap$set)
Hallmark_overlap$up_1000 <- as.numeric(Hallmark_overlap$up_1000)
Hallmark_overlap$down_1000 <- as.numeric(Hallmark_overlap$down_1000)

#populating Hallmark overlap matrix
for (set in hallmark_symbols)
{
  cont_table <- make_contingency(up_1000$SYMBOL, set)
  fisher_results <- fisher.test(cont_table)
  up <- fisher_results$p.value
  cont_table <- make_contingency(down_1000$SYMBOL, set)
  fisher_results <- fisher.test(cont_table)
  down <- fisher_results$p.value
  
  Hallmark_overlap[nrow(Hallmark_overlap) + 1,] = c(set@setName, up, down)
}

#PROBLEM 6: Create DFs containing enriched gene sets 

#splitting KEGG overlap into two subsets containing upregulated and downregulated genes, 
#calculating p adjusted values for both and sorting in ascending order of p-value 
k_up <- select(KEGG_overlap, c("set","up_1000")) %>% rename('p'='up_1000')
k_up['p_adj']<- stats::p.adjust(k_up$p, method = 'fdr')
k_up <- k_up[order(k_up$p),]

k_down <- select(KEGG_overlap, c("set","down_1000")) %>% rename('p'='down_1000')
k_down['p_adj']<- stats::p.adjust(k_down$p, method = 'fdr')
k_down <- k_down[order(k_down$p),]

#Repeat process for Hallmark gene set 
h_up <- select(Hallmark_overlap, c("set","up_1000")) %>% rename('p'='up_1000')
h_up['p_adj']<- stats::p.adjust(h_up$p, method = 'fdr')
h_up <- h_up[order(h_up$p),]

h_down <- select(Hallmark_overlap, c("set","down_1000")) %>% rename('p'='down_1000')
h_down['p_adj']<- stats::p.adjust(h_down$p, method = 'fdr')
h_down <- h_down[order(h_down$p),]

#Repeating process for GO gene set 
g_up <- select(GO_overlap, c("set","up_1000")) %>% rename('p'='up_1000')
g_up['p_adj']<- stats::p.adjust(g_up$p, method = 'fdr')
g_up <- g_up[order(g_up$p),] %>% head(1000)

g_down <- select(GO_overlap, c("set","down_1000")) %>% rename('p'='down_1000')
g_down['p_adj']<- stats::p.adjust(g_down$p, method = 'fdr')
g_down <- g_down[order(g_down$p),] %>% head(1000)

