install.packages("remotes")
remotes::install_github("omnideconv/immunedeconv")
library(immunedeconv)

#-------【Bulk RNA seq data(TPM)】-------

#{Discovery cohort-TPM}
TPM_df<-read.csv("./RNA_DATA/TPM_gene_samples_matrix.csv")
#dim(TPM_df) gene x sample: 78724 x 22('gene_id+samples)

#{Independent cohort-TPM}
TPM_df<-read.csv("./RNA_DATA/External_TPM_gene_samples_matrix(intersection gene).csv")

# 篩選出欄位名稱含 "NA" 或 "no"
keep_cols <- grep("NA|no", colnames(TPM_df), value = TRUE)
TPM_df<- TPM_df[, c("gene_id", keep_cols)]
#dim(TPM_df) 19722 x 24


#-----【ENSEMBL ID -> HGNC】-------- 
library(biomaRt)

rownames(TPM_df)<-TPM_df$gene_id
TPM_df <- TPM_df[,colnames(TPM_df) != "gene_id"]
ens_id <- rownames(TPM_df)
ens_id_clean <- sub("\\..*$", "", ens_id)

mart <- useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl"
)

map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = ens_id_clean,
  mart = mart
)
colnames(map) <- c("ENSEMBL", "SYMBOL", "gene_biotype")

expr_df <- as.data.frame(TPM_df)
expr_df$ENSEMBL <- ens_id_clean

library(dplyr)

expr_symbol <- expr_df %>%
  left_join(map, by = "ENSEMBL") %>%
  filter(!is.na(SYMBOL), SYMBOL != "") #刪除沒有symbol
#僅剩41954 rows

#{把symbol一樣的做相加}
expr_symbol_mat <- expr_symbol %>%
  select(-ENSEMBL, -gene_biotype) %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(), sum), .groups = "drop")

expr_tpm_symbol <- as.data.frame(expr_symbol_mat)
rownames(expr_tpm_symbol) <- expr_tpm_symbol$SYMBOL
expr_tpm_symbol$SYMBOL <- NULL
expr_tpm_symbol <- as.matrix(expr_tpm_symbol)

#----【計算Cell type fraction】-----
cell_type_df<-immunedeconv::deconvolute(expr_tpm_symbol, "quantiseq")
#discovery_Bcells_comp<-cell_type_df[cell_type_df$cell_type=='B cell',]
#discovery_Bcells_comp
#write.csv(discovery_Bcells_comp,
#          file = "./RNA_DATA/discovery_Bcells_comp.csv",
#          row.names = FALSE)

independent_Bcells_comp<-cell_type_df[cell_type_df$cell_type=='B cell',]
independent_Bcells_comp
write.csv(independent_Bcells_comp,
  file = "./RNA_DATA/independent_Bcells_comp.csv",
  row.names = FALSE)


#----【scRNA B cell types比例】----
celltype_metadata <- read.csv("C:/Users/Jane/Desktop/Wang實驗室/NMOSD研究計畫/scRNA/scRNA_DATA/My_merged_protein_coding_genes/My_merged_Azimuth_sub_celltype_CellChat.csv")
head(celltype_metadata)

B_celltypes <- c("B intermediate", "B memory", "B naive")
B_ratio_per_sample <- celltype_metadata %>%
  mutate(sample = sub("_.*", "", cell_id)) %>%
  group_by(sample) %>%
  summarise(
    total_cells = n(),
    B_cell_number = sum(sub_celltype %in% B_celltypes),
    B_cell_fraction = B_cell_number / total_cells,
    .groups = "drop"
  )
B_ratio_per_sample

#轉成列是B-cell，欄是samples
library(tidyr)
B_fraction_wide <- B_ratio_per_sample %>%
  select(sample, B_cell_fraction) %>%
  mutate(cell_type = "B cell") %>%
  pivot_wider(
    names_from = sample,
    values_from = B_cell_fraction
  ) %>%
  select(cell_type, everything())

scRNA_B_comp<-B_fraction_wide
write.csv(scRNA_B_comp,
          file = "./RNA_DATA/scRNA_Bcells_comp.csv",
          row.names = FALSE)

#----【Bulk 和 scRNA作統計檢定】----

#合併discovery & independent B-cell composition
bulk_B_comp <- full_join(
  discovery_Bcells_comp,
  independent_Bcells_comp,
  by = "cell_type"
)
bulk_B_comp

#轉成長格式
bulk_B_long <- bulk_B_comp %>%
  filter(cell_type == "B cell") %>%
  pivot_longer(
    cols = -cell_type,
    names_to = "sample",
    values_to = "B_fraction"
  ) %>%
  mutate(source = "bulk")

scrna_B_long <- scRNA_B_comp %>%
  filter(cell_type == "B cell") %>%
  pivot_longer(
    cols = -cell_type,
    names_to = "sample",
    values_to = "B_fraction"
  ) %>%
  mutate(source = "scRNA")

#合併bulk和scrna長格式
B_compare <- bind_rows(bulk_B_long, scrna_B_long)
B_compare
B_compare_Pvalue<-wilcox.test(B_fraction ~ source, data = B_compare)$p.value #unpaired test
print(B_compare_Pvalue)

#畫boxplot
boxplot_Bcell <- ggplot(B_compare, aes(x = source, y =B_fraction, fill =  source)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.7) +  # 添加點狀數據
  #theme_minimal() +
  theme_bw()+
  labs(title =" B cell composition between bulk and scRNA" ,
       x=NULL,y = "B_fraction",fill="source") +
  theme(plot.title = element_text(size = 15, hjust = 0.5,face='bold'),
        axis.text.x = element_text(size=13, angle=0,hjust = 0.5,vjust =1), #hjust=0.5置中對齊
        axis.title.y = element_text(size = 13),
        legend.background=element_rect(
          colour='white',fill="white"),
        legend.key.width = unit(0.8, "cm"),
        legend.key.height = unit(0.6, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))+
  annotate(
    "text",
    x = 1.5,
    y = max(B_compare$B_fraction, na.rm = TRUE) * 1.1,
    label = paste0("Padj = ",signif(B_compare_Pvalue, 3)),
    size = 4
  ) +
  scale_fill_manual(values = c("bulk" = "#DAC9EE", "scRNA" = "#A8F6BB"))

print(boxplot_Bcell)

