setwd("C:/Users/JANE/Desktop/Wang實驗室/NMOSD研究計畫/RNA/NMOSD_RNA_analysis")

#========================================【My data "edgeR" DEGs】========================================
library(edgeR)
library(ggplot2)

########################{原始gene matrix}########################
Raw_count_merged <- read.csv("./RNA_DATA/Raw_count_merged_matrix.csv",
                             row.names = 1, header = TRUE, check.names = FALSE)

# 轉 matrix 並移除基因 ID 尾端小數
count_df <- as.matrix(Raw_count_merged)
rownames(count_df) <- gsub("\\.\\d+$", "", rownames(count_df))

########################{Intersection gene list}########################
protein_coding_intersection_genes<-read.csv('./RNA_DATA/External_with_My_protein_coding_genes.csv')$x #19722

inter_count_df<-count_df[rownames(count_df) %in% intersection_genes,] #56657
protein_inter_count_df<-count_df[rownames(count_df) %in% protein_coding_intersection_genes,] #19722

########################{Sample metadata}########################
#Sample & Group
sample_df <- data.frame(Sample = colnames(protein_inter_count_df))
sample_df$Group <- ifelse(grepl("^SRR", sample_df$Sample), "Control", "NMOSD")
sample_df$Group <- factor(sample_df$Group, levels = c("Control", "NMOSD"))
rownames(sample_df) <- sample_df$Sample

# ----【建 DGEList】----
dge <- DGEList(counts = protein_inter_count_df, group = sample_df$Group)

# ----【過濾低表達genes】 ----
# 依組別與j文庫大小自動估門檻（約等 CPM>=1 在足夠多樣本）
keep <- filterByExpr(dge, group = sample_df$Group)
dge  <- dge[keep, , keep.lib.sizes = FALSE]

# ---- 【TMM 標準化】 ----
dge <- calcNormFactors(dge, method = "TMM")

# ---- 【設計矩陣（Control 為 baseline）】----
design <- model.matrix(~ Group, data = sample_df)
colnames(design)
# 係數解讀：Intercept = Control，GroupNMOSD = NMOSD vs Control 對比

# ---- 【估計離散度】----
dge <- estimateDisp(dge, design)

# ----【 GLM 擬合 + QLF 檢定】----
#QLF適合"樣本數不大、資料變異較大（異質性高）"的情況，假陽性率控制得比較好
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = "GroupNMOSD")  # NMOSD vs Control


# ---- 結果表 ----
table <- topTags(qlf, n = Inf)$table
# 欄位包含：logFC, logCPM, F, PValue, FDR（Benjamini-Hochberg）


res <- data.frame(
  log2FoldChange = table$logFC,
  pvalue = table$PValue,
  padj = table$FDR,
  row.names = rownames(table),
  check.names = FALSE
)

# 排序（依 padj 由小到大）
res <- res[order(res$padj), ]

# 總結
summary(decideTests(qlf))  # 上下調基因數量概覽（以 FDR<0.05 默認）

# ------------------{Volcano Plot}------------------
op <- par(no.readonly = TRUE)
par(mar = c(5, 4, 4, 6), xpd = TRUE)

with(res, plot(
  log2FoldChange, -log10(pvalue),
  pch  = 20,
  main = "NMOSD vs Healthy DEGs (edgeR)",
  xlim = c(-15, 15),
  col  = "grey",
  xlab = "log2 Fold Change",
  ylab = "-log10(P Value)"
))

#----【Significants genes criteria: padj < 0.05 & |log2FoldChange|>1】-----
with(subset(res, padj < 0.05 & log2FoldChange < -1),
     points(log2FoldChange, -log10(pvalue), pch = 20, col = "#8FA4FF"))

with(subset(res, padj < 0.05 & log2FoldChange > 1),
     points(log2FoldChange, -log10(pvalue), pch = 20, col = "#FF8080"))

legend("topright",
       inset  = c(-0.2, 0),
       legend = c("Up", "Down"),
       title  = "Change",
       pch    = 20,
       col    = c("#FF8080", "#8FA4FF"),
       pt.cex = 1.4,
       bty    = "n")
par(op)

# ---- 顯著基因----
Up_DEGs   <- row.names(subset(res, padj < 0.05 & log2FoldChange >  1))
Down_DEGs <- row.names(subset(res, padj < 0.05 & log2FoldChange < -1))
DEGS      <- row.names(subset(res, padj < 0.05 & abs(log2FoldChange) > 1))

length(Up_DEGs); length(Down_DEGs); length(DEGS)

# 輸出 CSV
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS(my data)(edgeR).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP(my data)(edgeR).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN(my data)(edgeR).csv',row.names = F)

#----【Significants genes criteria: padj < 0.01 & |log2FoldChange|>2】-----
op <- par(no.readonly = TRUE)
par(mar = c(5, 4, 4, 6), xpd = TRUE)

with(res, plot(
  log2FoldChange, -log10(pvalue),
  pch  = 20,
  main = "NMOSD vs Healthy DEGs (edgeR)",
  xlim = c(-15, 15),
  col  = "grey",
  xlab = "log2 Fold Change",
  ylab = "-log10(P Value)"
))

with(subset(res, padj < 0.01 & log2FoldChange < -2),
     points(log2FoldChange, -log10(pvalue), pch = 20, col = "#8FA4FF"))

with(subset(res, padj < 0.01 & log2FoldChange > 2),
     points(log2FoldChange, -log10(pvalue), pch = 20, col = "#FF8080"))

legend("topright",
       inset  = c(-0.2, 0),
       legend = c("Up", "Down"),
       title  = "Change",
       pch    = 20,
       col    = c("#FF8080", "#8FA4FF"),
       pt.cex = 1.4,
       bty    = "n")
par(op)

# ---- 顯著基因----
Up_DEGs   <- row.names(subset(res, padj < 0.01 & log2FoldChange >  2))
Down_DEGs <- row.names(subset(res, padj < 0.01 & log2FoldChange < -2))
DEGS      <- row.names(subset(res, padj < 0.01 & abs(log2FoldChange) > 2))

length(Up_DEGs); length(Down_DEGs); length(DEGS)

#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS_strictly(my data)(edgeR).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP_strictly(my data)(edgeR).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN_strictly(my data)(edgeR).csv',row.names = F)



#========================================【External data "edgeR" DEGs】========================================
#----{Package}----
library(edgeR)
library(ggplot2)

########################{原始gene matrix}########################
Raw_count_merged<-read.csv('./RNA_DATA/External_RNA_datasets_DEG_RESULT/RNA_gene_matrix_v2.csv',row.names = 1,header=T) #dim(Raw_count_merged):58304 x 23
Raw_count_merged <- Raw_count_merged[!rownames(Raw_count_merged) %in% c("Antibody status","Ensembl_Gene_ID"), ]
#dim(Raw_count_merged) 58302 x 23

#Raw count 
count_df<-as.matrix(Raw_count_merged)
storage.mode(count_df) <- "numeric"

########################{Intersection gene list}########################
protein_coding_intersection_genes<-read.csv('./RNA_DATA/External_with_My_protein_coding_genes.csv')$x #19722
protein_inter_count_df<-count_df[rownames(count_df) %in% protein_coding_intersection_genes,] #19722

########################{Sample metadata}########################
#Sample & Group
sample_df=data.frame(Sample=colnames(protein_inter_count_df))
sample_df$Group <- ifelse(
  grepl("NA", sample_df$Sample),  # 如果 Sample 以 "NA" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
# 把 sample_df 裡的 Group 欄轉成 factor，Control 當作 baseline
sample_df$Group <- factor(sample_df$Group,
                          levels = c("Control","NMOSD"))
rownames(sample_df) <- sample_df$Sample

# ----【建 DGEList】----
dge <- DGEList(counts = protein_inter_count_df, group = sample_df$Group)

# ----【過濾低表達genes】 ----
# 依組別與j文庫大小自動估門檻（約等 CPM>=1 在足夠多樣本）
keep <- filterByExpr(dge, group = sample_df$Group)
dge  <- dge[keep, , keep.lib.sizes = FALSE]

# ---- 【TMM 標準化】 ----
dge <- calcNormFactors(dge, method = "TMM")

# ---- 【設計矩陣（Control 為 baseline）】----
design <- model.matrix(~ Group, data = sample_df)
colnames(design)
# 係數解讀：Intercept = Control，GroupNMOSD = NMOSD vs Control 對比

# ---- 【估計離散度】----
dge <- estimateDisp(dge, design)

# ----【 GLM 擬合 + QLF 檢定】----
#QLF適合"樣本數不大、資料變異較大（異質性高）"的情況，假陽性率控制得比較好
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = "GroupNMOSD")  # NMOSD vs Control


# ---- 結果表 ----
table <- topTags(qlf, n = Inf)$table
# 欄位包含：logFC, logCPM, F, PValue, FDR（Benjamini-Hochberg）


res <- data.frame(
  log2FoldChange = table$logFC,
  pvalue = table$PValue,
  padj = table$FDR,
  row.names = rownames(table),
  check.names = FALSE
)

# 排序（依 padj 由小到大）
res <- res[order(res$padj), ]

# 總結
summary(decideTests(qlf))  # 上下調基因數量概覽（以 FDR<0.05 默認）

# ------------------{Volcano Plot}------------------
op <- par(no.readonly = TRUE)
par(mar = c(5, 4, 4, 6), xpd = TRUE)

with(res, plot(
  log2FoldChange, -log10(pvalue),
  pch  = 20,
  main = "NMOSD vs Healthy edgeR DEGs  (External Dataset)",
  xlim = c(-15, 15),
  col  = "grey",
  xlab = "log2 Fold Change",
  ylab = "-log10(P Value)"
))

#----【Significants genes criteria: padj < 0.05 & |log2FoldChange|>1】-----
with(subset(res, padj < 0.05 & log2FoldChange < -1),
     points(log2FoldChange, -log10(pvalue), pch = 20, col = "#8FA4FF"))

with(subset(res, padj < 0.05 & log2FoldChange > 1),
     points(log2FoldChange, -log10(pvalue), pch = 20, col = "#FF8080"))

legend("topright",
       inset  = c(-0.2, 0),
       legend = c("Up", "Down"),
       title  = "Change",
       pch    = 20,
       col    = c("#FF8080", "#8FA4FF"),
       pt.cex = 1.4,
       bty    = "n")
par(op)

# ---- 顯著基因----
Up_DEGs   <- row.names(subset(res, padj < 0.05 & log2FoldChange >  1))
Down_DEGs <- row.names(subset(res, padj < 0.05 & log2FoldChange < -1))
DEGS      <- row.names(subset(res, padj < 0.05 & abs(log2FoldChange) > 1))

length(Up_DEGs); length(Down_DEGs); length(DEGS)

# 輸出 CSV
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS(external data)(edgeR).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP(external data)(edgeR).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN(external data)(edgeR).csv',row.names = F)

#----【Significants genes criteria: padj < 0.01 & |log2FoldChange|>2】-----
op <- par(no.readonly = TRUE)
par(mar = c(5, 4, 4, 6), xpd = TRUE)

with(res, plot(
  log2FoldChange, -log10(pvalue),
  pch  = 20,
  main = "NMOSD vs Healthy edgeR DEGs  (External Dataset)",
  xlim = c(-15, 15),
  col  = "grey",
  xlab = "log2 Fold Change",
  ylab = "-log10(P Value)"
))

with(subset(res, padj < 0.01 & log2FoldChange < -2),
     points(log2FoldChange, -log10(pvalue), pch = 20, col = "#8FA4FF"))

with(subset(res, padj < 0.01 & log2FoldChange > 2),
     points(log2FoldChange, -log10(pvalue), pch = 20, col = "#FF8080"))

legend("topright",
       inset  = c(-0.2, 0),
       legend = c("Up", "Down"),
       title  = "Change",
       pch    = 20,
       col    = c("#FF8080", "#8FA4FF"),
       pt.cex = 1.4,
       bty    = "n")
par(op)

# ---- 顯著基因----
Up_DEGs   <- row.names(subset(res, padj < 0.01 & log2FoldChange >  2))
Down_DEGs <- row.names(subset(res, padj < 0.01 & log2FoldChange < -2))
DEGS      <- row.names(subset(res, padj < 0.01 & abs(log2FoldChange) > 2))

length(Up_DEGs); length(Down_DEGs); length(DEGS)

#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS_strictly(external data)(edgeR).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP_strictly(external data)(edgeR).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN_strictly(external data)(edgeR).csv',row.names = F)

#########################################【Function 計算IFN數量】#########################################
library(clusterProfiler)

#=====================================================<My Dataset>=====================================================
##################{My dataset "edgeR" protein coding-DEGs (criteria: padjust<0.05 & |log2FC|>1}##################
#---{Interferon count}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(my data)(edgeR).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:2886個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:50個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE)
interferon_df<-data.frame(Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
print(interferon_df)

#---{type.I.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(my data)(edgeR).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
length(interferon_in_DEGs)
print(interferon_in_DEGs)
#write.csv(interferon_in_DEGs,file = "./RNA_DATA/Type_I_IFN_DEG_genes(my data)(edgeR).csv",row.names = FALSE)

#---{hallmark.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(my data)(edgeR).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
length(hallmark_interferon_in_DEGs)
print(hallmark_interferon_in_DEGs)
#write.csv(hallmark_interferon_in_DEGs,file= "./RNA_DATA/Hallmark_IFN_DEG_genes(my data)(edgeR).csv",row.names = FALSE)


##################{My dataset "edgeR" protein coding-DEGs (criteria: padjust<0.01 & |log2FC|>2}##################
#---{Interferon count}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(my data)(edgeR).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:2886個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:50個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE)
interferon_df<-data.frame(Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
print(interferon_df)

#---{type.I.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(my data)(edgeR).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
length(interferon_in_DEGs)
print(interferon_in_DEGs)
#write.csv(interferon_in_DEGs,file = "./RNA_DATA/Type_I_IFN_DEG_genes_strictly(my data)(edgeR).csv",row.names = FALSE)

#---{hallmark.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(my data)(edgeR).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
length(hallmark_interferon_in_DEGs)
print(hallmark_interferon_in_DEGs)
#write.csv(hallmark_interferon_in_DEGs,file= "./RNA_DATA/Hallmark_IFN_DEG_genes_strictly(my data)(edgeR).csv",row.names = FALSE)


#=====================================================<External Dataset>=====================================================
##################{External dataset "edgeR" protein coding-DEGs (criteria: padjust<0.05 & |log2FC|>1}##################
#---{Interferon count}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external data)(edgeR).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:2886個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:50個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE)
interferon_df<-data.frame(Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
print(interferon_df)

#---{type.I.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external data)(edgeR).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
length(interferon_in_DEGs)
print(interferon_in_DEGs)
#write.csv(interferon_in_DEGs,file = "./RNA_DATA/Type_I_IFN_DEG_genes(external data)(edgeR).csv",row.names = FALSE)

#---{hallmark.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external data)(edgeR).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
length(hallmark_interferon_in_DEGs)
print(hallmark_interferon_in_DEGs)
#write.csv(hallmark_interferon_in_DEGs,file= "./RNA_DATA/Hallmark_IFN_DEG_genes(external data)(edgeR).csv",row.names = FALSE)


##################{External dataset "edgeR" protein coding-DEGs (criteria: padjust<0.01 & |log2FC|>2}##################
#---{Interferon count}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(external data)(edgeR).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:50個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE)
interferon_df<-data.frame(Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
print(interferon_df)

#---{type.I.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(external data)(edgeR).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
length(interferon_in_DEGs)
print(interferon_in_DEGs)
#write.csv(interferon_in_DEGs,file = "./RNA_DATA/Type_I_IFN_DEG_genes_strictly(external data)(edgeR).csv",row.names = FALSE)

#---{hallmark.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(external data)(edgeR).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
length(hallmark_interferon_in_DEGs)
print(hallmark_interferon_in_DEGs)
#write.csv(hallmark_interferon_in_DEGs,file= "./RNA_DATA/Hallmark_IFN_DEG_genes_strictly(external data)(edgeR).csv",row.names = FALSE)

