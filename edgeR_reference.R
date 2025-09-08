setwd("C:/Users/JANE/Desktop/Wang實驗室/NMOSD研究計畫/RNA/NMOSD_RNA_analysis")

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

########################{Sample metadata}########################
#Sample & Group
sample_df=data.frame(Sample=colnames(count_df))
sample_df$Group <- ifelse(
  grepl("NA", sample_df$Sample),  # 如果 Sample 以 "NA" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
# 把 sample_df 裡的 Group 欄轉成 factor，Control 當作 baseline
sample_df$Group <- factor(sample_df$Group,
                          levels = c("Control","NMOSD"))
rownames(sample_df) <- sample_df$Sample

# ---- 建 DGEList ----
dge <- DGEList(counts = count_df, group = sample_df$Group)

#----edgeR過濾低expression genes----
keep <- filterByExpr(dge, group = sample_df$Group)
dge  <- dge[keep, , keep.lib.sizes = FALSE]

# ---- TMM 標準化 ----
dge <- calcNormFactors(dge, method = "TMM")

# ---- 設計矩陣（Control 為 baseline）----
design <- model.matrix(~ Group, data = sample_df)
colnames(design)
# 係數解讀：Intercept = Control，GroupNMOSD = NMOSD vs Control 對比

# ---- 估計離散度----
dge <- estimateDisp(dge, design)

# ---- GLM 擬合 + QLF 檢定----
#QLF適合"樣本數不大、資料變異較大（異質性高）"的情況，假陽性率控制得比較好
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = "GroupNMOSD")  # NMOSD vs Control


# ---- 完整結果表 ----
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
Up_DEGs   <- row.names(subset(res, padj < 0.05 & log2FoldChange >  1))#19
Down_DEGs <- row.names(subset(res, padj < 0.05 & log2FoldChange < -1))#24
DEGs      <- row.names(subset(res, padj < 0.05 & abs(log2FoldChange) > 1))#43

length(Up_DEGs); length(Down_DEGs); length(DEGs)

# 輸出 CSV
#write.csv(DEGs,      file = "RNA_DATA/NMOSD_RNA_DEGs(external)(edgeR).csv", row.names = FALSE)
#write.csv(Up_DEGs,   file = "RNA_DATA/NMOSD_RNA_DEGs_UP(external)(edgeR).csv", row.names = FALSE)
#write.csv(Down_DEGs, file = "RNA_DATA/NMOSD_RNA_DEGs_DOWN(external)(edgeR).csv", row.names = FALSE)

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
Up_DEGs   <- row.names(subset(res, padj < 0.01 & log2FoldChange >  2))#4
Down_DEGs <- row.names(subset(res, padj < 0.01 & log2FoldChange < -2))#1
DEGs      <- row.names(subset(res, padj < 0.01 & abs(log2FoldChange) > 2))#5

length(Up_DEGs); length(Down_DEGs); length(DEGs)

# 輸出 CSV
#write.csv(DEGs,      file = "RNA_DATA/NMOSD_RNA_DEGs_strictly(external)(edgeR).csv", row.names = FALSE)
#write.csv(Up_DEGs,   file = "RNA_DATA/NMOSD_RNA_DEGs_UP_strictly(external)(edgeR).csv", row.names = FALSE)
#write.csv(Down_DEGs, file = "RNA_DATA/NMOSD_RNA_DEGs_DOWN_strictly(external)(edgeR).csv", row.names = FALSE)



#################################【篩"proteing coding genes"的 DEGs】#################################
Raw_count_merged<-read.csv('./RNA_DATA/External_RNA_datasets_DEG_RESULT/RNA_gene_matrix_v2.csv',row.names = 1,header=T) #dim(Raw_count_merged):58304 x 23
Raw_count_merged <- Raw_count_merged[!rownames(Raw_count_merged) %in% c("Antibody status","Ensembl_Gene_ID"), ]
#dim(Raw_count_merged) 58302 x 23

#Raw count 
count_df<-as.matrix(Raw_count_merged)
storage.mode(count_df) <- "numeric"

#----連到Ensembl
library(biomaRt)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

#----My gene matrix
All_genes<-rownames(count_df)


#----抓註解
annot <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters    = "ensembl_gene_id",
  values     = All_genes,
  mart       = ensembl
)

#----篩 protein coding genes
protein_coding_genes <- annot[annot$gene_biotype == "protein_coding", ] #20091個

#----proteing coding df
protein_coding_count_df<-count_df[rownames(count_df) %in% protein_coding_genes$ensembl_gene_id,]

########################{Sample metadata}########################
#Sample & Group
sample_df=data.frame(Sample=colnames(protein_coding_count_df))
sample_df$Group <- ifelse(
  grepl("NA", sample_df$Sample),  # 如果 Sample 以 "NA" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
# 把 sample_df 裡的 Group 欄轉成 factor，Control 當作 baseline
sample_df$Group <- factor(sample_df$Group,
                          levels = c("Control","NMOSD"))
rownames(sample_df) <- sample_df$Sample

# ---- 建 DGEList ----
dge <- DGEList(counts = protein_coding_count_df, group = sample_df$Group)

#----edgeR過濾低expression genes----
keep <- filterByExpr(dge, group = sample_df$Group)
dge  <- dge[keep, , keep.lib.sizes = FALSE]

# ---- TMM 標準化 ----
dge <- calcNormFactors(dge, method = "TMM")

# ---- 設計矩陣（Control 為 baseline）----
design <- model.matrix(~ Group, data = sample_df)
colnames(design)
# 係數解讀：Intercept = Control，GroupNMOSD = NMOSD vs Control 對比

# ---- 估計離散度----
dge <- estimateDisp(dge, design)

# ---- GLM 擬合 + QLF 檢定----
#QLF適合"樣本數不大、資料變異較大（異質性高）"的情況，假陽性率控制得比較好
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = "GroupNMOSD")  # NMOSD vs Control


# ---- 完整結果表 ----
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
Up_DEGs   <- row.names(subset(res, padj < 0.05 & log2FoldChange >  1))#19
Down_DEGs <- row.names(subset(res, padj < 0.05 & log2FoldChange < -1))#22
DEGs      <- row.names(subset(res, padj < 0.05 & abs(log2FoldChange) > 1))#41

length(Up_DEGs); length(Down_DEGs); length(DEGs)

# 輸出 CSV
#write.csv(DEGs,      file = "RNA_DATA/NMOSD_RNA_DEGs(external)(edgeR)(protein coding genes).csv", row.names = FALSE)
#write.csv(Up_DEGs,   file = "RNA_DATA/NMOSD_RNA_DEGs_UP(external)(edgeR)(protein coding genes).csv", row.names = FALSE)
#write.csv(Down_DEGs, file = "RNA_DATA/NMOSD_RNA_DEGs_DOWN(external)(edgeR)(protein coding genes).csv", row.names = FALSE)

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
Up_DEGs   <- row.names(subset(res, padj < 0.01 & log2FoldChange >  2))#4
Down_DEGs <- row.names(subset(res, padj < 0.01 & log2FoldChange < -2))#1
DEGs      <- row.names(subset(res, padj < 0.01 & abs(log2FoldChange) > 2))#5

length(Up_DEGs); length(Down_DEGs); length(DEGs)

# 輸出 CSV
#write.csv(DEGs,      file = "RNA_DATA/NMOSD_RNA_DEGs_strictly(external)(edgeR)(protein coding genes).csv", row.names = FALSE)
#write.csv(Up_DEGs,   file = "RNA_DATA/NMOSD_RNA_DEGs_UP_strictly(external)(edgeR)(protein coding genes).csv", row.names = FALSE)
#write.csv(Down_DEGs, file = "RNA_DATA/NMOSD_RNA_DEGs_DOWN_strictly(external)(edgeR)(protein coding genes).csv", row.names = FALSE)




#-------------------------------------------------------------------------------------------------------------------------


#########################################【Function 計算IFN數量】#########################################
library(clusterProfiler)
#---------{Origin "edgeR" DEGs (criteria: padjust<0.05 & |log2FC|>1}---------

#---{Interferon count}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGs(external)(edgeR).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:2886個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:50個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE)
interferon_df<-data.frame(Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
print(interferon_df)

#---{type.I.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGs(external)(edgeR).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
print(interferon_in_DEGs)
length(interferon_in_DEGs)
#write.csv(interferon_in_DEGs,file = "./RNA_DATA/Type_I_IFN_DEG_genes(external)(edgeR).csv",row.names = FALSE)

#---{hallmark.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGs(external)(edgeR).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
print(hallmark_interferon_in_DEGs)
length(hallmark_interferon_in_DEGs)
#write.csv(hallmark_interferon_in_DEGs,file= "./RNA_DATA/Hallmark_IFN_DEG_genes(external)(edgeR).csv",row.names = FALSE)

#---------{Origin protein coding-DEGs (criteria: padjust<0.05 & |log2FC|>1}---------

#---{Interferon count}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGs(external)(edgeR)(protein coding genes).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE)
interferon_df<-data.frame(Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
print(interferon_df)

#---{type.I.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGs(external)(edgeR)(protein coding genes).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
print(interferon_in_DEGs) 
length(interferon_in_DEGs)
#write.csv(interferon_in_DEGs,file = "./RNA_DATA/Type_I_IFN_DEG_genes(external)(protein codin genes)(edgeR).csv",row.names = FALSE)

#---{hallmark.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGs(external)(edgeR)(protein coding genes).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
print(hallmark_interferon_in_DEGs)
length(hallmark_interferon_in_DEGs)
#write.csv(hallmark_interferon_in_DEGs,file= "./RNA_DATA/Hallmark_IFN_DEG_genes(external)(protein codin genes)(edgeR).csv",row.names = FALSE)
#---------------------------------------------------------------------------------------------------------------------


#############################################【edgeR - DEGs interection】#############################################
#---{ORIGIN & EXTERNAL}
library(VennDiagram)
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGs(external)(edgeR).csv',header=T)
DEG_genes<-DEGs_df$x 
DEGs_ori<-read.csv('RNA_DATA/NMOSD_RNA_DEGs(edgeR).csv',header=T)$x 
external_with_ori_DEGs<-intersect(DEG_genes,DEGs_ori)
venn.list2<-list(`Origin dataset`=DEGs_ori,`External dataset`=DEG_genes)
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Origin DEGs with External Dataset DEGs(edgeR)',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.61), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFDF80","#C0CBDF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFDF80","#C0CBDF"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 1.5,#每個圓圈裡數字大小
  cat.dist   = 0.04, #圓圈標題距離圓圈
  cat.pos    = c(-30,30), #圓圈標題角度 #往右正
  margin     = 1, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot2)


#---{ORIGIN-proteing coding genes & EXTERNAL}
#Criteria:p.adjust<0.05
library(VennDiagram)
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGs(external)(edgeR)(protein coding genes).csv',header=T)
DEG_genes<-DEGs_df$x 
DEGs_ori<-read.csv('RNA_DATA/NMOSD_RNA_DEGs(edgeR)(protein coding genes).csv',header=T)$x 
external_with_ori_DEGs<-intersect(DEG_genes,DEGs_ori)
venn.list2<-list(`Origin dataset`=DEGs_ori,`External dataset`=DEG_genes)
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Origin protein coding-DEGs with External Dataset DEGs(edgeR)',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.61), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFDF80","#C0CBDF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFDF80","#C0CBDF"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 1.5,#每個圓圈裡數字大小
  cat.dist   = 0.04, #圓圈標題距離圓圈
  cat.pos    = c(-30,30), #圓圈標題角度 #往右正
  margin     = 1, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot2)



#--------------------------------------------------------------------------------------------------------------
######################【篩"Protein coding genes" :External DEGs in Type I IFN】#######################

#---------------【Type I IFN related-DEGs (criteria:padj<0.05 & |log2FC|>1)】------------------
####################【畫2圓交集圖】################################
#----{Origin data(protein coding genes) & External data (protein coding genes)}
library(VennDiagram)
type.I.IFN.DEGs<-read.csv("./RNA_DATA/Type_I_IFN_DEG_genes(origin data)(protein codin genes)(edgeR).csv",check.names = FALSE)
IFN_ori<-type.I.IFN.DEGs$x
external_interferon_in_protein_coding_DEGs<-read.csv("./RNA_DATA/Type_I_IFN_DEG_genes(external)(protein codin genes)(edgeR).csv",check.names = FALSE)$x
venn.list2<-list(`Origin dataset`=IFN_ori,`External dataset`=external_interferon_in_protein_coding_DEGs)
ori_external_IFN_DEGs<-intersect(IFN_ori,external_interferon_in_protein_coding_DEGs)
IFN_DEG_genes_df<-bitr(ori_external_IFN_DEGs,fromType = 'SYMBOL',toType=c('ENTREZID','ENSEMBL'),OrgDb='org.Hs.eg.db')
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Origin dataset type I IFN related-protein coding DEGs with External Dataset(edgeR)',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.61), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFDF80","#C0CBDF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFDF80","#C0CBDF"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 1.8,#每個圓圈裡數字大小
  cat.dist   = 0.04, #圓圈標題距離圓圈
  cat.pos    = c(-30,30), #圓圈標題角度 #往右正
  margin     = 1, #空白區域(越大越多)
) 

# 顯示
grid.newpage()
grid.draw(venn.plot2) #save:1200 x 1200

#---------------【Hallmark Type I IFN related-DEGs (criteria:padj<0.05 & |log2FC|>1)】------------------
####################【畫2圓交集圖】#################################
#----{Origin data(protein coding gene) & External data(protein coding gene)}
library(VennDiagram)
hallmark.IFN.DEGs<-read.csv("./RNA_DATA/Hallmark_IFN_DEG_genes(origin data)(protein codin genes)(edgeR).csv",check.names = FALSE)
hallmark_IFN_ori<-hallmark.IFN.DEGs$x
external_hallmark_interferon_in_protein_coding_DEGs<-read.csv("./RNA_DATA/Hallmark_IFN_DEG_genes(external)(protein codin genes)(edgeR).csv",check.names = FALSE)$x
venn.list2<-list(`Origin dataset`=hallmark_IFN_ori,`External dataset`=external_hallmark_interferon_in_protein_coding_DEGs)
ori_external_hallmark.IFN_DEGs<-intersect(hallmark_IFN_ori,external_hallmark_interferon_in_protein_coding_DEGs)
hallmark.IFN_DEG_genes_df<-bitr(ori_external_hallmark.IFN_DEGs,fromType = 'SYMBOL',toType=c('ENTREZID','ENSEMBL'),OrgDb='org.Hs.eg.db')
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Origin dataset Hallmark Type I IFN protein coding DEGs with External Dataset(edgeR)',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.61), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFDF80","#C0CBDF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFDF80","#C0CBDF"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 2,#每個圓圈裡數字大小
  cat.dist   = 0.04, #圓圈標題距離圓圈
  cat.pos    = c(-30,30), #圓圈標題角度 #往右正
  margin     = 1, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot2)
