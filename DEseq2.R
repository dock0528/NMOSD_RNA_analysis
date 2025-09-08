setwd("C:/Users/JANE/Desktop/Wang實驗室/NMOSD研究計畫/RNA/NMOSD_RNA_analysis")
#BiocManager::install('DESeq2')
library(DESeq2)
library(edgeR)
library(ggplot2)
#========================================【My data "Deseq2" DEGs】========================================
########################{原始gene matrix}########################
Raw_count_merged<-read.csv('./RNA_DATA/Raw_count_merged_matrix.csv',row.names = 1,header=T)

#Raw count 
count_df<-as.matrix(Raw_count_merged)
rownames(count_df)<-gsub("\\.\\d+$", "",rownames(count_df)) #去除小數點以後的值

########################{Intersection gene list}########################
protein_coding_intersection_genes<-read.csv('./RNA_DATA/External_with_My_protein_coding_genes.csv')$x #19722
protein_inter_count_df<-count_df[rownames(count_df) %in% protein_coding_intersection_genes,] #19722

########################{Sample metadata}########################
#Sample & Group
sample_df=data.frame(Sample=colnames(protein_inter_count_df))
sample_df$Group <- ifelse(
  grepl("^SRR", sample_df$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
# 把 sample_df 裡的 Group 欄轉成 factor，Control 當作 baseline
sample_df$Group <- factor(sample_df$Group,
                          levels = c("Control","NMOSD"))
rownames(sample_df) <- sample_df$Sample

########################{edgeR filterByExpr}########################
dge0 <- DGEList(counts = protein_inter_count_df, group = sample_df$Group)

keep <- filterByExpr(dge0, group = sample_df$Group)   # 過濾低表達genes
dge  <- dge0[keep, , keep.lib.sizes = FALSE] #keep.lib.sizes = FALSE:重新計算每個樣本的 library size

counts_filter <- dge$counts    # 過濾後matrix

########################{DESeq2 分析}########################
# 確保樣本順序一致
counts_filter <- counts_filter[, rownames(sample_df)]

# ---- {建 DESeqDataSet} 
dds <- DESeqDataSetFromMatrix(countData = counts_filter,
                              colData   = sample_df,
                              design    = ~ Group)
# ----{執行DESeq }
dds <- DESeq(dds)


res<-results(dds)
head(results(dds,tidy=T)) #tidy=T整齊回傳
#表達量太低>>>padj=NA


#summary of differential gene expression
summary(res)

#sort summary list by padjust
res<-res[order(res$padj),]
head(res)



###########################{Volcano Plot}###########################

# res 欄位含 log2FoldChange, pvalue, padj
# 先把原本的 par 設定存起來，最後再還原
op <- par(no.readonly = TRUE) # no.readonly = TRUE:full list of parameters 

# 調整邊界，並允許在外部繪圖
par(mar = c(5, 4, 4, 6) , xpd = TRUE) #c(bottom, left, top, right)

# 基本火山圖：logFC vs. -log10(P.Value)
with(res, plot(
  log2FoldChange , -log10(pvalue),
  pch   = 20,
  main  = "NMOSD vs Healthy DEGs",
  xlim=c(-15,15),
  col   = "grey",
  xlab  = "log2 Fold Change",
  ylab  = "-log10(P Value)"
))


#----【Significants genes criteria: padj < 0.05 & |log2FoldChange|>1】-----
# 下調基因 (藍色)
with(subset(res,  padj < 0.05 & log2FoldChange   < -1),
     points(log2FoldChange, -log10(pvalue), pch=20, col='#8FA4FF')
)

# 上調基因 (紅色)
with(subset(res, padj < 0.05 & log2FoldChange   >1),
     points(log2FoldChange, -log10(pvalue), pch=20, col='#FF8080')
)

# 在圖外加圖例，往右外推 20%
legend("topright",
       inset   = c(-0.2, 0),
       legend  = c("Up", "Down"),
       title   = "Change",
       pch     = 20,
       col     = c("#FF8080", "#8FA4FF"),
       pt.cex  = 1.4, #點符號放大倍數
       bty     = "n" #the type of box 
)

# 還原原本的 par 設定
par(op)


#顯著基因
Up_DEGs<-row.names(subset(res,padj<0.05 & log2FoldChange>1)) #length(Up_DEGs):4491
Down_DEGs<-row.names(subset(res,padj<0.05 & log2FoldChange < -1)) #length(Down_DEGs):2917
DEGS<-row.names(subset(res,padj<0.05 & abs(log2FoldChange)>1)) #length(DEGS):7408
length(Up_DEGs);length(Down_DEGs);length(DEGS)
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS(my data)(Deseq2).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP(my data)(Deseq2).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN(my data)(Deseq2).csv',row.names = F)

#----【Significants genes criteria: padj < 0.01 & |log2FoldChange|>2】-----

op <- par(no.readonly = TRUE) # no.readonly = TRUE:full list of parameters 

# 調整邊界，並允許在外部繪圖
par(mar = c(5, 4, 4, 6) , xpd = TRUE) #c(bottom, left, top, right)

# 基本火山圖：logFC vs. -log10(P.Value)
with(res, plot(
  log2FoldChange , -log10(pvalue),
  pch   = 20,
  main  = "NMOSD vs Healthy DEGs",
  xlim=c(-15,15),
  col   = "grey",
  xlab  = "log2 Fold Change",
  ylab  = "-log10(P Value)"
))


# 下調基因 (藍色)
with(subset(res,  padj < 0.01 & log2FoldChange < -2),
     points(log2FoldChange, -log10(pvalue), pch=20, col='#8FA4FF')
)

# 上調基因 (紅色)
with(subset(res, padj < 0.01 & log2FoldChange >2),
     points(log2FoldChange, -log10(pvalue), pch=20, col='#FF8080')
)

# 在圖外加圖例，往右外推 20%
legend("topright",
       inset   = c(-0.2, 0),
       legend  = c("Up", "Down"),
       title   = "Change",
       pch     = 20,
       col     = c("#FF8080", "#8FA4FF"),
       pt.cex  = 1.4, #點符號放大倍數
       bty     = "n" #the type of box 
)

# 還原原本的 par 設定
par(op)


#顯著基因
Up_DEGs<-row.names(subset(res,padj<0.01 & log2FoldChange>2)) #length(Up_DEGs):1959
Down_DEGs<-row.names(subset(res,padj<0.01 & log2FoldChange < -2)) #length(Down_DEGs):1056
DEGS<-row.names(subset(res,padj<0.01 & abs(log2FoldChange)>2)) #length(DEGS):3015
length(Up_DEGs);length(Down_DEGs);length(DEGS)
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS_strictly(my data)(Deseq2).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP_strictly(my data)(Deseq2).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN_strictly(my data)(Deseq2).csv',row.names = F)
#--------------------------------------------------------------------------------------------------------------------------



#========================================【External data "Deseq2" DEGs】========================================
library(DESeq2)
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

########################{edgeR filterByExpr}########################
dge0 <- DGEList(counts = protein_inter_count_df, group = sample_df$Group)

keep <- filterByExpr(dge0, group = sample_df$Group)   # 過濾低表達genes
dge  <- dge0[keep, , keep.lib.sizes = FALSE] #keep.lib.sizes = FALSE:重新計算每個樣本的 library size

counts_filter <- dge$counts    # 過濾後matrix

########################{DESeq2 分析}########################
# 確保樣本順序一致
counts_filter <- counts_filter[, rownames(sample_df)]

# ---- {建 DESeqDataSet} 
dds <- DESeqDataSetFromMatrix(countData = counts_filter,
                              colData   = sample_df,
                              design    = ~ Group)
# ----{執行DESeq }
dds <- DESeq(dds)


res<-results(dds)
head(results(dds,tidy=T)) #tidy=T整齊回傳
#表達量太低>>>padj=NA


#summary of differential gene expression
summary(res)

#sort summary list by padjust
res<-res[order(res$padj),]
head(res)



###########################{Volcano Plot}###########################

# res 欄位含 log2FoldChange, pvalue, padj
# 先把原本的 par 設定存起來，最後再還原
op <- par(no.readonly = TRUE) # no.readonly = TRUE:full list of parameters 

# 調整邊界，並允許在外部繪圖
par(mar = c(5, 4, 4, 6) , xpd = TRUE) #c(bottom, left, top, right)

# 基本火山圖：logFC vs. -log10(P.Value)
with(res, plot(
  log2FoldChange , -log10(pvalue),
  pch   = 20,
  main  = "NMOSD vs Healthy DEGs (External Dataset)",
  xlim=c(-15,15),
  col   = "grey",
  xlab  = "log2 Fold Change",
  ylab  = "-log10(P Value)"
))


#----【Significants genes criteria: padj < 0.05 & |log2FoldChange|>1】-----
# 下調基因 (藍色)
with(subset(res,  padj < 0.05 & log2FoldChange   < -1),
     points(log2FoldChange, -log10(pvalue), pch=20, col='#8FA4FF')
)

# 上調基因 (紅色)
with(subset(res, padj < 0.05 & log2FoldChange   >1),
     points(log2FoldChange, -log10(pvalue), pch=20, col='#FF8080')
)

# 在圖外加圖例，往右外推 20%
legend("topright",
       inset   = c(-0.2, 0),
       legend  = c("Up", "Down"),
       title   = "Change",
       pch     = 20,
       col     = c("#FF8080", "#8FA4FF"),
       pt.cex  = 1.4, #點符號放大倍數
       bty     = "n" #the type of box 
)

# 還原原本的 par 設定
par(op)


#顯著基因
Up_DEGs<-row.names(subset(res,padj<0.05 & log2FoldChange>1)) #length(Up_DEGs):25
Down_DEGs<-row.names(subset(res,padj<0.05 & log2FoldChange < -1)) #length(Down_DEGs):103
DEGS<-row.names(subset(res,padj<0.05 & abs(log2FoldChange)>1)) #length(DEGS):128
length(Up_DEGs);length(Down_DEGs);length(DEGS)
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS(external data)(Deseq2).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP(external data)(Deseq2).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN(external data)(Deseq2).csv',row.names = F)

#----【Significants genes criteria: padj < 0.01 & |log2FoldChange|>2】-----

op <- par(no.readonly = TRUE) # no.readonly = TRUE:full list of parameters 

# 調整邊界，並允許在外部繪圖
par(mar = c(5, 4, 4, 6) , xpd = TRUE) #c(bottom, left, top, right)

# 基本火山圖：logFC vs. -log10(P.Value)
with(res, plot(
  log2FoldChange , -log10(pvalue),
  pch   = 20,
  main  = "NMOSD vs Healthy DEGs (External Dataset)",
  xlim=c(-15,15),
  col   = "grey",
  xlab  = "log2 Fold Change",
  ylab  = "-log10(P Value)"
))


# 下調基因 (藍色)
with(subset(res,  padj < 0.01 & log2FoldChange < -2),
     points(log2FoldChange, -log10(pvalue), pch=20, col='#8FA4FF')
)

# 上調基因 (紅色)
with(subset(res, padj < 0.01 & log2FoldChange >2),
     points(log2FoldChange, -log10(pvalue), pch=20, col='#FF8080')
)

# 在圖外加圖例，往右外推 20%
legend("topright",
       inset   = c(-0.2, 0),
       legend  = c("Up", "Down"),
       title   = "Change",
       pch     = 20,
       col     = c("#FF8080", "#8FA4FF"),
       pt.cex  = 1.4, #點符號放大倍數
       bty     = "n" #the type of box 
)

# 還原原本的 par 設定
par(op)


#顯著基因
Up_DEGs<-row.names(subset(res,padj<0.01 & log2FoldChange>2)) #length(Up_DEGs):6
Down_DEGs<-row.names(subset(res,padj<0.01 & log2FoldChange < -2)) #length(Down_DEGs):14
DEGS<-row.names(subset(res,padj<0.01 & abs(log2FoldChange)>2)) #length(DEGS):20
length(Up_DEGs);length(Down_DEGs);length(DEGS)
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS_strictly(external data)(Deseq2).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP_strictly(external data)(Deseq2).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN_strictly(external data)(Deseq2).csv',row.names = F)

#########################################【Function 計算IFN數量】#########################################
library(clusterProfiler)

#=====================================================<My Dataset>=====================================================
##################{My dataset "Deseq2" protein coding-DEGs (criteria: padjust<0.05 & |log2FC|>1}##################
#---{Interferon count}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(my data)(Deseq2).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:2886個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:50個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE)
interferon_df<-data.frame(Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
print(interferon_df)

#---{type.I.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(my data)(Deseq2).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
length(interferon_in_DEGs)
print(interferon_in_DEGs)
#write.csv(interferon_in_DEGs,file = "./RNA_DATA/Type_I_IFN_DEG_genes(my data)(Deseq2).csv",row.names = FALSE)

#---{hallmark.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(my data)(Deseq2).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
length(hallmark_interferon_in_DEGs)
print(hallmark_interferon_in_DEGs)
#write.csv(hallmark_interferon_in_DEGs,file= "./RNA_DATA/Hallmark_IFN_DEG_genes(my data)(Deseq2).csv",row.names = FALSE)


##################{My dataset "Deseq2" protein coding-DEGs (criteria: padjust<0.01 & |log2FC|>2}##################
#---{Interferon count}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(my data)(Deseq2).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:2886個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:50個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE)
interferon_df<-data.frame(Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
print(interferon_df)

#---{type.I.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(my data)(Deseq2).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
length(interferon_in_DEGs)
print(interferon_in_DEGs)
#write.csv(interferon_in_DEGs,file = "./RNA_DATA/Type_I_IFN_DEG_genes_strictly(my data)(Deseq2).csv",row.names = FALSE)

#---{hallmark.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(my data)(Deseq2).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
length(hallmark_interferon_in_DEGs)
print(hallmark_interferon_in_DEGs)
#write.csv(hallmark_interferon_in_DEGs,file= "./RNA_DATA/Hallmark_IFN_DEG_genes_strictly(my data)(Deseq2).csv",row.names = FALSE)


#=====================================================<External Dataset>=====================================================
##################{External dataset "Deseq2" protein coding-DEGs (criteria: padjust<0.05 & |log2FC|>1)}##################
#---{Interferon count}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external data)(Deseq2).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:2886個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:50個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE)
interferon_df<-data.frame(Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
print(interferon_df)

#---{type.I.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external data)(Deseq2).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
length(interferon_in_DEGs)
print(interferon_in_DEGs)
#write.csv(interferon_in_DEGs,file = "./RNA_DATA/Type_I_IFN_DEG_genes(external data)(Deseq2).csv",row.names = FALSE)

#---{hallmark.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external data)(Deseq2).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
length(hallmark_interferon_in_DEGs)
print(hallmark_interferon_in_DEGs)
#write.csv(hallmark_interferon_in_DEGs,file= "./RNA_DATA/Hallmark_IFN_DEG_genes(external data)(Deseq2).csv",row.names = FALSE)


##################{External dataset "Deseq2" protein coding-DEGs (criteria: padjust<0.01 & |log2FC|>2}##################
#---{Interferon count}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(external data)(Deseq2).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:2886個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:50個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE)
interferon_df<-data.frame(Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
print(interferon_df)

#---{type.I.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(external data)(Deseq2).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
length(interferon_in_DEGs)
print(interferon_in_DEGs)
#write.csv(interferon_in_DEGs,file = "./RNA_DATA/Type_I_IFN_DEG_genes_strictly(external data)(Deseq2).csv",row.names = FALSE)

#---{hallmark.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(external data)(Deseq2).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
length(hallmark_interferon_in_DEGs)
print(hallmark_interferon_in_DEGs)
#write.csv(hallmark_interferon_in_DEGs,file= "./RNA_DATA/Hallmark_IFN_DEG_genes_strictly(external data)(Deseq2).csv",row.names = FALSE)












###【Enrichment analysis】
BiocManager::install('DOSE')
BiocManager::install('org.Hs.eg.db')
BiocManager::install('topGO')
BiocManager::install('clusterProfiler')
BiocManager::install('pathview')

library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

###【GO】
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS.csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')

###DEG ENTREZID資料集
DEG_genes<-DEG_genes_df$ENTREZID

#--------------【GO(MF)】------------------
#進行 enrichGO 富集分析
ego_MF<-enrichGO(gene=DEG_genes,
                 OrgDb=org.Hs.eg.db,
                 keyType = 'ENTREZID',
                 ont='MF', #也可以是CC、BP
                 pAdjustMethod = 'BH', #Benjamini-Hochberg
                 pvalueCutoff = 1, #1為不過濾
                 qvalueCutoff = 1, #1為不過濾
                 readable = T) #將Gene ID轉成gene Symbol(易讀)
ego_MF_df<-ego_MF@result
View(ego_MF_df)

#dotplot
dotplot(ego_MF,title='EnrichmentGO_MF_dot',label_format = 100)
#label_format讓文字不會重疊

#barplot
barplot(ego_MF,showCategory=20,title='EnrichmentGO_MF_bar',label_format = 100)
#showCategory=20繪製前20個

#過濾
ego_fil<-simplify(ego_MF,cutoff=0.5,by='pvalue',select_fun = match.fun("min"))
View(ego_fil@result)

#---------------------------------------------

#--------------【GO(MF)】------------------
#進行 enrichGO 富集分析
ego_BP<-enrichGO(gene=DEG_genes,
                 OrgDb=org.Hs.eg.db,
                 keyType = 'ENTREZID',
                 ont='BP', #也可以是CC、BP
                 pAdjustMethod = 'BH', #Benjamini-Hochberg
                 pvalueCutoff = 1, #1為不過濾
                 qvalueCutoff = 1, #1為不過濾
                 readable = T) #將Gene ID轉成gene Symbol(易讀)
ego_BP_df<-ego_BP@result
View(ego_BP_df)

#dotplot
dotplot(ego_BP,showCategory=10,title='EnrichmentGO_BP',label_format = 100)
#label_format讓文字不會重疊

#barplot
barplot(ego_BP,showCategory=20,title='EnrichmentGO_BP_bar',label_format = 100)
#showCategory=20繪製前20個

#---【與 type I interferon 有關的 pathway】
typeI_IFN <- grepl("type I interferon", ego_BP@result$Description, ignore.case = TRUE)

# 建立新的資料框
typeI_IFN_df <- ego_BP@result[typeI_IFN, ]

# 檢視結果（這才是資料表）
View(typeI_IFN_df)
#--------------------------------------------
#KEGG pathway富集分析
kk<-enrichKEGG(gene=DEG_genes,#你的基因列表
               organism='hsa', #指定物種 #hsa人類
               pvalueCutoff = 1) #1不進行過濾
kk@result


#KEGG dotplot
dotplot(kk,title='Enrichment KEGG_dot')

###(特定)KEGG pahtway
genelist<-as.numeric(res[1:4289,'log2FoldChange']) #DEGs_gene有4289個
names(genelist)<-DEG_genes_df[1:4289,'ENTREZID'] #pathview默認是entrez id
select_pathway<-kk@result$ID[1] #hsa04610

pathview(gene.data = genelist,
         pathway.id = select_pathway,
         species='hsa',
         kegg.native = T, #是否使用KEGG原生圖像
         new.sigmature=F, #是否生成新的標誌(T=>生成新的標誌以顯示顯著變化的基因)
         limit=list(gene=2.5,cpd=1))#基因表達的限制是 2.5，對化合物的限制是 1

