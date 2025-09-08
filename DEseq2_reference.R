setwd("C:/Users/JANE/Desktop/Wang實驗室/NMOSD研究計畫/RNA/NMOSD_RNA_analysis")
#BiocManager::install('DESeq2')
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

# ---- {建 DESeqDataSet} 
dds <- DESeqDataSetFromMatrix(countData = protein_inter_count_df,
                              colData  = sample_df,
                              design   = ~ Group)
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
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS(external)(Deseq2).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP(external)(Deseq2).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN(external)(Deseq2).csv',row.names = F)

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
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS_strictly(external)(Deseq2).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP_strictly(external)(Deseq2).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN_strictly(external)(Deseq2).csv',row.names = F)


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

# ---- {建 DESeqDataSet} 
#先進行過濾，讓後續結果的log2FC不會過於極大極小 & 不易p-value不穩定或padj=NA
keep <- rowSums(protein_coding_count_df >= 10) >= 5   # 至少在 ≥5 個樣本 count >= 10（可依資料量調整）
#單個樣本有 10 個 reads，對應的 CPM / TPM 大概落在一個可檢測且穩定的表達量區間
protein_coding_count_df <- protein_coding_count_df[keep, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = protein_coding_count_df,
                              colData  = sample_df,
                              design   = ~ Group)
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
Up_DEGs<-row.names(subset(res,padj<0.05 & log2FoldChange>1)) #length(Up_DEGs):24
Down_DEGs<-row.names(subset(res,padj<0.05 & log2FoldChange < -1)) #length(Down_DEGs):82
DEGS<-row.names(subset(res,padj<0.05 & abs(log2FoldChange)>1)) #length(DEGS):106
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS(external)(Deseq2)(protein coding genes).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP(external)(Deseq2)(protein coding genes).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN(external)(Deseq2)(protein coding genes).csv',row.names = F)

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
Down_DEGs<-row.names(subset(res,padj<0.01 & log2FoldChange < -2)) #length(Down_DEGs):11
DEGS<-row.names(subset(res,padj<0.01 & abs(log2FoldChange)>2)) #length(DEGS):17
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS_strictly(external)(Deseq2)(protein coding genes).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP_strictly(external)(Deseq2)(protein coding genes).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN_strictly(external)(Deseq2)(protein coding genes).csv',row.names = F)
#-------------------------------------------------------------------------------------------------------------------------


#########################################【Function 計算IFN數量】#########################################
library(clusterProfiler)
#---------{Origin "Deseq2" DEGs (criteria: padjust<0.05 & |log2FC|>1}---------

#---{Interferon count}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external)(Deseq2).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:2886個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:50個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE)
interferon_df<-data.frame(Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
print(interferon_df)

#---{type.I.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external)(Deseq2).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
print(interferon_in_DEGs)
#write.csv(interferon_in_DEGs,file = "./RNA_DATA/Type_I_IFN_DEG_genes(external)(Deseq2).csv",row.names = FALSE)

#---{hallmark.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external)(Deseq2).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
print(hallmark_interferon_in_DEGs)
#write.csv(hallmark_interferon_in_DEGs,file= "./RNA_DATA/Hallmark_IFN_DEG_genes(external)(Deseq2).csv",row.names = FALSE)

#---------{Origin protein coding-DEGs (criteria: padjust<0.05 & |log2FC|>1}---------

#---{Interferon count}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external)(Deseq2)(protein coding genes).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE)
interferon_df<-data.frame(Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
print(interferon_df)

#---{type.I.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external)(Deseq2)(protein coding genes).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
print(interferon_in_DEGs) 
length(interferon_in_DEGs)
#write.csv(interferon_in_DEGs,file = "./RNA_DATA/Type_I_IFN_DEG_genes(external)(protein codin genes)(Deseq2).csv",row.names = FALSE)

#---{hallmark.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external)(Deseq2)(protein coding genes).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
print(hallmark_interferon_in_DEGs)
length(hallmark_interferon_in_DEGs)
#write.csv(hallmark_interferon_in_DEGs,file= "./RNA_DATA/Hallmark_IFN_DEG_genes(external)(protein codin genes)(Deseq2).csv",row.names = FALSE)
#---------------------------------------------------------------------------------------------------------------------


#############################################【Deseq2 - DEGs interection】#############################################
#---{ORIGIN & EXTERNAL}
library(VennDiagram)
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external)(Deseq2).csv',header=T)
DEG_genes<-DEGs_df$x 
DEGs_ori<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(Deseq2).csv',header=T)$x 
external_with_ori_DEGs<-intersect(DEG_genes,DEGs_ori)
venn.list2<-list(`Origin dataset`=DEGs_ori,`External dataset`=DEG_genes)
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Origin DEGs with External Dataset DEGs(Deseq2)',
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
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external)(Deseq2)(protein coding genes).csv',header=T)
DEG_genes<-DEGs_df$x 
DEGs_ori<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(Deseq2)(protein coding genes).csv',header=T)$x 
external_with_ori_DEGs<-intersect(DEG_genes,DEGs_ori)
venn.list2<-list(`Origin dataset`=DEGs_ori,`External dataset`=DEG_genes)
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Origin protein coding-DEGs with External Dataset DEGs(Deseq2)',
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

#Criteria:p.adjust<0.01
library(VennDiagram)
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(External dataset)(protein coding genes).csv',header=T)
DEG_genes<-DEGs_df$x #16個
DEGs_ori<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_origin_strictly(protein codin genes).csv',header=T)$x #DEGs count1177
external_with_ori_DEGs<-intersect(DEG_genes,DEGs_ori)
venn.list2<-list(`Origin dataset`=DEGs_ori,`External dataset`=DEG_genes)
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Origin protein coding-DEGs with External Dataset DEGs',
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
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
###【GO】
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external)(Deseq2)(protein coding genes).csv',header=T)

DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
#轉成ENTREZID & SYMBOL (從37DEGs -> 37DEGs)

#消失的ENSG (70個)
DEG_genes[!DEG_genes %in% DEG_genes_df$ENSEMBL]

#重複的ENSG(3個ENSG)
dup_ENSGs<-DEG_genes_df[duplicated(DEG_genes_df$ENSEMBL),]$ENSEMBL

#重複的ENTREZID(1個)
dup_ENTREZIDs<-DEG_genes_df[duplicated(DEG_genes_df$ENTREZID),]$ENTREZID
DEG_genes_df[DEG_genes_df$ENTREZID %in% dup_ENTREZIDs,]

#重複的SYMBOL(1個)
dup_SYMBOLs<-DEG_genes_df[duplicated(DEG_genes_df$SYMBOL),]$SYMBOL
DEG_genes_df[DEG_genes_df$SYMBOL %in% dup_SYMBOLs,]

###DEG ENTREZID資料集
DEG_genes<-DEG_genes_df$ENTREZID #300個DEG

#取和MsigDB有交集的type l interferon genes
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes) #和DEGs交集有102
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有20

external_interferon_in_protein_coding_DEGs<-c(
  "IFI27", "EIF1AY", "USP18", "LY6E", "NETO1", "HLA-DQA1", "LARGE1", "TNFRSF13C", "ITGB8", "CPA3", "OAS1",
  "AKAP12", "IFI44L", "IL4", "OAS2", "ISG15", "IFI44", "RSAD2", "GATM", "KDM5D", "OASL", "HRH4",
  "CA8", "IFIT1", "COBLL1", "XKRX", "PTPRS", "IFI6", "ADAM12", "DDX60", "XAF1", "UTY", "FCRL2",
  "RTP4", "TNFRSF13B", "MX1", "LAMC1", "WNT16", "FGFR2", "DDX3Y", "POU2AF1", "SLC7A8", "IFIT2", "CMPK2",
  "HPGDS", "FCRLA", "RAPGEF3", "OAS3", "THSD7A", "CD79A", "UGT2B17", "BHLHE41", "IL17RB", "HERC5", "IFIT3",
  "SPATS2L"
)#56

external_hallmark_interferon_in_protein_coding_DEGs<-c(
  "IFI27", "USP18", "LY6E", "OAS1", "IFI44L", "ISG15", "IFI44", "RSAD2", 
  "OASL", "DDX60", "RTP4", "MX1", "IFIT2", "CMPK2", "IFIT3"
)#15

#---------------【Type I IFN related-DEGs (criteria:padj<0.05 & |log2FC|>1)】------------------
####################【畫2圓交集圖】################################
#----{Origin data(protein coding genes) & External data (protein coding genes)}
library(VennDiagram)
type.I.IFN.DEGs<-read.csv("./RNA_DATA/Type_I_IFN_DEG_genes(origin data)(protein codin genes)(Deseq2).csv",check.names = FALSE)
IFN_ori<-type.I.IFN.DEGs$x
external_interferon_in_protein_coding_DEGs<-c(
  "IFI27", "EIF1AY", "USP18", "LY6E", "NETO1", "HLA-DQA1", "LARGE1", "TNFRSF13C", "ITGB8", "CPA3", "OAS1",
  "AKAP12", "IFI44L", "IL4", "OAS2", "ISG15", "IFI44", "RSAD2", "GATM", "KDM5D", "OASL", "HRH4",
  "CA8", "IFIT1", "COBLL1", "XKRX", "PTPRS", "IFI6", "ADAM12", "DDX60", "XAF1", "UTY", "FCRL2",
  "RTP4", "TNFRSF13B", "MX1", "LAMC1", "WNT16", "FGFR2", "DDX3Y", "POU2AF1", "SLC7A8", "IFIT2", "CMPK2",
  "HPGDS", "FCRLA", "RAPGEF3", "OAS3", "THSD7A", "CD79A", "UGT2B17", "BHLHE41", "IL17RB", "HERC5", "IFIT3",
  "SPATS2L"
)#56
venn.list2<-list(`Origin dataset`=IFN_ori,`External dataset`=external_interferon_in_protein_coding_DEGs)
ori_external_IFN_DEGs<-intersect(IFN_ori,external_interferon_in_protein_coding_DEGs)
IFN_DEG_genes_df<-bitr(ori_external_IFN_DEGs,fromType = 'SYMBOL',toType=c('ENTREZID','ENSEMBL'),OrgDb='org.Hs.eg.db')
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Origin dataset type I IFN related-protein coding DEGs with External Dataset(Deseq2)',
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
hallmark.IFN.DEGs<-read.csv("./RNA_DATA/Hallmark_IFN_DEG_genes(origin data)(protein codin genes)(Deseq2).csv",check.names = FALSE)
hallmark_IFN_ori<-hallmark.IFN.DEGs$x
external_hallmark_interferon_in_protein_coding_DEGs<-c(
  "IFI27", "USP18", "LY6E", "OAS1", "IFI44L", "ISG15", "IFI44", "RSAD2", 
  "OASL", "DDX60", "RTP4", "MX1", "IFIT2", "CMPK2", "IFIT3"
)#15
venn.list2<-list(`Origin dataset`=hallmark_IFN_ori,`External dataset`=external_hallmark_interferon_in_protein_coding_DEGs)
ori_external_hallmark.IFN_DEGs<-intersect(hallmark_IFN_ori,external_hallmark_interferon_in_protein_coding_DEGs)
hallmark.IFN_DEG_genes_df<-bitr(ori_external_hallmark.IFN_DEGs,fromType = 'SYMBOL',toType=c('ENTREZID','ENSEMBL'),OrgDb='org.Hs.eg.db')
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Origin dataset Hallmark Type I IFN protein coding DEGs with External Dataset(Deseq2)',
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
