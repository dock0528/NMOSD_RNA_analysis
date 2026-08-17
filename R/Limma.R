setwd("C:/Users/JANE/Desktop/Wang實驗室/NMOSD研究計畫/RNA/NMOSD_RNA_analysis")

#========================================【My data "Limma+Voom" DEGs】========================================
#載入套件
#BiocManager::install(c("edgeR","limma"))
library(limma)
library(edgeR)

########################{原始gene matrix}########################
Raw_count_merged<-read.csv('./RNA_DATA/Raw_count_merged_matrix.csv',row.names = 1,header=T) #dim(Raw_count_merged):78724 x 21

#Raw count 
count_df<-as.matrix(Raw_count_merged)
rownames(count_df)<-gsub("\\.\\d+$", "",rownames(count_df)) #去除小數點以後的值

########################{Intersection gene list}########################
intersection_genes<-row.names(read.csv('./RNA_DATA/External_TPM_gene_samples_matrix(intersection gene).csv',row.names = 1,header=TRUE))
#56657 intersection genes (My & External)

protein_coding_intersection_genes<-read.csv('./RNA_DATA/External_with_My_protein_coding_genes.csv')$x #19722

inter_count_df<-count_df[rownames(count_df) %in% intersection_genes,] #56657
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

#----【Create DGEList】----
dge <- DGEList(counts=protein_inter_count_df, group=sample_df$Group) #counts:row->genes,col->samples


# ----【過濾低表達genes】 ----
# 依組別與j文庫大小自動估門檻（約等 CPM>=1 在足夠多樣本）
keep <- filterByExpr(dge, group = sample_df$Group)
dge  <- dge[keep, , keep.lib.sizes = FALSE]

#----【TMM（Trimmed Mean of M-values）】----
#目的:校正不同樣本間因測序深度或組成偏差所帶來的系統性差異
#Step1.消除極端表達基因的影響
#Step2.估算每個樣本的「有效 library size」
#Step3.使樣本間的基因表達量可直接比較
dge <- calcNormFactors(dge,method = 'TMM')    # TMM 標準化


########################{Voom兩步驟}########################
#轉成log2CPM -> 讓基因可以在不同cell做比較
#將低變異度的genes變異度變小(RNA-seq:低豐度基因的變異通常比高豐度基因大)

#----【Voom轉換】----
#建立 design matrix（不含截距，以 Control 作 baseline）
design <- model.matrix(~ 0 + sample_df$Group)
colnames(design) <- levels(sample_df$Group)  #Control vs NMOSD


#----【將 DGEList 轉成 voom 物件（log2-CPM + 權重）】----
v <- voom(dge, design, plot = TRUE)
# plot = TRUE 會畫出 mean-variance trend，可用來檢查變異度是否隨平均表達量下降

#---【取log2cpm Matrix】
log2cpm_matrix  <- v$E                           # voom 轉出的 log2-CPM matrix 
#write.csv(log2cpm_matrix,file='RNA_DATA/log2cpm_matrix(My dataset)(19722genes).csv')

#############{Limma三步驟:lmFit + contrasts.fit + eBayes}############
#----【線性模型fitting】----
# fit <- lmFit(log2-CPM 與 權重)
fit <- lmFit(v, design)

#----【設定對比（contrast）】----
# NMOSD vs Control
contrast.matrix <- makeContrasts(
  NMOSD_vs_Control = NMOSD - Control,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)

#----【Empirical Bayes 調整】----
fit2 <- eBayes(fit2)

######################################################################

# 產生全基因結果（不排序）
Allgene <- topTable(fit2,
                   coef         = "NMOSD_vs_Control",
                   adjust.method= "BH",
                   number       = Inf,        # 回傳所有基因
                   sort.by      = "none")     # 不排序，保留原始順序

###########################{Volcano Plot}###########################

# Allgene 欄位含 logFC, P.Value, adj.P.Val
# 先把原本的 par 設定存起來，最後再還原
op <- par(no.readonly = TRUE) # no.readonly = TRUE:full list of parameters 

# 調整邊界，並允許在外部繪圖
par(mar = c(5, 4, 4, 6) , xpd = TRUE) #c(bottom, left, top, right)

# 基本火山圖：logFC vs. -log10(P.Value)
with(Allgene, plot(
  logFC, -log10(P.Value),
  pch   = 20,
  main  = "NMOSD vs Healthy DEGs",
  xlim=c(-15,15),
  col   = "grey",
  xlab  = "log2 Fold Change",
  ylab  = "-log10(P Value)"
))


#----【Significants genes criteria: adj.P.Val < 0.05 & |logFC|>1】-----
# 下調基因 (藍色)
with(subset(Allgene, adj.P.Val < 0.05 & logFC < -1),
     points(logFC, -log10(P.Value), pch=20, col='#8FA4FF')
)

# 上調基因 (紅色)
with(subset(Allgene, adj.P.Val < 0.05 & logFC >1),
     points(logFC, -log10(P.Value), pch=20, col='#FF8080')
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
Up_DEGs<-row.names(subset(Allgene,adj.P.Val < 0.05 & logFC > 1)) 
Down_DEGs<-row.names(subset(Allgene,adj.P.Val < 0.05 & logFC < -1)) 
DEGS<-row.names(subset(Allgene,adj.P.Val<0.05 & abs(logFC)>1)) 
length(Up_DEGs);length(Down_DEGs);length(DEGS)
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS(my data)(limma).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP(my data)(limma).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN(my data)(limma).csv',row.names = F)

#----【Significants genes criteria: adj.P.Val < 0.01 & |logFC|>2】----
# Allgene 欄位含 logFC, P.Value, adj.P.Val
# 先把原本的 par 設定存起來，最後再還原
op <- par(no.readonly = TRUE) # no.readonly = TRUE:full list of parameters 

# 調整邊界，並允許在外部繪圖
par(mar = c(5, 4, 4, 6) , xpd = TRUE) #c(bottom, left, top, right)

# 基本火山圖：logFC vs. -log10(P.Value)
with(Allgene, plot(
  logFC, -log10(P.Value),
  pch   = 20,
  main  = "NMOSD vs Healthy DEGs",
  xlim=c(-15,15),
  col   = "grey",
  xlab  = "log2 Fold Change",
  ylab  = "-log10(P Value)"
))
# 下調基因 (藍色)
with(subset(Allgene, adj.P.Val < 0.01 & logFC < -2),
     points(logFC, -log10(P.Value), pch=20, col='#8FA4FF')
)

# 上調基因 (紅色)
with(subset(Allgene, adj.P.Val < 0.01 & logFC >2),
     points(logFC, -log10(P.Value), pch=20, col='#FF8080')
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
Up_DEGs<-row.names(subset(Allgene,adj.P.Val < 0.01 & logFC > 2)) 
Down_DEGs<-row.names(subset(Allgene,adj.P.Val < 0.01 & logFC < -2)) 
DEGS<-row.names(subset(Allgene,adj.P.Val<0.01 & abs(logFC)>2)) 
length(Up_DEGs);length(Down_DEGs);length(DEGS)
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS_strictly(my data)(limma).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP_strictly(my data)(limma).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN_strictly(my data)(limma).csv',row.names = F)


#################################【示範連到Ensembl抓protein coding gene】#################################
Raw_count_merged<-read.csv('./RNA_DATA/Raw_count_merged_matrix.csv',row.names = 1,header=T) #dim(Raw_count_merged):78724 x 21

#Raw count 
count_df<-as.matrix(Raw_count_merged)
rownames(count_df)<-gsub("\\.\\d+$", "",rownames(count_df)) #去除小數點以後的值

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

#-----------------------------------------------------------------------------------------------------------------------------------------------

#========================================【External data "Limma+Voom" DEGs】========================================
library(limma)
library(edgeR)

########################{原始gene matrix}########################
Raw_count_merged<-read.csv('./RNA_DATA/External_RNA_datasets_DEG_RESULT/RNA_gene_matrix_v2.csv',row.names = 1,header=T) #dim(Raw_count_merged):58304 x 23
Raw_count_merged <- Raw_count_merged[!rownames(Raw_count_merged) %in% c("Antibody status","Ensembl_Gene_ID"), ]
#dim(Raw_count_merged) 58302 x 23

#Raw count 
count_df<-as.matrix(Raw_count_merged)
storage.mode(count_df) <- "numeric"

########################{Intersection gene list}########################
intersection_genes<-row.names(read.csv('./RNA_DATA/External_TPM_gene_samples_matrix(intersection gene).csv',row.names = 1,header=TRUE))
#56657 intersection genes (My & External)

protein_coding_intersection_genes<-read.csv('./RNA_DATA/External_with_My_protein_coding_genes.csv')$x #19722

inter_count_df<-count_df[rownames(count_df) %in% intersection_genes,] #56657
protein_inter_count_df<-count_df[rownames(count_df) %in% protein_coding_intersection_genes,] #19722

########################{Sample metadata}########################
#Sample & Group
sample_df=data.frame(Sample=colnames(protein_inter_count_df))
sample_df$Group <- ifelse(
  grepl("^NA", sample_df$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
# 把 sample_df 裡的 Group 欄轉成 factor，Control 當作 baseline
sample_df$Group <- factor(sample_df$Group,
                          levels = c("Control","NMOSD"))

#----【Create DGEList】----
dge <- DGEList(counts=protein_inter_count_df, group=sample_df$Group) #counts:row->genes,col->samples

# ----【過濾低表達genes】 ----
# 依組別與j文庫大小自動估門檻（約等 CPM>=1 在足夠多樣本）
keep <- filterByExpr(dge, group = sample_df$Group)
dge  <- dge[keep, , keep.lib.sizes = FALSE]

#----【TMM（Trimmed Mean of M-values）】----
#目的:校正不同樣本間因測序深度或組成偏差所帶來的系統性差異
#Step1.消除極端表達基因的影響
#Step2.估算每個樣本的「有效 library size」
#Step3.使樣本間的基因表達量可直接比較
dge <- calcNormFactors(dge,method = 'TMM')    # TMM 標準化


########################{Voom兩步驟}########################
#轉成log2CPM -> 讓基因可以在不同cell做比較
#將低變異度的genes變異度變小(RNA-seq:低豐度基因的變異通常比高豐度基因大)

#----【Voom轉換】----
#建立 design matrix（不含截距，以 Control 作 baseline）
design <- model.matrix(~ 0 + sample_df$Group)
colnames(design) <- levels(sample_df$Group)  #Control vs NMOSD


#----【將 DGEList 轉成 voom 物件（log2-CPM + 權重）】----
v <- voom(dge, design, plot = TRUE)
# plot = TRUE 會畫出 mean-variance trend，可用來檢查變異度是否隨平均表達量下降

#---【取log2cpm Matrix】
log2cpm_matrix  <- v$E                           # voom 轉出的 log2-CPM matrix 
#write.csv(log2cpm_matrix,file='RNA_DATA/log2cpm_matrix(External dataset)(19722genes).csv')

#############{Limma三步驟:lmFit + contrasts.fit + eBayes}############
#----【線性模型fitting】----
# fit <- lmFit(log2-CPM 與 權重)
fit <- lmFit(v, design)

#----【設定對比（contrast）】----
# NMOSD vs Control
contrast.matrix <- makeContrasts(
  NMOSD_vs_Control = NMOSD - Control,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)

#----【Empirical Bayes 調整】----
fit2 <- eBayes(fit2)

######################################################################

# 產生全基因結果（不排序）
Allgene <- topTable(fit2,
                    coef         = "NMOSD_vs_Control",
                    adjust.method= "BH",
                    number       = Inf,        # 回傳所有基因
                    sort.by      = "P")     # 不排序，保留原始順序

###########################{Volcano Plot}###########################

# Allgene 欄位含 logFC, P.Value, adj.P.Val
# 先把原本的 par 設定存起來，最後再還原
op <- par(no.readonly = TRUE) # no.readonly = TRUE:full list of parameters 

# 調整邊界，並允許在外部繪圖
par(mar = c(5, 4, 4, 6) , xpd = TRUE) #c(bottom, left, top, right)

# 基本火山圖：logFC vs. -log10(P.Value)
with(Allgene, plot(
  logFC, -log10(P.Value),
  pch   = 20,
  main  = "NMOSD vs Healthy DEGs (External data)",
  xlim=c(-15,15),
  col   = "grey",
  xlab  = "log2 Fold Change",
  ylab  = "-log10(P Value)"
))


#----【Significants genes criteria: adj.P.Val < 0.05 & |logFC|>1】-----
# 下調基因 (藍色)
with(subset(Allgene, adj.P.Val < 0.05 & logFC < -1),
     points(logFC, -log10(P.Value), pch=20, col='#8FA4FF')
)

# 上調基因 (紅色)
with(subset(Allgene, adj.P.Val < 0.05 & logFC >1),
     points(logFC, -log10(P.Value), pch=20, col='#FF8080')
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
Up_DEGs<-row.names(subset(Allgene,adj.P.Val < 0.05 & logFC > 1)) 
Down_DEGs<-row.names(subset(Allgene,adj.P.Val < 0.05 & logFC < -1)) 
DEGS<-row.names(subset(Allgene,adj.P.Val<0.05 & abs(logFC)>1)) 
length(Up_DEGs);length(Down_DEGs);length(DEGS)
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS(external data)(limma).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP(external data)(limma).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN(external data)(limma).csv',row.names = F)

#----【Significants genes criteria: adj.P.Val < 0.01 & |logFC|>2】----
# Allgene 欄位含 logFC, P.Value, adj.P.Val
# 先把原本的 par 設定存起來，最後再還原
op <- par(no.readonly = TRUE) # no.readonly = TRUE:full list of parameters 

# 調整邊界，並允許在外部繪圖
par(mar = c(5, 4, 4, 6) , xpd = TRUE) #c(bottom, left, top, right)

# 基本火山圖：logFC vs. -log10(P.Value)
with(Allgene, plot(
  logFC, -log10(P.Value),
  pch   = 20,
  main  = "NMOSD vs Healthy DEGs (External data)",
  xlim=c(-15,15),
  col   = "grey",
  xlab  = "log2 Fold Change",
  ylab  = "-log10(P Value)"
))
# 下調基因 (藍色)
with(subset(Allgene, adj.P.Val < 0.01 & logFC < -2),
     points(logFC, -log10(P.Value), pch=20, col='#8FA4FF')
)

# 上調基因 (紅色)
with(subset(Allgene, adj.P.Val < 0.01 & logFC >2),
     points(logFC, -log10(P.Value), pch=20, col='#FF8080')
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
Up_DEGs<-row.names(subset(Allgene,adj.P.Val < 0.01 & logFC > 2)) 
Down_DEGs<-row.names(subset(Allgene,adj.P.Val < 0.01 & logFC < -2)) 
DEGS<-row.names(subset(Allgene,adj.P.Val<0.01 & abs(logFC)>2)) 
length(Up_DEGs);length(Down_DEGs);length(DEGS)
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS_strictly(external data)(limma).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP_strictly(external data)(limma).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN_strictly(external data)(limma).csv',row.names = F)

#########################################【Function 計算IFN數量】#########################################
library(clusterProfiler)

#=====================================================<My Dataset>=====================================================
##################{My dataset "Limma+Voom" protein coding-DEGs (criteria: padjust<0.05 & |log2FC|>1}##################
#---{Interferon count}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(my data)(limma).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:2886個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:50個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE)
interferon_df<-data.frame(Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
print(interferon_df)

#---{type.I.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(my data)(limma).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
length(interferon_in_DEGs)
print(interferon_in_DEGs)
#write.csv(interferon_in_DEGs,file = "./RNA_DATA/Type_I_IFN_DEG_genes(my data)(limma).csv",row.names = FALSE)

#---{hallmark.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(my data)(limma).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
length(hallmark_interferon_in_DEGs)
print(hallmark_interferon_in_DEGs)
#write.csv(hallmark_interferon_in_DEGs,file= "./RNA_DATA/Hallmark_IFN_DEG_genes(my data)(limma).csv",row.names = FALSE)


##################{My dataset "Limma+Voom" protein coding-DEGs (criteria: padjust<0.01 & |log2FC|>2}##################
#---{Interferon count}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(my data)(limma).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:2886個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:50個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE)
interferon_df<-data.frame(Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
print(interferon_df)

#---{type.I.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(my data)(limma).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
length(interferon_in_DEGs)
print(interferon_in_DEGs)
#write.csv(interferon_in_DEGs,file = "./RNA_DATA/Type_I_IFN_DEG_genes_strictly(my data)(limma).csv",row.names = FALSE)

#---{hallmark.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(my data)(limma).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
length(hallmark_interferon_in_DEGs)
print(hallmark_interferon_in_DEGs)
#write.csv(hallmark_interferon_in_DEGs,file= "./RNA_DATA/Hallmark_IFN_DEG_genes_strictly(my data)(limma).csv",row.names = FALSE)

#=====================================================<External Dataset>=====================================================
##################{My dataset "Limma+Voom" protein coding-DEGs (criteria: padjust<0.05 & |log2FC|>1}##################
#---{Interferon count}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external data)(limma).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:2886個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:50個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE)
interferon_df<-data.frame(Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
print(interferon_df)

#---{type.I.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external data)(limma).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
length(interferon_in_DEGs)
print(interferon_in_DEGs)
#write.csv(interferon_in_DEGs,file = "./RNA_DATA/Type_I_IFN_DEG_genes(external data)(limma).csv",row.names = FALSE)

#---{hallmark.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(external data)(limma).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
length(hallmark_interferon_in_DEGs)
print(hallmark_interferon_in_DEGs)
#write.csv(hallmark_interferon_in_DEGs,file= "./RNA_DATA/Hallmark_IFN_DEG_genes(external data)(limma).csv",row.names = FALSE)



#========================================【TEST :External dataset take out NMO#0020 & VIMS-119】========================================
library(limma)
library(edgeR)

########################{原始gene matrix}########################
Raw_count_merged<-read.csv('./RNA_DATA/External_RNA_datasets_DEG_RESULT/RNA_gene_matrix_v3.csv',row.names = 1,header=T) #dim(Raw_count_merged):58304 x 23
Raw_count_merged <- Raw_count_merged[!rownames(Raw_count_merged) %in% c("Antibody status","Ensembl_Gene_ID"), ]
#dim(Raw_count_merged) 58302 x 21

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
  grepl("^NA", sample_df$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
# 把 sample_df 裡的 Group 欄轉成 factor，Control 當作 baseline
sample_df$Group <- factor(sample_df$Group,
                          levels = c("Control","NMOSD"))

#----【Create DGEList】----
dge <- DGEList(counts=protein_inter_count_df, group=sample_df$Group) #counts:row->genes,col->samples

# ----【過濾低表達genes】 ----
# 依組別與j文庫大小自動估門檻（約等 CPM>=1 在足夠多樣本）
keep <- filterByExpr(dge, group = sample_df$Group)
dge  <- dge[keep, , keep.lib.sizes = FALSE]

#----【TMM（Trimmed Mean of M-values）】----
#目的:校正不同樣本間因測序深度或組成偏差所帶來的系統性差異
#Step1.消除極端表達基因的影響
#Step2.估算每個樣本的「有效 library size」
#Step3.使樣本間的基因表達量可直接比較
dge <- calcNormFactors(dge,method = 'TMM')    # TMM 標準化


########################{Voom兩步驟}########################
#轉成log2CPM -> 讓基因可以在不同cell做比較
#將低變異度的genes變異度變小(RNA-seq:低豐度基因的變異通常比高豐度基因大)

#----【Voom轉換】----
#建立 design matrix（不含截距，以 Control 作 baseline）
design <- model.matrix(~ 0 + sample_df$Group)
colnames(design) <- levels(sample_df$Group)  #Control vs NMOSD


#----【將 DGEList 轉成 voom 物件（log2-CPM + 權重）】----
v <- voom(dge, design, plot = TRUE)
# plot = TRUE 會畫出 mean-variance trend，可用來檢查變異度是否隨平均表達量下降

#---【取log2cpm Matrix】
log2cpm_matrix  <- v$E                           # voom 轉出的 log2-CPM matrix 
#write.csv(log2cpm_matrix,file='RNA_DATA/log2cpm_matrix(External dataset)(19722genes).csv')

#############{Limma三步驟:lmFit + contrasts.fit + eBayes}############
#----【線性模型fitting】----
# fit <- lmFit(log2-CPM 與 權重)
fit <- lmFit(v, design)

#----【設定對比（contrast）】----
# NMOSD vs Control
contrast.matrix <- makeContrasts(
  NMOSD_vs_Control = NMOSD - Control,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)

#----【Empirical Bayes 調整】----
fit2 <- eBayes(fit2)

######################################################################

# 產生全基因結果（不排序）
Allgene <- topTable(fit2,
                    coef         = "NMOSD_vs_Control",
                    adjust.method= "BH",
                    number       = Inf,        # 回傳所有基因
                    sort.by      = "P")     # 不排序，保留原始順序

###########################{Volcano Plot}###########################

# Allgene 欄位含 logFC, P.Value, adj.P.Val
# 先把原本的 par 設定存起來，最後再還原
op <- par(no.readonly = TRUE) # no.readonly = TRUE:full list of parameters 

# 調整邊界，並允許在外部繪圖
par(mar = c(5, 4, 4, 6) , xpd = TRUE) #c(bottom, left, top, right)

# 基本火山圖：logFC vs. -log10(P.Value)
with(Allgene, plot(
  logFC, -log10(P.Value),
  pch   = 20,
  main  = "NMOSD vs Healthy DEGs (External data)",
  xlim=c(-15,15),
  col   = "grey",
  xlab  = "log2 Fold Change",
  ylab  = "-log10(P Value)"
))


#----【Significants genes criteria: adj.P.Val < 0.05 & |logFC|>1】-----
# 下調基因 (藍色)
with(subset(Allgene, adj.P.Val < 0.05 & logFC < -1),
     points(logFC, -log10(P.Value), pch=20, col='#8FA4FF')
)

# 上調基因 (紅色)
with(subset(Allgene, adj.P.Val < 0.05 & logFC >1),
     points(logFC, -log10(P.Value), pch=20, col='#FF8080')
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
Up_DEGs<-row.names(subset(Allgene,adj.P.Val < 0.05 & logFC > 1)) 
Down_DEGs<-row.names(subset(Allgene,adj.P.Val < 0.05 & logFC < -1)) 
DEGS<-row.names(subset(Allgene,adj.P.Val<0.05 & abs(logFC)>1)) 
length(Up_DEGs);length(Down_DEGs);length(DEGS)
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS(external data)(limma).csv',row.names = F)
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP(external data)(limma).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN(external data)(limma).csv',row.names = F)



