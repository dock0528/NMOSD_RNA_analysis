setwd("C:/Users/JANE/Desktop/Wang實驗室/NMOSD研究計畫/RNA/NMOSD_RNA_analysis")

###############################【DEGs:Limma+Vomm】###############################
#載入套件
BiocManager::install(c("edgeR","limma"))
library(limma)
library(edgeR)

########################{原始gene matrix}########################
Raw_count_merged<-read.csv('./RNA_DATA/Raw_count_merged_matrix.csv',row.names = 1,header=T) #dim(Raw_count_merged):78724 x 21

#Raw count 
count_df<-as.matrix(Raw_count_merged)
rownames(count_df)<-gsub("\\.\\d+$", "",rownames(count_df)) #去除小數點以後的值

########################{Sample metadata}########################
#Sample & Group
sample_df=data.frame(Sample=colnames(count_df))
sample_df$Group <- ifelse(
  grepl("^SRR", sample_df$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
# 把 sample_df 裡的 Group 欄轉成 factor，Control 當作 baseline
sample_df$Group <- factor(sample_df$Group,
                          levels = c("Control","NMOSD"))

#----【Create DGEList】----
dge <- DGEList(counts=count_df, group=sample_df$Group) #counts:row->genes,col->samples

#----【Filter low expression genes]----
#filterByExpr -> Determine which genes have sufficiently large counts to be retained 
keep <- filterByExpr(dge)      # 自動依 group 大小和 library size 計算過濾閾值
dge  <- dge[keep, , keep.lib.sizes=FALSE] #[保留TRUE的列，所有欄,參數]
#keep.lib.sizes=FALSE不沿用原本的gene sizes，改用filter過後的

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

#【Significants genes criteria: adj.P.Val < 0.01 & |logFC|>2】
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
Up_DEGs<-row.names(subset(Allgene,adj.P.Val < 0.01 & logFC > 2)) #length(Up_DEGs):732
Down_DEGs<-row.names(subset(Allgene,adj.P.Val < 0.01 & logFC < -2)) #length(Down_DEGs):735
DEGS<-row.names(subset(Allgene,adj.P.Val<0.01 & abs(logFC)>2)) #DEGs numbers:1467
write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS(limma).csv',row.names = F)

#-----------------------------------------------------------------------------------------------------------------------------------------------
###############【Housekeeping Genes(11) Map to Original Gene Matrix】###############
#----取RefSeq
housekeeping_genes_df<-read.csv('./RNA_DATA/Housekeeping_genes_for_calibration(11).csv',header=T)
housekeeping_genes<-housekeeping_genes_df$RefSeq.accession.number 

#----RefSeq(NM)轉ENSG
library(biomaRt)
mart <- useEnsembl(
  biomart  = "genes",
  dataset  = "hsapiens_gene_ensembl"
)

mapping <- getBM(
  attributes = c("refseq_mrna",      # NM_ accession
                 "ensembl_gene_id"), # ENSG
  filters    = "refseq_mrna",
  values     = housekeeping_genes,
  mart       = mart
)
head(mapping) #轉換後的ENSG:4267 (Refseq:11 -> ENSG:12)


housekeeping_genes_ENSG_df <- merge(
  housekeeping_genes_df,
  mapping,
  by.x   = "RefSeq.accession.number", #原始欄位
  by.y   = "refseq_mrna", #mapping欄位
  all.x  = TRUE #保留所有原始列
)
#length(housekeeping_genes_ENSG_df$ensembl_gene_id):12(含duplicate ENSG、RefSeq.accession.number)
#write.csv(housekeeping_genes_ENSG_df,file='RNA_DATA/Housekeeping_genes_for_claibration_ENSG.csv',row.names = F)

housekeeping_genes_ENSG_df<-read.csv('RNA_DATA/Housekeeping_genes_for_claibration_ENSG.csv')

#----重複RefSeq的row
#找出重複出現的 RefSeq（只會回傳每個 ID 一次）
dup_refseq <- unique(housekeeping_genes_ENSG_df$RefSeq.accession.number[ duplicated(housekeeping_genes_ENSG_df$RefSeq.accession.number) ])

# 取出重複RefSeq的列
dup_refseq_rows<- housekeeping_genes_ENSG_df[housekeeping_genes_ENSG_df$RefSeq.accession.number %in% dup_refseq, ]



#####################{Control genes(11)}#####################
#----housekeeping ENSG(無NA,無duplicate ENSG)
housekeeping_genes_ENSG_df<-read.csv('RNA_DATA/Housekeeping_genes_for_claibration_ENSG.csv')
housekeeping_ENSG <- housekeeping_genes_ENSG_df$ensembl_gene_id #12個ENSG

#----原始 Gene Matrix
Raw_count_merged<-read.csv('./RNA_DATA/Raw_count_merged_matrix.csv',row.names = 1,header=T) #dim(Raw_count_merged):78724 x 21
#Raw count 
count_df<-as.matrix(Raw_count_merged)
rownames(count_df)<-gsub("\\.\\d+$", "",rownames(count_df)) #去除小數點以後的值 #gene:78724

#----control genes
control_genes<-intersect(rownames(count_df),housekeeping_ENSG) #length(control_genes):11
control_genes_df <- housekeeping_genes_ENSG_df[housekeeping_genes_ENSG_df$ensembl_gene_id %in%control_genes, ] #nrow(control_genes_df):11
#write.csv(control_genes_df,file = 'RNA_DATA/control_genes_df_for_calibration.csv',row.names =F)


###############【DEGs:Limma+Vomm (Batch Effect Correction)】###############
BiocManager::install("sva")
library(sva)
library(limma)
library(edgeR)
#Sample & Group
sample_df=data.frame(Sample=colnames(count_df))
sample_df$Group <- ifelse(
  grepl("^SRR", sample_df$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
# 把 sample_df 裡的 Group 欄轉成 factor，Control 當作 baseline
sample_df$Group <- factor(sample_df$Group,
                          levels = c("Control","NMOSD"))

#----【Create DGEList】----
dge <- DGEList(counts=count_df, group=sample_df$Group) #counts:row->genes,col->samples

#----【Filter low expression genes]----
#filterByExpr -> Determine which genes have sufficiently large counts to be retained 
keep <- filterByExpr(dge)      # 自動依 group 大小和 library size 計算過濾閾值
dge  <- dge[keep, , keep.lib.sizes=FALSE] #[保留TRUE的列，所有欄,參數]
#keep.lib.sizes=FALSE不沿用原本的gene sizes，改用filter過後的

#----【TMM（Trimmed Mean of M-values）】----
#Step1.消除極端表達基因的影響
#Step2.估算每個樣本的「有效 library size」
#Step3.使樣本間的基因表達量可直接比較
dge <- calcNormFactors(dge,method = 'TMM')    # TMM 標準化


#----【Voom轉換】----
#建立 design matrix（不含截距，以 Control 作 baseline）
design <- model.matrix(~ 0 + sample_df$Group)
colnames(design) <- levels(sample_df$Group)  #Control vs NMOSD


#----【將 DGEList 轉成 voom 物件（log2-CPM + 權重）】----
v <- voom(dge, design, plot = TRUE)
# plot = TRUE 會畫出 mean-variance trend，可用來檢查變異度是否隨平均表達量下降

# ----準備 dat & design
log2cpm_matrix  <- v$E                           # voom 轉出的 log2-CPM matrix  
mod  <- model.matrix(~ 0 + sample_df$Group)   # 你的原 design  
colnames(mod) <- levels(sample_df$Group)      # Control, NMOSD  
mod0 <- model.matrix(~ 1, sample_df)          # null model（只含截距）  



#----載入control genes
control_genes_df<-read.csv('RNA_DATA/control_genes_df_for_calibration.csv')
control_genes_ENSG<-control_genes_df$ensembl_gene_id #11
control_genes <- rownames(log2cpm_matrix ) %in% control_genes_ENSG


###########{畫11個housekeeping genes在NMOSD log2CPM的分布}###########
library(reshape2)
library(ggplot2)

log2cpm_control_genes_matrix<-log2cpm_matrix[rownames(log2cpm_matrix) %in% control_genes_ENSG,]
NMOSD_cols <- grepl("Patient", colnames(log2cpm_control_genes_matrix))
NM0SD_log2cpm_control_genes_matrix<-log2cpm_control_genes_matrix[,NMOSD_cols]

#----寬Matrix >>> 長Matrix
long_NMOSD_log2cpm <- melt(
  log2cpm_control_genes_matrix[, NMOSD_cols, drop = FALSE],
  varnames    = c("Gene", "Sample"),
  value.name  = "Log2CPM"
)

#----boxplot distribution in NMOSD patients
ggplot(long_NMOSD_log2cpm, aes(x = Sample, y = Log2CPM)) +
  geom_boxplot(fill = "#EDABC0", outlier.size = 0.5, width = 0.6) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(hjust = 0.5, size = 8),
    axis.title.x = element_blank(),
    plot.title   = element_text(hjust = 0.5, face = "bold",size=14)  # 置中且加粗
  ) +
  labs(
    y     = "log2(CPM)",
    title = "11 Housekeeping genes log2CPM acorss NMOSD patients"
  )

#----Heatmap of wilcox rank sum test between each other
pw <- pairwise.wilcox.test(long_NMOSD_log2cpm$Log2CPM,
                           long_NMOSD_log2cpm$Sample,
                           paired = TRUE,
                           p.adjust.method = "BH")
####################################################
#取pavlue
p.mat <- pw$p.value

#對角設NA
diag(p.mat) <- 1

#畫Heatmap
library(pheatmap)
pheatmap(p.mat,
         display_numbers = TRUE,      #顯示數字
         number_format   = "%.3f",    #數字format
         cluster_rows    = FALSE,
         cluster_cols    = FALSE,
         na_col          = "white",       # the color of the NA cell
         main            = "Heatmap of 11 Housekeeping genes log2CPM acorss NMOSD patients",
         color           = colorRampPalette(c("#74C7EA","white","#E67884"))(50))


library(pheatmap)

# 原始 p 值矩阵
# diag(p.mat) <- NA  # 如果你还想把对角也设成 NA

# 1) 构造标签矩阵：非 NA 就格式化，NA 就 ""
label_mat <- matrix(
  ifelse(is.na(p.mat), "",
         sprintf("%.3f", p.mat)),
  nrow      = nrow(p.mat),
  ncol      = ncol(p.mat),
  dimnames  = dimnames(p.mat)
)

# 2) 画图，用 display_numbers=label_mat
pheatmap(
  p.mat,
  display_numbers = label_mat,               # 自定义标签
  number_color    = "black",
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  na_col          = "white",                 # NA 格子的底色
  main            = "Heatmap of 11 Housekeeping genes log2CPM across NMOSD patients",
  color           = colorRampPalette(c("#74C7EA","white","#E67884"))(50),
  fontsize_number = 8
)
#####################{Svaseq+SSVA}#####################
#----自動估 surrogate variables 的個數
n.sv <- num.sv(log2cpm_matrix , mod, method = "leek")   #多少個batch effect要校正

#----執行 supervised SVA(ssva)
sv_result <- ssva(  
  log2cpm_matrix,  
  controls = control_genes,  
  n.sv      = n.sv  
)  
colnames(sv_result$sv) <- paste0("SV", seq_len(ncol(sv_result$sv))) #加SV欄位名 'SV1' & 'SV2'
#----把 sv 加回 design 
design_sv <- cbind(mod, sv_result$sv)  

#-----重新跑 voom + lmFit（用含 sv 的 design） 
v2   <- voom(dge, design_sv, plot = TRUE)  
fit2 <- lmFit(v2, design_sv)  

#----再做 contrasts + eBayes
contrast.matrix <- makeContrasts(  
  NMOSD_vs_Control = NMOSD - Control,  
  levels           = design_sv  
)  
fit2 <- contrasts.fit(fit2, contrast.matrix)  
fit2 <- eBayes(fit2)  

#----產生全基因結果（不排序）
Allgene_sva <- topTable(  
  fit2,  
  coef         = "NMOSD_vs_Control",  
  adjust.method= "BH",  
  number       = Inf,  
  sort.by      = "none"  
)

#################【DEGs criteria: adjust.P<0.05,abs(log2FC)>1】#################
# Allgene 欄位含 logFC, P.Value, adj.P.Val
# 先把原本的 par 設定存起來，最後再還原
op <- par(no.readonly = TRUE) # no.readonly = TRUE:full list of parameters 

# 調整邊界，並允許在外部繪圖
par(mar = c(5, 4, 4, 6) , xpd = TRUE) #c(bottom, left, top, right)

# 基本火山圖：logFC vs. -log10(P.Value)
with(Allgene_sva, plot(
  logFC, -log10(P.Value),
  pch   = 20,
  main  = "NMOSD vs Healthy DEGs (Batch Effect Correction)",
  #  xlim=c(-15,15),
  col   = "grey",
  xlab  = "log2 Fold Change",
  ylab  = "-log10(P Value)"
))

#【Significants genes criteria: adj.P.Val < 0.01 & |logFC|>2】
# 下調基因 (藍色)
with(subset(Allgene_sva, P.Value < 0.05 & logFC < -1),
     points(logFC, -log10(P.Value), pch=20, col='#8FA4FF')
)

# 上調基因 (紅色)
with(subset(Allgene_sva, P.Value < 0.05 & logFC >1),
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
Up_DEGs<-row.names(subset(Allgene_sva,P.Value < 0.05 & logFC > 1)) #length(Up_DEGs):197
Down_DEGs<-row.names(subset(Allgene_sva,P.Value < 0.05 & logFC < -1)) #length(Down_DEGs):111
DEGS<-row.names(subset(Allgene_sva,P.Value<0.05 & abs(logFC)>1)) #DEGs numbers:308

#################{DEGs進行'BH'校正}#################
#----原始pvalue進行BH校正確保為DEGs
#用原始pvalue取DEGs
DEG_sva<- subset(
  Allgene_sva,
  P.Value  < 0.05 & abs(logFC) > 1
)
#DEG_sva重新做 BH 校正
DEG_sva$adj.P.Val <- p.adjust(DEG_sva$P.Value, method="BH")

#顯著基因(BH校正後)
Up_DEGs<-row.names(subset(DEG_sva,adj.P.Val< 0.05 & logFC > 1)) #length(Up_DEGs):197
Down_DEGs<-row.names(subset(DEG_sva,P.Value < 0.05 & logFC < -1)) #length(Down_DEGs):111
DEGS<-row.names(subset(Allgene_sva,P.Value<0.05 & abs(logFC)>1)) #DEGs numbers:308
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS(Batch Effect Correction).csv',row.names = F)
#--------------------------------------------------


#################【DEGs criteria: adjust.P<0.01,abs(log2FC)>2】#################

op <- par(no.readonly = TRUE) # no.readonly = TRUE:full list of parameters 

# 調整邊界，並允許在外部繪圖
par(mar = c(5, 4, 4, 6) , xpd = TRUE) #c(bottom, left, top, right)

# 基本火山圖：logFC vs. -log10(P.Value)
with(Allgene_sva, plot(
  logFC, -log10(P.Value),
  pch   = 20,
  main  = "NMOSD vs Healthy DEGs (Batch Effect Correction)",
  #  xlim=c(-15,15),
  col   = "grey",
  xlab  = "log2 Fold Change",
  ylab  = "-log10(P Value)"
))

#【Significants genes criteria: adj.P.Val < 0.01 & |logFC|>2】
# 下調基因 (藍色)
with(subset(Allgene_sva, P.Value < 0.01 & logFC < -2),
     points(logFC, -log10(P.Value), pch=20, col='#8FA4FF')
)

# 上調基因 (紅色)
with(subset(Allgene_sva, P.Value < 0.01 & logFC >2),
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
Up_DEGs<-row.names(subset(Allgene_sva,P.Value < 0.01 & logFC > 2)) #length(Up_DEGs):39
Down_DEGs<-row.names(subset(Allgene_sva,P.Value < 0.01 & logFC < -2)) #length(Down_DEGs):12
DEGS<-row.names(subset(Allgene_sva,P.Value<0.01 & abs(logFC)>2)) #DEGs numbers:51

#################{DEGs進行'BH'校正}#################
#----原始pvalue進行BH校正確保為DEGs
#用原始pvalue取DEGs
DEG_sva<- subset(
  Allgene_sva,
  P.Value  < 0.01 & abs(logFC) > 2
)
#DEG_sva重新做 BH 校正
DEG_sva$adj.P.Val <- p.adjust(DEG_sva$P.Value, method="BH")

#顯著基因(BH校正後)
Up_DEGs<-row.names(subset(DEG_sva,adj.P.Val< 0.01 & logFC > 2)) #length(Up_DEGs):39
Down_DEGs<-row.names(subset(DEG_sva,P.Value < 0.01 & logFC < -2)) #length(Down_DEGs):12
DEGS<-row.names(subset(Allgene_sva,P.Value<0.01 & abs(logFC)>2)) #DEGs numbers:51
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS_strictly(Batch Effect Correction).csv',row.names = F)












#---------------------------------------------------------------------------------------------------------------------------------------------------
###############【Housekeeping Genes(3804) Map to Original Gene Matrix】###############
#----取RefSeq
housekeeping_genes_df<-read.csv('./RNA_DATA/Housekeeping_genes(3804).csv',header=T)
housekeeping_genes<-housekeeping_genes_df$RefSeq.accession.number 

#----RefSeq(NM)轉ENSG
library(biomaRt)
mart <- useEnsembl(
  biomart  = "genes",
  dataset  = "hsapiens_gene_ensembl"
)

mapping <- getBM(
  attributes = c("refseq_mrna",      # NM_ accession
                 "ensembl_gene_id"), # ENSG
  filters    = "refseq_mrna",
  values     = housekeeping_genes,
  mart       = mart
)
head(mapping) #轉換後的ENSG:4267 (Refseq:3804 -> ENSG:4267)


housekeeping_genes_ENSG_df <- merge(
  housekeeping_genes_df,
  mapping,
  by.x   = "RefSeq.accession.number", #原始欄位
  by.y   = "refseq_mrna", #mapping欄位
  all.x  = TRUE #保留所有原始列
)
#length(housekeeping_genes_ENSG_df$ensembl_gene_id):4292 (含NA、duplicate ENSG)
#write.csv(housekeeping_genes_ENSG_df,file='RNA_DATA/Housekeeping_genes_ENSG.csv',row.names = F)

#-----{ENSEMBL:處理NA的row}
housekeeping_genes_ENSG_df<- read.csv('RNA_DATA/Housekeeping_genes_ENSG.csv')
ensembl_NA_row<-read.csv('RNA_DATA/Ensembl_housekeeping_genes_NA.csv')
colnames(ensembl_NA_row)[3]<-"ensembl_gene_id"

# 去除空白
housekeeping_genes_ENSG_df$Gene.Name<- trimws(housekeeping_genes_ENSG_df$Gene.Name)
ensembl_NA_row$HGNC.symbol<- trimws(ensembl_NA_row$HGNC.symbol)

ensembl_sub  <- setNames(
  ensembl_NA_row$ensembl_gene_id,
  ensembl_NA_row$HGNC.symbol
)

#找出原 table 哪些列是 NA
na_idx <- is.na(housekeeping_genes_ENSG_df$ensembl_gene_id)

#只把 NA 那些行，用 Hugo Symbol 去取 ensembl_sub 裡的值
housekeeping_genes_ENSG_df$ensembl_gene_id[na_idx] <-
  ensembl_sub [housekeeping_genes_ENSG_df$Gene.Name[na_idx]]

head(housekeeping_genes_ENSG_df)
housekeeping_genes_ENSG_df[is.na(housekeeping_genes_ENSG_df$ensembl_gene_id),]
#------------------------

#---去除NA的row
housekeeping_genes_ENSG_dropna <- housekeeping_genes_ENSG_df[ !is.na(housekeeping_genes_ENSG_df$ensembl_gene_id), ]
# length(housekeeping_genes_ENSG_dropna$ensembl_gene_id):4267 (含duplicate ENSG)
#write.csv(housekeeping_genes_ENSG_dropna,file='RNA_DATA/Housekeeping_genes_ENSG_dropna.csv',row.names = F)

#----重複ENSG的row
#找出重複出現的 ENSG（只會回傳每個 ID 一次）
dup_ENSG <- unique(housekeeping_genes_ENSG_dropna$ensembl_gene_id[ duplicated(housekeeping_genes_ENSG_dropna$ensembl_gene_id) ])

# 取出重複ENSG 的列
dup_ENSG_rows <- housekeeping_genes_ENSG_dropna[housekeeping_genes_ENSG_dropna$ensembl_gene_id %in% dup_ENSG, ]

#----重複RefSeq的row
#找出重複出現的 RefSeq（只會回傳每個 ID 一次）
dup_refseq <- unique(housekeeping_genes_ENSG_dropna$ensembl_gene_id[ duplicated(housekeeping_genes_ENSG_dropna$RefSeq.accession.number) ])

# 取出重複RefSeq的列
dup_refseq_rows<- housekeeping_genes_ENSG_dropna[housekeeping_genes_ENSG_dropna$ensembl_gene_id %in% dup_refseq, ]



#####################{Control genes(3804)}#####################
#----housekeeping ENSG(無NA,無duplicate ENSG)
housekeeping_genes_ENSG_dropna<- read.csv('RNA_DATA/Housekeeping_genes_ENSG_dropna.csv') #row:4267 (有2個重複的ENSG)
housekeeping_ENSG <- housekeeping_genes_ENSG_dropna$ensembl_gene_id[! duplicated(housekeeping_genes_ENSG_dropna$ensembl_gene_id)] #length(housekeeping_ENSG):4265
housekeeping_ENSG_dropna_dropdup <- housekeeping_genes_ENSG_dropna[! duplicated(housekeeping_genes_ENSG_dropna$ensembl_gene_id),]


#----原始 Gene Matrix
Raw_count_merged<-read.csv('./RNA_DATA/Raw_count_merged_matrix.csv',row.names = 1,header=T) #dim(Raw_count_merged):78724 x 21
#Raw count 
count_df<-as.matrix(Raw_count_merged)
rownames(count_df)<-gsub("\\.\\d+$", "",rownames(count_df)) #去除小數點以後的值 #gene:78724

#----control genes
control_genes<-intersect(rownames(count_df),housekeeping_ENSG) #length(control_genes):3776
control_genes_df <- housekeeping_ENSG_dropna_dropdup[housekeeping_ENSG_dropna_dropdup$ensembl_gene_id %in%control_genes, ] #nrow(control_genes_df):3776
#write.csv(control_genes_df,file = 'RNA_DATA/control_genes_df.csv',row.names =F)

###############【DEGs:Limma+Vomm (Batch Effect Correction)】###############
BiocManager::install("sva")
library(sva)
library(limma)
library(edgeR)
#Sample & Group
sample_df=data.frame(Sample=colnames(count_df))
sample_df$Group <- ifelse(
  grepl("^SRR", sample_df$Sample),  # 如果 Sample 以 "SRR" 開頭
  "Control",                         # 就標成 Control
  "NMOSD"                            # 否則都當 NMOSD
)
# 把 sample_df 裡的 Group 欄轉成 factor，Control 當作 baseline
sample_df$Group <- factor(sample_df$Group,
                          levels = c("Control","NMOSD"))

#----【Create DGEList】----
dge <- DGEList(counts=count_df, group=sample_df$Group) #counts:row->genes,col->samples

#----【Filter low expression genes]----
#filterByExpr -> Determine which genes have sufficiently large counts to be retained 
keep <- filterByExpr(dge)      # 自動依 group 大小和 library size 計算過濾閾值
dge  <- dge[keep, , keep.lib.sizes=FALSE] #[保留TRUE的列，所有欄,參數]
#keep.lib.sizes=FALSE不沿用原本的gene sizes，改用filter過後的

#----【TMM（Trimmed Mean of M-values）】----
#Step1.消除極端表達基因的影響
#Step2.估算每個樣本的「有效 library size」
#Step3.使樣本間的基因表達量可直接比較
dge <- calcNormFactors(dge,method = 'TMM')    # TMM 標準化


#----【Voom轉換】----
#建立 design matrix（不含截距，以 Control 作 baseline）
design <- model.matrix(~ 0 + sample_df$Group)
colnames(design) <- levels(sample_df$Group)  #Control vs NMOSD


#----【將 DGEList 轉成 voom 物件（log2-CPM + 權重）】----
v <- voom(dge, design, plot = TRUE)
# plot = TRUE 會畫出 mean-variance trend，可用來檢查變異度是否隨平均表達量下降

# ----準備 dat & design
log2cpm_matrix  <- v$E                           # voom 轉出的 log2-CPM matrix  
mod  <- model.matrix(~ 0 + sample_df$Group)   # 你的原 design  
colnames(mod) <- levels(sample_df$Group)      # Control, NMOSD  
mod0 <- model.matrix(~ 1, sample_df)          # null model（只含截距）  

#----載入control genes
control_genes_df<-read.csv('RNA_DATA/control_genes_df.csv')
control_genes_ENSG<-control_genes_df$ensembl_gene_id #3776
control_genes <- rownames(log2cpm_matrix ) %in% control_genes_ENSG

#####################{Svaseq+SSVA}#####################
#----自動估 surrogate variables 的個數
n.sv <- num.sv(log2cpm_matrix , mod, method = "leek")   #多少個batch effect要校正

#----執行 supervised SVA(ssva)
sv_result <- ssva(  
  log2cpm_matrix,  
  controls = control_genes,  
  n.sv      = n.sv  
)  
colnames(sv_result$sv) <- paste0("SV", seq_len(ncol(sv_result$sv))) #加SV欄位名 'SV1' & 'SV2'
#----把 sv 加回 design 
design_sv <- cbind(mod, sv_result$sv)  

#-----重新跑 voom + lmFit（用含 sv 的 design） 
v2   <- voom(dge, design_sv, plot = TRUE)  
fit2 <- lmFit(v2, design_sv)  

#----再做 contrasts + eBayes
contrast.matrix <- makeContrasts(  
  NMOSD_vs_Control = NMOSD - Control,  
  levels           = design_sv  
)  
fit2 <- contrasts.fit(fit2, contrast.matrix)  
fit2 <- eBayes(fit2)  

#----產生全基因結果（不排序）
Allgene_sva <- topTable(  
  fit2,  
  coef         = "NMOSD_vs_Control",  
  adjust.method= "BH",  
  number       = Inf,  
  sort.by      = "none"  
)

#################【DEGs criteria: adjust.P<0.05,abs(log2FC)>1】#################
# Allgene 欄位含 logFC, P.Value, adj.P.Val
# 先把原本的 par 設定存起來，最後再還原
op <- par(no.readonly = TRUE) # no.readonly = TRUE:full list of parameters 

# 調整邊界，並允許在外部繪圖
par(mar = c(5, 4, 4, 6) , xpd = TRUE) #c(bottom, left, top, right)

# 基本火山圖：logFC vs. -log10(P.Value)
with(Allgene_sva, plot(
  logFC, -log10(P.Value),
  pch   = 20,
  main  = "NMOSD vs Healthy DEGs (Batch Effect Correction)",
#  xlim=c(-15,15),
  col   = "grey",
  xlab  = "log2 Fold Change",
  ylab  = "-log10(P Value)"
))

#【Significants genes criteria: adj.P.Val < 0.01 & |logFC|>2】
# 下調基因 (藍色)
with(subset(Allgene_sva, P.Value < 0.05 & logFC < -1),
     points(logFC, -log10(P.Value), pch=20, col='#8FA4FF')
)

# 上調基因 (紅色)
with(subset(Allgene_sva, P.Value < 0.05 & logFC >1),
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
Up_DEGs<-row.names(subset(Allgene_sva,P.Value < 0.05 & logFC > 1)) #length(Up_DEGs):197
Down_DEGs<-row.names(subset(Allgene_sva,P.Value < 0.05 & logFC < -1)) #length(Down_DEGs):111
DEGS<-row.names(subset(Allgene_sva,P.Value<0.05 & abs(logFC)>1)) #DEGs numbers:308

#################{DEGs進行'BH'校正}#################
#----原始pvalue進行BH校正確保為DEGs
#用原始pvalue取DEGs
DEG_sva<- subset(
  Allgene_sva,
  P.Value  < 0.05 & abs(logFC) > 1
)
#DEG_sva重新做 BH 校正
DEG_sva$adj.P.Val <- p.adjust(DEG_sva$P.Value, method="BH")

#顯著基因(BH校正後)
Up_DEGs<-row.names(subset(DEG_sva,adj.P.Val< 0.05 & logFC > 1)) #length(Up_DEGs):197
Down_DEGs<-row.names(subset(DEG_sva,P.Value < 0.05 & logFC < -1)) #length(Down_DEGs):111
DEGS<-row.names(subset(Allgene_sva,P.Value<0.05 & abs(logFC)>1)) #DEGs numbers:308
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS(Batch Effect Correction).csv',row.names = F)
#--------------------------------------------------


#################【DEGs criteria: adjust.P<0.01,abs(log2FC)>2】#################

op <- par(no.readonly = TRUE) # no.readonly = TRUE:full list of parameters 

# 調整邊界，並允許在外部繪圖
par(mar = c(5, 4, 4, 6) , xpd = TRUE) #c(bottom, left, top, right)

# 基本火山圖：logFC vs. -log10(P.Value)
with(Allgene_sva, plot(
  logFC, -log10(P.Value),
  pch   = 20,
  main  = "NMOSD vs Healthy DEGs (Batch Effect Correction)",
  #  xlim=c(-15,15),
  col   = "grey",
  xlab  = "log2 Fold Change",
  ylab  = "-log10(P Value)"
))

#【Significants genes criteria: adj.P.Val < 0.01 & |logFC|>2】
# 下調基因 (藍色)
with(subset(Allgene_sva, P.Value < 0.01 & logFC < -2),
     points(logFC, -log10(P.Value), pch=20, col='#8FA4FF')
)

# 上調基因 (紅色)
with(subset(Allgene_sva, P.Value < 0.01 & logFC >2),
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
Up_DEGs<-row.names(subset(Allgene_sva,P.Value < 0.01 & logFC > 2)) #length(Up_DEGs):39
Down_DEGs<-row.names(subset(Allgene_sva,P.Value < 0.01 & logFC < -2)) #length(Down_DEGs):12
DEGS<-row.names(subset(Allgene_sva,P.Value<0.01 & abs(logFC)>2)) #DEGs numbers:51

#################{DEGs進行'BH'校正}#################
#----原始pvalue進行BH校正確保為DEGs
#用原始pvalue取DEGs
DEG_sva<- subset(
  Allgene_sva,
  P.Value  < 0.01 & abs(logFC) > 2
)
#DEG_sva重新做 BH 校正
DEG_sva$adj.P.Val <- p.adjust(DEG_sva$P.Value, method="BH")

#顯著基因(BH校正後)
Up_DEGs<-row.names(subset(DEG_sva,adj.P.Val< 0.01 & logFC > 2)) #length(Up_DEGs):39
Down_DEGs<-row.names(subset(DEG_sva,P.Value < 0.01 & logFC < -2)) #length(Down_DEGs):12
DEGS<-row.names(subset(Allgene_sva,P.Value<0.01 & abs(logFC)>2)) #DEGs numbers:51
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS_strictly(Batch Effect Correction).csv',row.names = F)
