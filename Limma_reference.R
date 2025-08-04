setwd("C:/Users/Jane/Desktop/Wang實驗室/NMOSD研究計畫/RNA/NMOSD_RNA_analysis")
library(limma)
library(edgeR)

########################{原始gene matrix}########################
Raw_count_merged<-read.csv('./RNA_DATA/External_RNA_datasets_DEG_RESULT/RNA_gene_matrix_v2.csv',row.names = 1,header=T) #dim(Raw_count_merged):58304 x 23
Raw_count_merged <- Raw_count_merged[!rownames(Raw_count_merged) %in% c("Antibody status","Ensembl_Gene_ID"), ]

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
#################【DEGs criteria: adjust.P<0.05,abs(log2FC)>1】#################
# Allgene 欄位含 logFC, P.Value, adj.P.Val
# 先把原本的 par 設定存起來，最後再還原
op <- par(no.readonly = TRUE) # no.readonly = TRUE:full list of parameters 

# 調整邊界，並允許在外部繪圖
par(mar = c(5, 4, 4, 6) , xpd = TRUE) #c(bottom, left, top, right)

# 基本火山圖：logFC vs. -log10(P.Value)
with(Allgene, plot(
  logFC, -log10(P.Value),
  pch   = 20,
  main  = "NMOSD vs Healthy DEGs (External Dataset)",
  #  xlim=c(-15,15),
  col   = "grey",
  xlab  = "log2 Fold Change",
  ylab  = "-log10(P Value)"
))

#【Significants genes criteria: adj.P.Val < 0.05 & |logFC|>1】
# 下調基因 (藍色)
with(subset(Allgene, P.Value < 0.05 & logFC < -1),
     points(logFC, -log10(P.Value), pch=20, col='#8FA4FF')
)

# 上調基因 (紅色)
with(subset(Allgene, P.Value < 0.05 & logFC >1),
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
Up_DEGs<-row.names(subset(Allgene,P.Value < 0.05 & logFC > 1)) #length(Up_DEGs):90
Down_DEGs<-row.names(subset(Allgene,P.Value < 0.05 & logFC < -1)) #length(Down_DEGs):277
DEGS<-row.names(subset(Allgene,P.Value<0.05 & abs(logFC)>1)) #length(DEGS):367

#################{DEGs進行'BH'校正}#################
#----原始pvalue進行BH校正確保為DEGs
#用原始pvalue取DEGs
DEG_BH<- subset(
  Allgene,
  P.Value  < 0.05 & abs(logFC) > 1
)
#DEG_sva重新做 BH 校正
DEG_BH$adj.P.Val <- p.adjust(DEG_BH$P.Value, method="BH")

#顯著基因(BH校正後)
Up_DEGs<-row.names(subset(DEG_BH,adj.P.Val< 0.05 & logFC > 1)) #length(Up_DEGs):90
Down_DEGs<-row.names(subset(DEG_BH,adj.P.Val < 0.05 & logFC < -1)) #length(Down_DEGs):277
DEGS<-row.names(subset(DEG_BH,adj.P.Val<0.05 & abs(logFC)>1)) #length(DEGS):367
#write.csv(Up_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_UP(External dataset).csv',row.names = F)
#write.csv(Down_DEGs,file='RNA_DATA/NMOSD_RNA_DEGS_DOWN(External dataset).csv',row.names = F)

#################【DEGs criteria: adjust.P<0.01,abs(log2FC)>2】#################

op <- par(no.readonly = TRUE) # no.readonly = TRUE:full list of parameters 

# 調整邊界，並允許在外部繪圖
par(mar = c(5, 4, 4, 6) , xpd = TRUE) #c(bottom, left, top, right)

# 基本火山圖：logFC vs. -log10(P.Value)
with(Allgene, plot(
  logFC, -log10(P.Value),
  pch   = 20,
  main  = "NMOSD vs Healthy DEGs (External Dataset)",
  #  xlim=c(-15,15),
  col   = "grey",
  xlab  = "log2 Fold Change",
  ylab  = "-log10(P Value)"
))

#【Significants genes criteria: adj.P.Val < 0.01 & |logFC|>2】
# 下調基因 (藍色)
with(subset(Allgene, P.Value < 0.01 & logFC < -2),
     points(logFC, -log10(P.Value), pch=20, col='#8FA4FF')
)

# 上調基因 (紅色)
with(subset(Allgene, P.Value < 0.01 & logFC >2),
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
Up_DEGs<-row.names(subset(Allgene,P.Value < 0.01 & logFC > 2)) #length(Up_DEGs):7
Down_DEGs<-row.names(subset(Allgene,P.Value < 0.01 & logFC < -2)) #length(Down_DEGs):13
DEGS<-row.names(subset(Allgene,P.Value<0.01 & abs(logFC)>2)) #length(DEGS):20

#################{DEGs進行'BH'校正}#################
#----原始pvalue進行BH校正確保為DEGs
#用原始pvalue取DEGs
DEG_BH<- subset(
  Allgene,
  P.Value  < 0.01 & abs(logFC) > 2
)
#DEG_sva重新做 BH 校正
DEG_BH$adj.P.Val <- p.adjust(DEG_BH$P.Value, method="BH")

#顯著基因(BH校正後)
Up_DEGs<-row.names(subset(DEG_BH,adj.P.Val< 0.01 & logFC > 2)) #length(Up_DEGs):7
Down_DEGs<-row.names(subset(DEG_BH,adj.P.Val < 0.01 & logFC < -2)) #length(Down_DEGs):13
DEGS<-row.names(subset(DEG_BH,adj.P.Val<0.01 & abs(logFC)>2)) #length(DEGS):20
#write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS_strictly(External dataset).csv',row.names = F)


######################【External DEGs in Type I IFN】#######################
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
###【GO】
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(External dataset).csv',header=T)

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
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes) #和DEGs交集有99
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有20


external_hallmark_interferon_in_DEGs<-c("ISG15", "IFI44L", "IFI44", "CMPK2", "RSAD2", "LAMP3", "RTP4", "CXCL10",
"DDX60", "LY6E", "IFITM3", "BATF2", "IFIT2", "IFIT3", "OAS1", "OASL", "EPSTI1", "IFI27", "USP18", "MX1") #20

#-----------------【External DEGs 交集 Control genes DEGs】---------------------------
#################【DEGs criteria: adjust.P<0.05,abs(log2FC)>1】#################
#----{10 & 11 control genes & External dataset}
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
###【GO】
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(External dataset).csv',header=T)

DEG_genes<-DEGs_df$x #367個
DEGs_10<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(10Control Genes Batch Effect Correction).csv',header=T)$x #DEGs count3149
DEGs_11<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(11Control Genes Batch Effect Correction).csv',header=T)$x #DEGs count3247
external_with10_DEGs<-intersect(DEG_genes,DEGs_10)
external_with11_DEGs<-intersect(DEG_genes,DEGs_11)
library(VennDiagram)
venn.list3<-list(`10 Control genes`=DEGs_10,`11 Control genes`=DEGs_11
                 ,`External dataset`=DEG_genes)

venn.plot3 <- venn.diagram(
  x          = venn.list3,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Different housekeeping DEGs with External Dataset DEGs',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.52), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFB3C0","#92D5FF","#C0CBDF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFB3C0","#92D5FF","#C0CBDF"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 1.4,#每個圓圈裡數字大小
  cat.dist   = c(0.05,0.11,0.08), #圓圈標題距離圓圈
  cat.default.pos='outer',
  margin     = 2, #空白區域(越大越多)
  cat.pos=c(-45,45,0)# 顯示
)
grid.newpage()
grid.draw(venn.plot3)

#----{Origin & 10 control genes & External dataset}
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(External dataset).csv',header=T)
DEG_genes<-DEGs_df$x #367個
DEGs_10<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(10Control Genes Batch Effect Correction).csv',header=T)$x #DEGs count3149
DEGs_ori<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_origin(limma).csv',header=T)$x #DEGs count3149
external_with10_DEGs<-intersect(DEG_genes,DEGs_10)
external_with_ori_DEGs<-intersect(DEG_genes,DEGs_ori)
library(VennDiagram)
venn.list3<-list(`10 Control genes`=DEGs_10,`Origin dataset`=DEGs_ori
                 ,`External dataset`=DEG_genes)

venn.plot3 <- venn.diagram(
  x          = venn.list3,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Different housekeeping DEGs with External Dataset DEGs',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.58), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFB3C0","#FFDF80","#C0CBDF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFB3C0","#FFDF80","#C0CBDF"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 2,#每個圓圈裡數字大小
  cat.dist   = c(0.07,0.07,0.03), #圓圈標題距離圓圈
  cat.default.pos='outer',
  margin     = 2, #空白區域(越大越多)
  cat.pos=c(-45,45,-180)# 顯示
)
grid.newpage()
grid.draw(venn.plot3)

#----{1285 & 2268 & 3371 control genes & External dataset}
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(External dataset).csv',header=T)

DEG_genes<-DEGs_df$x #367個
DEGs_1285<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(1285Control Genes Batch Effect Correction).csv',header=T)$x #DEGs count355
DEGs_2268<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(2268Control Genes Batch Effect Correction).csv',header=T)$x #DEGs count322
DEGs_3371<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(3371Control Genes Batch Effect Correction).csv',header=T)$x #DEGs count314
external_with1285_DEGs<-intersect(DEG_genes,DEGs_1285)
external_with2268_DEGs<-intersect(DEG_genes,DEGs_2268)
external_with3371_DEGs<-intersect(DEG_genes,DEGs_3371)
library(VennDiagram)
venn.list4<-list(`1285 Control genes`=DEGs_1285,`2268 Control genes`=DEGs_2268
                 ,`3371 Control genes`=DEGs_3371,`External dataset`=DEG_genes)

venn.plot4 <- venn.diagram(
  x          = venn.list4,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Different housekeeping DEGs with External Dataset DEGs',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.6), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#ffc18c","#D3C3FF","#8EFFE4","#C0CBDF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#ffc18c","#D3C3FF","#8EFFE4","#C0CBDF"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 1.4,#每個圓圈裡數字大小
  #  cat.dist   = c(0.075,0.075,0.03), #圓圈標題距離圓圈
  cat.default.pos='outer',
  margin     = 2, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot4)

####################【畫2圓交集圖】#################################
#---{10 & EXTERNAL}
library(VennDiagram)
###【GO】
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(External dataset).csv',header=T)

DEG_genes<-DEGs_df$x #367個
DEGs_10<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(10Control Genes Batch Effect Correction).csv',header=T)$x #DEGs count3149
external_with10_DEGs<-intersect(DEG_genes,DEGs_10)
venn.list2<-list(`10 Control genes`=DEGs_10,`External dataset`=DEG_genes)
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='10 housekeeping DEGs with External Dataset DEGs',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.61), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFB3C0","#C0CBDF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFB3C0","#C0CBDF"),
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



#---{ORIGIN & EXTERNAL}
library(VennDiagram)
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(External dataset).csv',header=T)
DEG_genes<-DEGs_df$x #367個
DEGs_ori<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_origin(limma).csv',header=T)$x #DEGs count3149
external_with_ori_DEGs<-intersect(DEG_genes,DEGs_ori)
venn.list2<-list(`Origin dataset`=DEGs_ori,`External dataset`=DEG_genes)
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Origin DEGs with External Dataset DEGs',
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

#################【DEGs criteria: adjust.P<0.01,abs(log2FC)>2】#################
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
###【GO】
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(External dataset).csv',header=T)

DEG_genes<-DEGs_df$x #367個
DEGs_1285<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(1285Control Genes Batch Effect Correction).csv',header=T)$x #DEGs count355
DEGs_2268<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(2268Control Genes Batch Effect Correction).csv',header=T)$x #DEGs count322
DEGs_3371<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(3371Control Genes Batch Effect Correction).csv',header=T)$x #DEGs count314
library(VennDiagram)
venn.list4<-list(`1285 Control genes`=DEGs_1285,`2268 Control genes`=DEGs_2268
                 ,`3371 Control genes`=DEGs_3371,`External dataset`=DEG_genes)

venn.plot4 <- venn.diagram(
  x          = venn.list4,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Different housekeeping DEGs with External Dataset DEGs',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.6), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#ffc18c","#D3C3FF","#8EFFE4","#C0CBDF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#ffc18c","#D3C3FF","#8EFFE4","#C0CBDF"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 1.4,#每個圓圈裡數字大小
  #  cat.dist   = c(0.075,0.075,0.03), #圓圈標題距離圓圈
  cat.default.pos='outer',
  margin     = 2, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot4)

#---------------【Type I IFN related-DEGs (criteria:padj<0.05 & |log2FC|>1)】------------------
####################【畫2圓交集圖】#################################
#----{Origin data & External data}
library(VennDiagram)
type.I.IFN.DEGs<-read.csv("./RNA_DATA/Type_I_IFN_DEG_genes(origin data).csv",check.names = FALSE)
IFN_ori<-type.I.IFN.DEGs$x
external_interferon_in_DEGs<-c("ISG15",   "IFI6",    "IFI44L",  "IFI44",   "MCOLN3",  "GBP1",    "GBP5",    "FCRL2",   "FCRL1",   "FCER1A",
                               "FCRLA",   "LAMC1",   "NMNAT2",  "CMPK2",   "RSAD2",   "EFR3B",   "TNFAIP6", "COBLL1",  "SLC38A11", "PLCD4",
                               "KLHL30",  "IL17RE",  "IL17RB",  "CD200R1L","CPA3",    "RTP4",    "UGT2B17", "CXCL10",  "HERC5",    "DDX60",
                               "GPM6A",   "IL4",     "SYNPO",   "EBF1",    "RCAN2",   "LGSN",    "ELOVL4",  "AKAP12",  "GPR146",   "THSD7A",
                               "ITGB8",   "NRCAM",   "WNT16",   "ALAS2",   "XKRX",    "CA8",     "LY6E",    "CD274",   "AQP7",     "CD72",
                               "IFITM3",  "SLC22A8", "BATF2",   "MMP8",    "POU2AF1", "KIAA1217","IFIT2",   "IFIT3",   "IFIT1",    "IFIT5",
                               "FGFR2",   "ADAM12",  "BHLHE41", "RAPGEF3", "IGFBP6",  "SYT1",    "NTN4",    "OAS1",    "OAS3",     "OAS2",
                               "OASL",    "EPSTI1",  "SLC7A8",  "IFI27",   "IGHM",    "GATM",    "FAM81A",  "LIAT1",   "XAF1",     "TNFRSF13B",
                               "EPN3",    "HRH4",    "KLHL14",  "NETO1",   "SIGLEC1", "PKIG",    "SLC12A5", "PTPRS",   "FCER2",    "CD22",
                               "CD79A",   "CRX",     "LILRA4",  "KIR3DL1", "USP18",   "LARGE1",  "ELFN2",   "TNFRSF13C","MX1")
venn.list2<-list(`Origin dataset`=IFN_ori,`External dataset`=external_interferon_in_DEGs)
ori_external_IFN_DEGs<-intersect(IFN_ori,external_interferon_in_DEGs)
IFN_DEG_genes_df<-bitr(ori_external_IFN_DEGs,fromType = 'SYMBOL',toType=c('ENTREZID','ENSEMBL'),OrgDb='org.Hs.eg.db')
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Origin dataset type I IFN related-DEGs with External Dataset',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.61), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFDF80","#C0CBDF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFDF80","#C0CBDF"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 1.6,#每個圓圈裡數字大小
  cat.dist   = 0.04, #圓圈標題距離圓圈
  cat.pos    = c(-30,30), #圓圈標題角度 #往右正
  margin     = 1, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot2)
#----{10 & External data}
library(VennDiagram)
type.I.IFN.DEGs<-read.csv("./RNA_DATA/Type_I_IFN_DEG_genes_all.csv",check.names = FALSE)
control_gene10<-type.I.IFN.DEGs$`10`
external_interferon_in_DEGs<-c("ISG15",   "IFI6",    "IFI44L",  "IFI44",   "MCOLN3",  "GBP1",    "GBP5",    "FCRL2",   "FCRL1",   "FCER1A",
                               "FCRLA",   "LAMC1",   "NMNAT2",  "CMPK2",   "RSAD2",   "EFR3B",   "TNFAIP6", "COBLL1",  "SLC38A11", "PLCD4",
                               "KLHL30",  "IL17RE",  "IL17RB",  "CD200R1L","CPA3",    "RTP4",    "UGT2B17", "CXCL10",  "HERC5",    "DDX60",
                               "GPM6A",   "IL4",     "SYNPO",   "EBF1",    "RCAN2",   "LGSN",    "ELOVL4",  "AKAP12",  "GPR146",   "THSD7A",
                               "ITGB8",   "NRCAM",   "WNT16",   "ALAS2",   "XKRX",    "CA8",     "LY6E",    "CD274",   "AQP7",     "CD72",
                               "IFITM3",  "SLC22A8", "BATF2",   "MMP8",    "POU2AF1", "KIAA1217","IFIT2",   "IFIT3",   "IFIT1",    "IFIT5",
                               "FGFR2",   "ADAM12",  "BHLHE41", "RAPGEF3", "IGFBP6",  "SYT1",    "NTN4",    "OAS1",    "OAS3",     "OAS2",
                               "OASL",    "EPSTI1",  "SLC7A8",  "IFI27",   "IGHM",    "GATM",    "FAM81A",  "LIAT1",   "XAF1",     "TNFRSF13B",
                               "EPN3",    "HRH4",    "KLHL14",  "NETO1",   "SIGLEC1", "PKIG",    "SLC12A5", "PTPRS",   "FCER2",    "CD22",
                               "CD79A",   "CRX",     "LILRA4",  "KIR3DL1", "USP18",   "LARGE1",  "ELFN2",   "TNFRSF13C","MX1")

venn.list2<-list(`10 Control genes`=control_gene10,`External dataset`=external_interferon_in_DEGs)
intersect(control_gene10,external_interferon_in_DEGs)
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='10 housekeeping type I IFN related-DEGs with External Dataset',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.61), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFB3C0","#C0CBDF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFB3C0","#C0CBDF"),
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

####################【畫3圓交集圖】#################################
#install.packages("VennDiagram")
library(VennDiagram)
type.I.IFN.DEGs<-read.csv("./RNA_DATA/Type_I_IFN_DEG_genes_all.csv",check.names = FALSE)
control_gene10<-type.I.IFN.DEGs$`10`
control_gene11<-type.I.IFN.DEGs$`11`
external_interferon_in_DEGs<-c("ISG15",   "IFI6",    "IFI44L",  "IFI44",   "MCOLN3",  "GBP1",    "GBP5",    "FCRL2",   "FCRL1",   "FCER1A",
                      "FCRLA",   "LAMC1",   "NMNAT2",  "CMPK2",   "RSAD2",   "EFR3B",   "TNFAIP6", "COBLL1",  "SLC38A11", "PLCD4",
                      "KLHL30",  "IL17RE",  "IL17RB",  "CD200R1L","CPA3",    "RTP4",    "UGT2B17", "CXCL10",  "HERC5",    "DDX60",
                      "GPM6A",   "IL4",     "SYNPO",   "EBF1",    "RCAN2",   "LGSN",    "ELOVL4",  "AKAP12",  "GPR146",   "THSD7A",
                      "ITGB8",   "NRCAM",   "WNT16",   "ALAS2",   "XKRX",    "CA8",     "LY6E",    "CD274",   "AQP7",     "CD72",
                      "IFITM3",  "SLC22A8", "BATF2",   "MMP8",    "POU2AF1", "KIAA1217","IFIT2",   "IFIT3",   "IFIT1",    "IFIT5",
                      "FGFR2",   "ADAM12",  "BHLHE41", "RAPGEF3", "IGFBP6",  "SYT1",    "NTN4",    "OAS1",    "OAS3",     "OAS2",
                      "OASL",    "EPSTI1",  "SLC7A8",  "IFI27",   "IGHM",    "GATM",    "FAM81A",  "LIAT1",   "XAF1",     "TNFRSF13B",
                      "EPN3",    "HRH4",    "KLHL14",  "NETO1",   "SIGLEC1", "PKIG",    "SLC12A5", "PTPRS",   "FCER2",    "CD22",
                      "CD79A",   "CRX",     "LILRA4",  "KIR3DL1", "USP18",   "LARGE1",  "ELFN2",   "TNFRSF13C","MX1")

venn.list3<-list(`10 Control genes`=control_gene10,`11 Control genes`=control_gene11,`External dataset`=external_interferon_in_DEGs)

venn.plot3 <- venn.diagram(
  x          = venn.list3,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Type I IFN related-DEGs with External Dataset',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.52), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFB3C0","#92D5FF","#BFC0E3"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFB3C0","#92D5FF","#BFC0E3"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 1.4,#每個圓圈裡數字大小
  cat.dist   = c(0.075,0.075,0.03), #圓圈標題距離圓圈
  cat.default.pos='outer',
  margin     = 2, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot3)

#----{算出3個交集的}
Reduce(intersect,
       list(control_gene10,
            control_gene11,
            external_interferon_in_DEGs)) #9個

####################【畫4圓交集圖】#################################
control_gene1285<-type.I.IFN.DEGs$`1285`
control_gene2268<-type.I.IFN.DEGs$`2268`
control_gene3371<-type.I.IFN.DEGs$`3371`
external_interferon_in_DEGs<-c("ISG15",   "IFI6",    "IFI44L",  "IFI44",   "MCOLN3",  "GBP1",    "GBP5",    "FCRL2",   "FCRL1",   "FCER1A",
                               "FCRLA",   "LAMC1",   "NMNAT2",  "CMPK2",   "RSAD2",   "EFR3B",   "TNFAIP6", "COBLL1",  "SLC38A11", "PLCD4",
                               "KLHL30",  "IL17RE",  "IL17RB",  "CD200R1L","CPA3",    "RTP4",    "UGT2B17", "CXCL10",  "HERC5",    "DDX60",
                               "GPM6A",   "IL4",     "SYNPO",   "EBF1",    "RCAN2",   "LGSN",    "ELOVL4",  "AKAP12",  "GPR146",   "THSD7A",
                               "ITGB8",   "NRCAM",   "WNT16",   "ALAS2",   "XKRX",    "CA8",     "LY6E",    "CD274",   "AQP7",     "CD72",
                               "IFITM3",  "SLC22A8", "BATF2",   "MMP8",    "POU2AF1", "KIAA1217","IFIT2",   "IFIT3",   "IFIT1",    "IFIT5",
                               "FGFR2",   "ADAM12",  "BHLHE41", "RAPGEF3", "IGFBP6",  "SYT1",    "NTN4",    "OAS1",    "OAS3",     "OAS2",
                               "OASL",    "EPSTI1",  "SLC7A8",  "IFI27",   "IGHM",    "GATM",    "FAM81A",  "LIAT1",   "XAF1",     "TNFRSF13B",
                               "EPN3",    "HRH4",    "KLHL14",  "NETO1",   "SIGLEC1", "PKIG",    "SLC12A5", "PTPRS",   "FCER2",    "CD22",
                               "CD79A",   "CRX",     "LILRA4",  "KIR3DL1", "USP18",   "LARGE1",  "ELFN2",   "TNFRSF13C","MX1")

venn.list4<-list(`1285 Control genes`=control_gene1285,`2268 Control genes`=control_gene2268
                 ,`3371 Control genes`=control_gene3371,`External dataset`=external_interferon_in_DEGs)

venn.plot4 <- venn.diagram(
  x          = venn.list4,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Type I IFN related-DEGs with External Dataset',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.6), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#ffc18c","#D3C3FF","#8EFFE4","#BFC0E3"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#ffc18c","#D3C3FF","#8EFFE4","#BFC0E3"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 1.4,#每個圓圈裡數字大小
#  cat.dist   = c(0.075,0.075,0.03), #圓圈標題距離圓圈
  cat.default.pos='outer',
  margin     = 2, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot4)






#---------------【Hallmark Type I IFN related-DEGs (criteria:padj<0.05 & |log2FC|>1)】------------------
####################【畫2圓交集圖】#################################
#----{Origin data & External data}
library(VennDiagram)
hallmark.IFN.DEGs<-read.csv("./RNA_DATA/Hallmark_IFN_DEG_genes(origin data).csv",check.names = FALSE)
hallmark_IFN_ori<-hallmark.IFN.DEGs$x
external_hallmark_interferon_in_DEGs<-c("ISG15", "IFI44L", "IFI44", "CMPK2", "RSAD2", "LAMP3", "RTP4", "CXCL10",
                                        "DDX60", "LY6E", "IFITM3", "BATF2", "IFIT2", "IFIT3", "OAS1", "OASL", "EPSTI1", "IFI27", "USP18", "MX1") #20
venn.list2<-list(`Origin dataset`=hallmark_IFN_ori,`External dataset`=external_hallmark_interferon_in_DEGs)
ori_external_hallmark.IFN_DEGs<-intersect(hallmark_IFN_ori,external_hallmark_interferon_in_DEGs)
hallmark.IFN_DEG_genes_df<-bitr(ori_external_hallmark.IFN_DEGs,fromType = 'SYMBOL',toType=c('ENTREZID','ENSEMBL'),OrgDb='org.Hs.eg.db')
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Origin dataset Hallmark Type I IFN DEGs with External Dataset',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.61), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFDF80","#C0CBDF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFDF80","#C0CBDF"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 1.6,#每個圓圈裡數字大小
  cat.dist   = 0.04, #圓圈標題距離圓圈
  cat.pos    = c(-30,30), #圓圈標題角度 #往右正
  margin     = 1, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot2)
####################【畫3圓交集圖】#################################
library(VennDiagram)
hallmark.IFN.DEGs<-read.csv("./RNA_DATA/Hallmark_IFN_DEG_genes_all.csv",check.names = FALSE)
external_hallmark_interferon_in_DEGs<-c("ISG15", "IFI44L", "IFI44", "CMPK2", "RSAD2", "LAMP3", "RTP4", "CXCL10",
                                        "DDX60", "LY6E", "IFITM3", "BATF2", "IFIT2", "IFIT3", "OAS1", "OASL", "EPSTI1", "IFI27", "USP18", "MX1") #20
control_gene10<-hallmark.IFN.DEGs$`10`
control_gene11<-hallmark.IFN.DEGs$`11`
venn.list3<-list(`10 Control genes`=control_gene10,`11 Control genes`=control_gene11,`External dataset`=external_hallmark_interferon_in_DEGs)

venn.plot3 <- venn.diagram(
  x          = venn.list3,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Hallmark Type I IFN DEGs with External Dataset',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.58), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFB3C0","#92D5FF","#BFC0E3"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFB3C0","#92D5FF","#BFC0E3"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 1.4,#每個圓圈裡數字大小
  cat.dist   = c(0.075,0.075,0.03), #圓圈標題距離圓圈
  cat.default.pos='outer',
  margin     = 2, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot3)

#----{算出3個交集的}
Reduce(intersect,
       list(control_gene10,
            control_gene11,
            external_hallmark_interferon_in_DEGs)) #"LY6E"

####################【畫4圓交集圖】#################################
hallmark.IFN.DEGs<-read.csv("./RNA_DATA/Hallmark_IFN_DEG_genes_all.csv",check.names = FALSE)
control_gene1285<-hallmark.IFN.DEGs$`1285`
control_gene2268<-hallmark.IFN.DEGs$`2268`
control_gene3371<-hallmark.IFN.DEGs$`3371`
external_hallmark_interferon_in_DEGs<-c("ISG15", "IFI44L", "IFI44", "CMPK2", "RSAD2", "LAMP3", "RTP4", "CXCL10",
                                        "DDX60", "LY6E", "IFITM3", "BATF2", "IFIT2", "IFIT3", "OAS1", "OASL", "EPSTI1", "IFI27", "USP18", "MX1") #20

venn.list4<-list(`1285 Control genes`=control_gene1285,`2268 Control genes`=control_gene2268
                 ,`3371 Control genes`=control_gene3371,`External dataset`=external_hallmark_interferon_in_DEGs)

venn.plot4 <- venn.diagram(
  x          = venn.list4,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Hallmark Type I IFN DEGs with External Dataset',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.6), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#ffc18c","#D3C3FF","#8EFFE4","#BFC0E3"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#ffc18c","#D3C3FF","#8EFFE4","#BFC0E3"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 1.4,#每個圓圈裡數字大小
  #  cat.dist   = c(0.075,0.075,0.03), #圓圈標題距離圓圈
  cat.default.pos='outer',
  margin     = 2, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot4)

