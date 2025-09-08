setwd("C:/Users/Jane/Desktop/Wang實驗室/NMOSD研究計畫/RNA/NMOSD_RNA_analysis")
library(clusterProfiler)

#=====================================================<My dataset>=====================================================
######################【3 packages DEGs(criteria: padjust<0.05 & |log2FC|>1)】######################
#----【畫3圓交集圖】----
library(VennDiagram)
deseq2_degs<-read.csv('./RNA_DATA/NMOSD_RNA_DEGS(my data)(Deseq2).csv')$x
length(deseq2_degs) #3815

edgeR_degs<-read.csv('./RNA_DATA/NMOSD_RNA_DEGS(my data)(edgeR).csv')$x
length(edgeR_degs) #3811

limma_voom_degs<-read.csv('./RNA_DATA/NMOSD_RNA_DEGS(my data)(limma).csv')$x
length(limma_voom_degs) #3782

venn.list3<-list(`Limma-Voom`=limma_voom_degs,`Deseq2`=deseq2_degs
                 ,`edgeR`=edgeR_degs)
venn.plot3 <- venn.diagram(
  x          = venn.list3,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='3 protein coding DEGs intersection (adjp<0.05)',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.58), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#ffc18c","#D3C3FF","#8EFFE4"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#ffc18c","#D3C3FF","#8EFFE4"),
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

# ----取三者交集
all_intersect <- intersect(limma_voom_degs, intersect(deseq2_degs, edgeR_degs))
#write.csv(all_intersect,file="./RNA_DATA/3_protein_coding_DEGs_intersect(my data).csv", row.names = FALSE)
print(all_intersect)
length(all_intersect)  


# ----兩兩交集
LD <- intersect(limma_voom_degs, deseq2_degs)   # Limma ∩ Deseq2
LE <- intersect(limma_voom_degs, edgeR_degs)    # Limma ∩ edgeR
DE <- intersect(deseq2_degs,     edgeR_degs)    # DESeq2 ∩ edgeR

#合併兩兩交集，並去掉重複
degs_union_3packages <- unique(c(LD, LE, DE))

cat("兩兩交集（去重複）的基因數 =", length(degs_union_3packages), "\n")
#write.csv(degs_union_3packages ,file="./RNA_DATA/3_protein_coding_DEGs_two.two.intersect(my data).csv", row.names = FALSE)


######################【3 packages DEGs(criteria: padjust<0.01 & |log2FC|>2)】######################
#----【畫3圓交集圖】----
library(VennDiagram)
deseq2_degs<-read.csv('./RNA_DATA/NMOSD_RNA_DEGS_strictly(my data)(Deseq2).csv')$x
length(deseq2_degs) #1195

edgeR_degs<-read.csv('./RNA_DATA/NMOSD_RNA_DEGS_strictly(my data)(edgeR).csv')$x
length(edgeR_degs) #1188

limma_voom_degs<-read.csv('./RNA_DATA/NMOSD_RNA_DEGS_strictly(my data)(limma).csv')$x
length(limma_voom_degs) #1175

venn.list3<-list(`Limma-Voom`=limma_voom_degs,`Deseq2`=deseq2_degs
                 ,`edgeR`=edgeR_degs)
venn.plot3 <- venn.diagram(
  x          = venn.list3,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='3 protein coding DEGs intersection (adjp<0.01)',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.58), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#ffc18c","#D3C3FF","#8EFFE4"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#ffc18c","#D3C3FF","#8EFFE4"),
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

# 取三者交集
all_intersect <- intersect(limma_voom_degs, intersect(deseq2_degs, edgeR_degs))
#write.csv(all_intersect,file="./RNA_DATA/3_protein_coding_DEGs_intersect_strictly(my data).csv", row.names = FALSE)
print(all_intersect)
length(all_intersect) 


#=====================================================<External dataset>=====================================================
######################【3 packages DEGs(criteria: padjust<0.05 & |log2FC|>1)】######################
#----【畫3圓交集圖】----
library(VennDiagram)
deseq2_degs<-read.csv('./RNA_DATA/NMOSD_RNA_DEGS(external data)(Deseq2).csv')$x
length(deseq2_degs) #106

edgeR_degs<-read.csv('./RNA_DATA/NMOSD_RNA_DEGS(external data)(edgeR).csv')$x
length(edgeR_degs) #41

limma_voom_degs<-read.csv('./RNA_DATA/NMOSD_RNA_DEGS(external data)(limma).csv')$x
length(limma_voom_degs) #14

venn.list3<-list(`Limma-Voom`=limma_voom_degs,`Deseq2`=deseq2_degs
                 ,`edgeR`=edgeR_degs)
venn.plot3 <- venn.diagram(
  x          = venn.list3,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='External data 3 protein coding DEGs intersection (adjp<0.05)',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.52), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#ffc18c","#D3C3FF","#8EFFE4"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#ffc18c","#D3C3FF","#8EFFE4"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 1.4,#每個圓圈裡數字大小
  cat.dist   = c(0.09,0.09,0.05), #圓圈標題距離圓圈
  cat.default.pos='outer',
  margin     = 2, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot3)

# 取3者交集
all_intersect <- intersect(limma_voom_degs, intersect(deseq2_degs, edgeR_degs))
#write.csv(all_intersect,file="./RNA_DATA/3_protein_coding_DEGs_intersect(external data).csv", row.names = FALSE)
print(all_intersect)
length(all_intersect)  

# ----兩兩交集
LD <- intersect(limma_voom_degs, deseq2_degs)   # Limma ∩ Deseq2
LE <- intersect(limma_voom_degs, edgeR_degs)    # Limma ∩ edgeR
DE <- intersect(deseq2_degs,     edgeR_degs)    # DESeq2 ∩ edgeR

#合併兩兩交集，並去掉重複
degs_union_3packages <- unique(c(LD, LE, DE))

cat("兩兩交集（去重複）的基因數 =", length(degs_union_3packages), "\n")
#write.csv(degs_union_3packages ,file="./RNA_DATA/3_protein_coding_DEGs_two.two.intersect(external data).csv", row.names = FALSE)


######################【3 packages DEGs(criteria: padjust<0.01 & |log2FC|>2)】######################
#----【畫3圓交集圖】----
library(VennDiagram)
deseq2_degs<-read.csv('./RNA_DATA/NMOSD_RNA_DEGS_strictly(external data)(Deseq2).csv')$x
length(deseq2_degs) #17

edgeR_degs<-read.csv('./RNA_DATA/NMOSD_RNA_DEGS_strictly(external data)(edgeR).csv')$x
length(edgeR_degs) #7

limma_voom_degs<-list()
length(limma_voom_degs)

venn.list3<-list(`Limma-Voom`=limma_voom_degs,`Deseq2`=deseq2_degs
                 ,`edgeR`=edgeR_degs)
venn.plot3 <- venn.diagram(
  x          = venn.list3,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='External data 3 protein coding DEGs intersection (adjp<0.01)',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.58), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#ffc18c","#D3C3FF","#8EFFE4"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#ffc18c","#D3C3FF","#8EFFE4"),
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

# 取2者交集
two_intersect <- intersect(deseq2_degs, edgeR_degs)
#write.csv(all_intersect,file="./RNA_DATA/2_protein_coding_DEGs_intersect_strictly(external data).csv", row.names = FALSE)
print(two_intersect)
length(two_intersect)  # 聯集元素數量

#---------------------------------------------------------------------------------------------------------------------------
##############【Origin data & External data protein coding-DEGs(criteria: padjust<0.05 & |log2FC|>1) 】#############
#----{My & External three three intersect DEGs}
#-----------------【畫2圓交集圖】------------------
library(VennDiagram)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
my.intersect.DEGs<-read.csv("./RNA_DATA/3_protein_coding_DEGs_intersect(my data).csv",check.names = FALSE)$x
external.intersect.DEGs<-read.csv("./RNA_DATA/3_protein_coding_DEGs_intersect(external data).csv",check.names = FALSE)$x
venn.list2<-list(`My dataset`=my.intersect.DEGs,`External dataset`=external.intersect.DEGs)
my_external_protein.coding.DEGs<-intersect(my.intersect.DEGs,external.intersect.DEGs)
my_external_protein.coding.DEGs_df<-bitr(my_external_protein.coding.DEGs,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='My dataset & External dataset protein coding DEGs',
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


#----{My & External two two intersect DEGs}
#-----------------【畫2圓交集圖】------------------
library(VennDiagram)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
my.intersect.DEGs<-read.csv("./RNA_DATA/3_protein_coding_DEGs_two.two.intersect(my data).csv",check.names = FALSE)$x
external.intersect.DEGs<-read.csv("./RNA_DATA/3_protein_coding_DEGs_two.two.intersect(external data).csv",check.names = FALSE)$x
venn.list2<-list(`My dataset`=my.intersect.DEGs,`External dataset`=external.intersect.DEGs)
my_external_protein.coding.DEGs<-intersect(my.intersect.DEGs,external.intersect.DEGs)
my_external_protein.coding.DEGs_df<-bitr(my_external_protein.coding.DEGs,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='My dataset & External dataset protein coding DEGs',
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




###############################【Pathway analysis:ORA】###############################
#----{Package}----
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

#----{My & External two two intersect DEGs:8個}----
my_external_protein.coding.DEGs <- c(
  "ENSG00000187608",
  "ENSG00000126709",
  "ENSG00000134326",
  "ENSG00000131016",
  "ENSG00000160932",
  "ENSG00000135114",
  "ENSG00000167483",
  "ENSG00000159958"
)
my_external_protein.coding.DEGs_df<-bitr(my_external_protein.coding.DEGs,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')


#--------------【GO(MF)】------------------
#進行 enrichGO 富集分析
ego_MF<-enrichGO(gene=my_external_protein.coding.DEGs,
                 OrgDb=org.Hs.eg.db,
                 keyType = 'ENSEMBL',
                 ont='MF', #也可以是CC、BP
                 pAdjustMethod = 'BH', #Benjamini-Hochberg
                 pvalueCutoff = 0.05, #1為不過濾
                 qvalueCutoff = 1, #1為不過濾
                 readable = T) #將Gene ID轉成gene Symbol(易讀)
ego_MF_df<-ego_MF@result
View(ego_MF_df)

#dotplot
dotplot(ego_MF,title='EnrichmentGO_MF',label_format = 100)
#label_format讓文字不會重疊

#barplot
barplot(ego_MF,showCategory=20,title='EnrichmentGO_MF_bar',label_format = 100)
#showCategory=20繪製前20個

#過濾
ego_fil<-simplify(ego_MF,cutoff=0.5,by='pvalue',select_fun = match.fun("min"))
View(ego_fil@result)


#--------------【GO(BP)】------------------
#----{DEGs}
#進行 enrichGO 富集分析
ego_BP<-enrichGO(gene=my_external_protein.coding.DEGs,
                 OrgDb=org.Hs.eg.db,
                 keyType = 'ENSEMBL',
                 ont='BP', #也可以是CC、BP
                 pAdjustMethod = 'BH', #Benjamini-Hochberg
                 pvalueCutoff = 0.05, #adjusted pvalue cutoff #1為不過濾
                 readable = T) #將Gene ID轉成gene Symbol(易讀)
ego_BP_df<-ego_BP@result
View(ego_BP_df)
significant_BP_result<-ego_BP@result[ego_BP@result$p.adjust<0.05,]
nrow(significant_BP_result)#41
View(significant_BP_result)

#---【與 type I interferon 有關的 pathway】
typeI_IFN <- grepl("type I interferon|interferon-alpha|interferon-beta", significant_BP_result$Description, ignore.case = TRUE)

# 建立新的資料框
typeI_IFN_df <- significant_BP_result[typeI_IFN, ]
View(typeI_IFN_df)
nrow(typeI_IFN_df)#4

#dotplot
dotplot(ego_BP,title='Three DEanalysis protein coding DEGs Enrichment GO(BP)',font.size = 16,label_format = 100,showCategory = 41)
#label_format讓文字不會重疊

#barplot
barplot(ego_BP,showCategory=20,title='EnrichmentGO_BP_bar',label_format = 100)
#showCategory=20繪製前20個

#--------------【KEGG pathway富集分析】------------------
#----{DEGs}
kk<-enrichKEGG(gene=my_external_protein.coding.DEGs_df$ENTREZID,#你的基因列表
               organism='hsa', #指定物種 #hsa人類
               pvalueCutoff = 0.05 #1不進行過濾
)

View(kk@result)
nrow(kk@result[kk@result$p.adjust<0.05,])#0
View(kk@result[kk@result$p.adjust<0.05,]) 

#KEGG dotplot
dotplot(kk,title='Enrichment KEGG',font.size = 12,label_format = 100,showCategory = 30)

#------------------------------------------------------------------------------------------------------------------------


##################【Type.I.IFN-DEGs in 3 packages(criteria: padjust<0.05 & |log2FC|>1)】##################
#=====================================================<My dataset>=====================================================
#----【畫3圓交集圖】----
library(VennDiagram)
deseq2_degs<-read.csv('./RNA_DATA/Type_I_IFN_DEG_genes(my data)(Deseq2).csv')$x
length(deseq2_degs) #2102

edgeR_degs<-read.csv('./RNA_DATA/Type_I_IFN_DEG_genes(my data)(edgeR).csv')$x
length(edgeR_degs) #2105

limma_voom_degs<-read.csv('./RNA_DATA/Type_I_IFN_DEG_genes(my data)(limma).csv')$x
length(limma_voom_degs) #2090

venn.list3<-list(`Limma-Voom`=limma_voom_degs,`Deseq2`=deseq2_degs
                 ,`edgeR`=edgeR_degs)
venn.plot3 <- venn.diagram(
  x          = venn.list3,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='My data Type.I.IFN related-DEGs intersection (adjp<0.05)',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.58), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#ffc18c","#D3C3FF","#8EFFE4"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#ffc18c","#D3C3FF","#8EFFE4"),
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

# 取3者交集
all_intersect <- intersect(limma_voom_degs, intersect(deseq2_degs, edgeR_degs))
#write.csv(all_intersect,file="./RNA_DATA/3_Type.I.IFN.DEGs_intersect(my data).csv", row.names = FALSE)
print(all_intersect)
length(all_intersect)  

#=====================================================<External dataset>=====================================================
#----【畫3圓交集圖】----
library(VennDiagram)
deseq2_degs<-read.csv('./RNA_DATA/Type_I_IFN_DEG_genes(external data)(Deseq2).csv')$x
length(deseq2_degs) #56

edgeR_degs<-read.csv('./RNA_DATA/Type_I_IFN_DEG_genes(external data)(edgeR).csv')$x
length(edgeR_degs) #26

limma_voom_degs<-read.csv('./RNA_DATA/Type_I_IFN_DEG_genes(external data)(limma).csv')$x
length(limma_voom_degs) #10

venn.list3<-list(`Limma-Voom`=limma_voom_degs,`Deseq2`=deseq2_degs
                 ,`edgeR`=edgeR_degs)
venn.plot3 <- venn.diagram(
  x          = venn.list3,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='External data Type.I.IFN related-DEGs intersection (adjp<0.05)',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.58), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#ffc18c","#D3C3FF","#8EFFE4"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#ffc18c","#D3C3FF","#8EFFE4"),
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

# 取3者交集
all_intersect <- intersect(limma_voom_degs, intersect(deseq2_degs, edgeR_degs))
#write.csv(all_intersect,file="./RNA_DATA/3_Type.I.IFN.DEGs_intersect(external data).csv", row.names = FALSE)
print(all_intersect)
length(all_intersect)  

#======================================{My & External Type.I.IFN intersect DEGs}======================================
#-----------------【畫2圓交集圖】------------------
library(VennDiagram)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
my.intersect.type.I.IFN.DEGs<-read.csv("./RNA_DATA/3_Type.I.IFN.DEGs_intersect(my data).csv",check.names = FALSE)$x
external.intersect.type.I.IFN.DEGs<-read.csv("./RNA_DATA/3_Type.I.IFN.DEGs_intersect(external data).csv",check.names = FALSE)$x
venn.list2<-list(`My dataset`=my.intersect.type.I.IFN.DEGs,`External dataset`=external.intersect.type.I.IFN.DEGs)
my_external_type.I.IFN.DEGs<-intersect(my.intersect.type.I.IFN.DEGs,external.intersect.type.I.IFN.DEGs)
my_external_type.I.IFN.DEGs_df<-bitr(my_external_type.I.IFN.DEGs,fromType = 'SYMBOL',toType=c('ENTREZID','ENSEMBL'),OrgDb='org.Hs.eg.db')
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='My dataset & External dataset Type I IFN DEGs',
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



###################################【8 DEGs in IFN-I list】###########################
#----{Package}----
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

#----{My & External two two intersect DEGs:8個}----
my_external_protein.coding.DEGs <- c(
  "ENSG00000187608",
  "ENSG00000126709",
  "ENSG00000134326",
  "ENSG00000131016",
  "ENSG00000160932",
  "ENSG00000135114",
  "ENSG00000167483",
  "ENSG00000159958"
)
my_external_protein.coding.DEGs_df<-bitr(my_external_protein.coding.DEGs,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')

type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(my_external_protein.coding.DEGs_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
length(interferon_in_DEGs)
print(interferon_in_DEGs)
