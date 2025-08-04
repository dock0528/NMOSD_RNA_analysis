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

######################【3371 Housekeeping Genes】#######################
###【GO】
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(Batch Effect Correction).csv',header=T)

DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
#轉成ENTREZID & SYMBOL (從308DEGs -> 300DEGs)

#消失的ENSG (9個)
DEG_genes[!DEG_genes %in% DEG_genes_df$ENSEMBL]

#重複的ENSG(1個ENSG:2個ENTREZID:2個SYMBOL)
dup_ENSGs<-DEG_genes_df[duplicated(DEG_genes_df$ENSEMBL),]$ENSEMBL

#重複的ENTREZID(0個)
dup_ENTREZIDs<-DEG_genes_df[duplicated(DEG_genes_df$ENTREZID),]$ENTREZID
DEG_genes_df[DEG_genes_df$ENTREZID %in% dup_ENTREZIDs,]

#重複的SYMBOL(0個)
dup_SYMBOLs<-DEG_genes_df[duplicated(DEG_genes_df$SYMBOL),]$SYMBOL
DEG_genes_df[DEG_genes_df$SYMBOL %in% dup_SYMBOLs,]

###DEG ENTREZID資料集
DEG_genes<-DEG_genes_df$ENTREZID

#取和MsigDB有交集的type l interferon genes
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes) #和DEGs交集有165
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有1

######################【3371 Housekeeping Genes(STRICTLY)】#######################
###【GO】
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(Batch Effect Correction).csv',header=T)

DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
#轉成ENTREZID & SYMBOL (從54DEGs -> 55DEGs)

#消失的ENSG (0個)
DEG_genes[!DEG_genes %in% DEG_genes_df$ENSEMBL]

#重複的ENSG(1個ENSG:2個ENTREZID:2個SYMBOL)
dup_ENSGs<-DEG_genes_df[duplicated(DEG_genes_df$ENSEMBL),]$ENSEMBL

#重複的ENTREZID(0個)
dup_ENTREZIDs<-DEG_genes_df[duplicated(DEG_genes_df$ENTREZID),]$ENTREZID
DEG_genes_df[DEG_genes_df$ENTREZID %in% dup_ENTREZIDs,]

#重複的SYMBOL(0個)
dup_SYMBOLs<-DEG_genes_df[duplicated(DEG_genes_df$SYMBOL),]$SYMBOL
DEG_genes_df[DEG_genes_df$SYMBOL %in% dup_SYMBOLs,]

###DEG ENTREZID資料集
DEG_genes<-DEG_genes_df$ENTREZID

#取和MsigDB有交集的type l interferon genes
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes) #和DEGs交集有27
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0



######################【1285 Housekeeping Genes】#######################
###【GO】
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(1285Control Genes Batch Effect Correction).csv',header=T)

DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
#轉成ENTREZID & SYMBOL (從355DEGs -> 352DEGs)

#消失的ENSG (5個)
DEG_genes[!DEG_genes %in% DEG_genes_df$ENSEMBL]

#重複的ENSG(2個ENSG)
dup_ENSGs<-DEG_genes_df[duplicated(DEG_genes_df$ENSEMBL),]$ENSEMBL

#重複的ENTREZID(0個)
dup_ENTREZIDs<-DEG_genes_df[duplicated(DEG_genes_df$ENTREZID),]$ENTREZID
DEG_genes_df[DEG_genes_df$ENTREZID %in% dup_ENTREZIDs,]

#重複的SYMBOL(0個)
dup_SYMBOLs<-DEG_genes_df[duplicated(DEG_genes_df$SYMBOL),]$SYMBOL
DEG_genes_df[DEG_genes_df$SYMBOL %in% dup_SYMBOLs,]

###DEG ENTREZID資料集
DEG_genes<-DEG_genes_df$ENTREZID

#取和MsigDB有交集的type l interferon genes
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes) #和DEGs交集有178
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0

######################【1285 Housekeeping Genes (STRICTLY)】#######################
###【GO】
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(1285Control Genes Batch Effect Correction).csv',header=T)

DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
#轉成ENTREZID & SYMBOL (從37DEGs -> 37DEGs)

#消失的ENSG (0個)
DEG_genes[!DEG_genes %in% DEG_genes_df$ENSEMBL]

#重複的ENSG(0個ENSG)
dup_ENSGs<-DEG_genes_df[duplicated(DEG_genes_df$ENSEMBL),]$ENSEMBL

#重複的ENTREZID(0個)
dup_ENTREZIDs<-DEG_genes_df[duplicated(DEG_genes_df$ENTREZID),]$ENTREZID
DEG_genes_df[DEG_genes_df$ENTREZID %in% dup_ENTREZIDs,]

#重複的SYMBOL(0個)
dup_SYMBOLs<-DEG_genes_df[duplicated(DEG_genes_df$SYMBOL),]$SYMBOL
DEG_genes_df[DEG_genes_df$SYMBOL %in% dup_SYMBOLs,]

###DEG ENTREZID資料集
DEG_genes<-DEG_genes_df$ENTREZID

#取和MsigDB有交集的type l interferon genes
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes) #和DEGs交集有19
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0


#########################################【Function 計算IFN數量】#########################################
#---------{Origin DEGs (criteria: padjust<0.05 & |log2FC|>1}---------

#---{Interferon count}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_origin(limma).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE)
interferon_df<-data.frame(Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
print(interferon_df)

#---{type.I.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_origin(limma).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
print(interferon_in_DEGs)
#write.csv(interferon_in_DEGs,file = "./RNA_DATA/Type_I_IFN_DEG_genes(origin data).csv",row.names = FALSE)

#---{hallmark.IFN.list}
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_origin(limma).csv',header=T)
DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
print(hallmark_interferon_in_DEGs)
#write.csv(hallmark_interferon_in_DEGs,file= "./RNA_DATA/Hallmark_IFN_DEG_genes(origin data).csv",row.names = FALSE)

#---------{Housekeeping genes (criteria: padjust<0.05 & |log2FC|>1}---------
interferon_count<-function(number){
  path<-paste0('RNA_DATA/NMOSD_RNA_DEGS(',number,'Control Genes Batch Effect Correction).csv')
  DEGs_df<-read.csv(path,header=T)
  DEG_genes<-DEGs_df$x
  DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
  type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
  interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
  hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
  hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
  interferon_df<-data.frame(Control.genes=number,Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
  return(interferon_df)
}

type.I.IFN.list<-function(number){
  path<-paste0('RNA_DATA/NMOSD_RNA_DEGS(',number,'Control Genes Batch Effect Correction).csv')
  DEGs_df<-read.csv(path,header=T)
  DEG_genes<-DEGs_df$x
  DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
  type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
  interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
  return(interferon_in_DEGs)
}

hallmark.IFN.list<-function(number){
  path<-paste0('RNA_DATA/NMOSD_RNA_DEGS(',number,'Control Genes Batch Effect Correction).csv')
  DEGs_df<-read.csv(path,header=T)
  DEG_genes<-DEGs_df$x
  DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
  hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
  hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
  return(hallmark_interferon_in_DEGs)
}
interferon_df<-data.frame()
Type.I.IFN_df<-data.frame(matrix(nrow = 1, ncol = 0))
Hallmark.IFN_df<-data.frame(matrix(nrow = 1, ncol = 0))
number<-c(10,11,1285,2268,3371)
for(num in number){
  interferon_df<-rbind(interferon_df,interferon_count(num))
  Type.I.IFN_df[[ as.character(num) ]] <- list(type.I.IFN.list(num))
  Hallmark.IFN_df[[ as.character(num) ]] <- list(hallmark.IFN.list(num))
}

print(interferon_df)
str(Type.I.IFN_df) #str()結構
str(Hallmark.IFN_df) 
#Type.I.IFN_df[["10"]][[1]] 取control genes=10的交集Type.I.IFN genes

#----{Type.I.IFN_df轉成可存的df}
#取每個control genes下的Type I IFN genes
list_cols <- lapply(Type.I.IFN_df, `[[`, 1)

# 找最長col
max_len <- max(lengths(list_cols))

#把短向量尾端填 NA
padded <- lapply(list_cols, function(v) {
  length(v) <- max_len #修改長度
  v #return修改後的v
})

#再轉成df
Type.I.IFN_saved_df<- as.data.frame(padded,
                         check.names = FALSE,
                         stringsAsFactors = FALSE)
#write.csv(Type.I.IFN_saved_df,file = "./RNA_DATA/Type_I_IFN_DEG_genes_all.csv",row.names = FALSE)
#-------------------------------------------

#----{Hallmark.IFN_df轉成可存的df}
#取每個control genes下的Type I IFN genes
list_cols <- lapply(Hallmark.IFN_df, `[[`, 1)

# 找最長col
max_len <- max(lengths(list_cols))

#把短向量尾端填 NA
padded <- lapply(list_cols, function(v) {
  length(v) <- max_len #修改長度
  v #return修改後的v
})

#再轉成df
Hallmark.IFN_df_saved_df<- as.data.frame(padded,
                                    check.names = FALSE,
                                    stringsAsFactors = FALSE)
#write.csv(Hallmark.IFN_df_saved_df,file= "./RNA_DATA/Hallmark_IFN_DEG_genes_all.csv",row.names = FALSE)

#---------------【Type I IFN related-DEGs】------------------
####################【畫2圓交集圖】#################################
install.packages("VennDiagram")
library(VennDiagram)
type.I.IFN.DEGs<-read.csv("./RNA_DATA/Type_I_IFN_DEG_genes_all.csv",check.names = FALSE)
control_gene10<-type.I.IFN.DEGs$`10`
control_gene11<-type.I.IFN.DEGs$`11`
venn.list2<-list(`10 Control genes`=control_gene10,`11 Control genes`=control_gene11)
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Type I IFN related-DEGs',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.61), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFB3C0","#92D5FF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFB3C0","#92D5FF"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 1.5,#每個圓圈裡數字大小
  cat.dist   = 0.05, #圓圈標題距離圓圈
  cat.pos    = c(30,-30), #圓圈標題角度 #往右正
  margin     = 1, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot2)
####################【畫3圓交集圖】#################################
library(VennDiagram)
type.I.IFN.DEGs<-read.csv("./RNA_DATA/Type_I_IFN_DEG_genes_all.csv",check.names = FALSE)
control_gene1285<-type.I.IFN.DEGs$`1285`
control_gene2268<-type.I.IFN.DEGs$`2268`
control_gene3371<-type.I.IFN.DEGs$`3371`
venn.list3<-list(`1285 Control genes`=control_gene1285,`2268 Control genes`=control_gene2268
                 ,`3371 Control genes`=control_gene3371)
venn.plot3 <- venn.diagram(
  x          = venn.list3,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Type I IFN related-DEGs',
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
####################【畫5圓交集圖】#################################
library(VennDiagram)
type.I.IFN.DEGs<-read.csv("./RNA_DATA/Type_I_IFN_DEG_genes_all.csv",check.names = FALSE)
control_gene10<-type.I.IFN.DEGs$`10`
control_gene11<-type.I.IFN.DEGs$`11`
control_gene1285<-type.I.IFN.DEGs$`1285`
control_gene2268<-type.I.IFN.DEGs$`2268`
control_gene3371<-type.I.IFN.DEGs$`3371`
venn.list5<-list(`10 Control genes`=control_gene10,`11 Control genes`=control_gene11,
                 `1285 Control genes`=control_gene1285,`2268 Control genes`=control_gene2268
                ,`3371 Control genes`=control_gene3371)
venn.plot5 <- venn.diagram(
  x          = venn.list5,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Type I IFN related-DEGs',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.64), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFB3C0","#92D5FF","#ffc18c","#D3C3FF","#8EFFE4"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFB3C0","#92D5FF","#ffc18c","#D3C3FF","#8EFFE4"),
  alpha      = 0.4, 
  cat.cex    = 1.4, #每個圓圈的標題大小
  cex        = 1.2,#每個圓圈裡數字大小
  cat.dist   = c(0.2,0.2,0.3,0.25,0.2), #圓圈標題距離圓圈
  cat.pos=c(0,-10,250,130,10),
  cat.default.pos='outer',
  margin     = 2, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot5)
#----------------------------------------------------------------


#---------------【Hallmark Type I IFN DEGs】------------------
####################【畫2圓交集圖】#################################
install.packages("VennDiagram")
library(VennDiagram)
hallmark.IFN.DEGs<-read.csv("./RNA_DATA/Hallmark_IFN_DEG_genes_all.csv",check.names = FALSE)
control_gene10<-hallmark.IFN.DEGs$`10`
control_gene11<-hallmark.IFN.DEGs$`11`
venn.list2<-list(`10 Control genes`=control_gene10,`11 Control genes`=control_gene11)
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Hallmark Type I IFN DEGs',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.61), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFB3C0","#92D5FF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFB3C0","#92D5FF"),
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
library(VennDiagram)
hallmark.IFN.DEGs<-read.csv("./RNA_DATA/Hallmark_IFN_DEG_genes_all.csv",check.names = FALSE)
control_gene1285<-hallmark.IFN.DEGs$`1285`
control_gene2268<-hallmark.IFN.DEGs$`2268`
control_gene3371<-hallmark.IFN.DEGs$`3371`
venn.list3<-list(`1285 Control genes`=control_gene1285,`2268 Control genes`=control_gene2268
                 ,`3371 Control genes`=control_gene3371)
venn.plot3 <- venn.diagram(
  x          = venn.list3,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Hallmark Type I IFN DEGs',
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
####################【畫5圓交集圖】#################################
library(VennDiagram)
hallmark.IFN.DEGs<-read.csv("./RNA_DATA/Hallmark_IFN_DEG_genes_all.csv",check.names = FALSE)
control_gene10<-hallmark.IFN.DEGs$`10`
control_gene11<-hallmark.IFN.DEGs$`11`
control_gene1285<-hallmark.IFN.DEGs$`1285`
control_gene2268<-hallmark.IFN.DEGs$`2268`
control_gene3371<-hallmark.IFN.DEGs$`3371`
venn.list5<-list(`10 Control genes`=control_gene10,`11 Control genes`=control_gene11,
                 `1285 Control genes`=control_gene1285,`2268 Control genes`=control_gene2268
                 ,`3371 Control genes`=control_gene3371)
venn.plot5 <- venn.diagram(
  x          = venn.list5,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Hallmark Type I IFN DEGs',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.64), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFB3C0","#92D5FF","#ffc18c","#D3C3FF","#8EFFE4"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFB3C0","#92D5FF","#ffc18c","#D3C3FF","#8EFFE4"),
  alpha      = 0.4, 
  cat.cex    = 1.4, #每個圓圈的標題大小
  cex        = 1.2,#每個圓圈裡數字大小
  cat.dist   = c(0.2,0.2,0.3,0.25,0.2), #圓圈標題距離圓圈
  cat.pos=c(0,-10,250,130,10),
  cat.default.pos='outer',
  margin     = 2, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot5)
#----------------------------------------------------------------



#---------{Housekeeping genes (criteria: padjust<0.01 & |log2FC|>2}---------
interferon_strictly_count<-function(number){
  path<-paste0('RNA_DATA/NMOSD_RNA_DEGS_strictly(',number,'Control Genes Batch Effect Correction).csv')
  DEGs_df<-read.csv(path,header=T)
  DEG_genes<-DEGs_df$x
  DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
  type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
  interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
  hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
  hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
  interferon_df<-data.frame(Control.genes=number,Type.I.IFN=length(interferon_in_DEGs),Hallmark.IFN=length(hallmark_interferon_in_DEGs))
  return(interferon_df)
}
type.I.IFN.list<-function(number){
  path<-paste0('RNA_DATA/NMOSD_RNA_DEGS_strictly(',number,'Control Genes Batch Effect Correction).csv')
  DEGs_df<-read.csv(path,header=T)
  DEG_genes<-DEGs_df$x
  DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
  type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
  interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes)
  return(interferon_in_DEGs)
}

hallmark.IFN.list<-function(number){
  path<-paste0('RNA_DATA/NMOSD_RNA_DEGS_strictly(',number,'Control Genes Batch Effect Correction).csv')
  DEGs_df<-read.csv(path,header=T)
  DEG_genes<-DEGs_df$x
  DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
  hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
  hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有0
  return(hallmark_interferon_in_DEGs)
}
interferon_df<-data.frame()
Type.I.IFN_df<-data.frame(matrix(nrow = 1, ncol = 0))
Hallmark.IFN_df<-data.frame(matrix(nrow = 1, ncol = 0))
number<-c(10,11,1285,2268,3371)
for(num in number){
  interferon_df<-rbind(interferon_df,interferon_strictly_count(num))
  Type.I.IFN_df[[ as.character(num) ]] <- list(type.I.IFN.list(num))
  Hallmark.IFN_df[[ as.character(num) ]] <- list(hallmark.IFN.list(num))
}

print(interferon_df)
str(Type.I.IFN_df) #str()結構
str(Hallmark.IFN_df) 
#Type.I.IFN_df[["10"]][[1]] 取control genes=10的交集Type.I.IFN genes

#----{Type.I.IFN_df轉成可存的df}
#取每個control genes下的Type I IFN genes
list_cols <- lapply(Type.I.IFN_df, `[[`, 1)

# 找最長col
max_len <- max(lengths(list_cols))

#把短向量尾端填 NA
padded <- lapply(list_cols, function(v) {
  length(v) <- max_len #修改長度
  v #return修改後的v
})

#再轉成df
Type.I.IFN_saved_df<- as.data.frame(padded,
                                    check.names = FALSE,
                                    stringsAsFactors = FALSE)
#write.csv(Type.I.IFN_saved_df,file = "./RNA_DATA/Type_I_IFN_DEG_genes_all_strictly.csv",row.names = FALSE)
#-------------------------------------------------
#----{Hallmark.IFN_df轉成可存的df}
#取每個control genes下的Type I IFN genes
list_cols <- lapply(Hallmark.IFN_df, `[[`, 1)

# 找最長col
max_len <- max(lengths(list_cols))

#把短向量尾端填 NA
padded <- lapply(list_cols, function(v) {
  length(v) <- max_len #修改長度
  v #return修改後的v
})

#再轉成df
Hallmark.IFN_df_saved_df<- as.data.frame(padded,
                                         check.names = FALSE,
                                         stringsAsFactors = FALSE)
#write.csv(Hallmark.IFN_df_saved_df,file= "./RNA_DATA/Hallmark_IFN_DEG_genes_all_strictly.csv",row.names = FALSE)

#---------------【Type I IFN related-DEGs】------------------
####################【畫2圓交集圖】#################################
install.packages("VennDiagram")
library(VennDiagram)
type.I.IFN.DEGs<-read.csv("./RNA_DATA/Type_I_IFN_DEG_genes_all_strictly.csv",check.names = FALSE)
control_gene10<-type.I.IFN.DEGs$`10`
control_gene11<-type.I.IFN.DEGs$`11`
venn.list2<-list(`10 Control genes`=control_gene10,`11 Control genes`=control_gene11)
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Type I IFN related-DEGs',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.61), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFB3C0","#92D5FF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFB3C0","#92D5FF"),
  alpha      = 0.5,
  cat.cex    = 1.6, #每個圓圈的標題大小
  cex        = 1.5,#每個圓圈裡數字大小
  cat.dist   = 0.05, #圓圈標題距離圓圈
  cat.pos    = c(30,-30), #圓圈標題角度 #往右正
  margin     = 1, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot2)
####################【畫3圓交集圖】#################################
library(VennDiagram)
type.I.IFN.DEGs<-read.csv("./RNA_DATA/Type_I_IFN_DEG_genes_all_strictly.csv",check.names = FALSE)
control_gene1285<-type.I.IFN.DEGs$`1285`
control_gene2268<-type.I.IFN.DEGs$`2268`
control_gene3371<-type.I.IFN.DEGs$`3371`
venn.list3<-list(`1285 Control genes`=control_gene1285,`2268 Control genes`=control_gene2268
                 ,`3371 Control genes`=control_gene3371)
venn.plot3 <- venn.diagram(
  x          = venn.list3,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Type I IFN related-DEGs',
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
####################【畫5圓交集圖】#################################
library(VennDiagram)
type.I.IFN.DEGs<-read.csv("./RNA_DATA/Type_I_IFN_DEG_genes_all_strictly.csv",check.names = FALSE)
control_gene10<-type.I.IFN.DEGs$`10`
control_gene11<-type.I.IFN.DEGs$`11`
control_gene1285<-type.I.IFN.DEGs$`1285`
control_gene2268<-type.I.IFN.DEGs$`2268`
control_gene3371<-type.I.IFN.DEGs$`3371`
venn.list5<-list(`10 Control genes`=control_gene10,`11 Control genes`=control_gene11,
                 `1285 Control genes`=control_gene1285,`2268 Control genes`=control_gene2268
                 ,`3371 Control genes`=control_gene3371)
venn.plot5 <- venn.diagram(
  x          = venn.list5,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Type I IFN related-DEGs',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.64), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFB3C0","#92D5FF","#ffc18c","#D3C3FF","#8EFFE4"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFB3C0","#92D5FF","#ffc18c","#D3C3FF","#8EFFE4"),
  alpha      = 0.4, 
  cat.cex    = 1.4, #每個圓圈的標題大小
  cex        = 1.2,#每個圓圈裡數字大小
  cat.dist   = c(0.2,0.2,0.3,0.25,0.2), #圓圈標題距離圓圈
  cat.pos=c(0,-10,250,130,10),
  cat.default.pos='outer',
  margin     = 2, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot5)
#----------------------------------------------------------------


#---------------【Hallmark Type I IFN DEGs】------------------
####################【畫2圓交集圖】#################################
install.packages("VennDiagram")
library(VennDiagram)
hallmark.IFN.DEGs<-read.csv("./RNA_DATA/Hallmark_IFN_DEG_genes_all_strictly.csv",check.names = FALSE)
control_gene10<-hallmark.IFN.DEGs$`10`
control_gene11<-hallmark.IFN.DEGs$`11`
venn.list2<-list(`10 Control genes`=control_gene10,`11 Control genes`=control_gene11)
venn.plot2 <- venn.diagram(
  x          = venn.list2,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Hallmark Type I IFN DEGs',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.61), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFB3C0","#92D5FF"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFB3C0","#92D5FF"),
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
library(VennDiagram)
hallmark.IFN.DEGs<-read.csv("./RNA_DATA/Hallmark_IFN_DEG_genes_all_strictly.csv",check.names = FALSE)
control_gene1285<-hallmark.IFN.DEGs$`1285`
control_gene2268<-hallmark.IFN.DEGs$`2268`
control_gene3371<-hallmark.IFN.DEGs$`3371`
venn.list3<-list(`1285 Control genes`=control_gene1285,`2268 Control genes`=control_gene2268
                 ,`3371 Control genes`=control_gene3371)
venn.plot3 <- venn.diagram(
  x          = venn.list3,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Hallmark Type I IFN DEGs',
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
####################【畫5圓交集圖】#################################
library(VennDiagram)
hallmark.IFN.DEGs<-read.csv("./RNA_DATA/Hallmark_IFN_DEG_genes_all_strictly.csv",check.names = FALSE)
control_gene10<-hallmark.IFN.DEGs$`10`
control_gene11<-hallmark.IFN.DEGs$`11`
control_gene1285<-hallmark.IFN.DEGs$`1285`
control_gene2268<-hallmark.IFN.DEGs$`2268`
control_gene3371<-hallmark.IFN.DEGs$`3371`
venn.list5<-list(`10 Control genes`=control_gene10,`11 Control genes`=control_gene11,
                 `1285 Control genes`=control_gene1285,`2268 Control genes`=control_gene2268
                 ,`3371 Control genes`=control_gene3371)
venn.plot5 <- venn.diagram(
  x          = venn.list5,
  na         ='remove',
  scale      =F, #圓圈大小不依照數字大小變化
  main       ='Hallmark Type I IFN DEGs',
  main.cex   =2,
  main.fontface=2, #粗
  main.pos=c(0.5,0.64), #(x.y) #x->0:左 1:右 #y->0:下 1:上
  lwd=2,lty=1,col=c("#FFB3C0","#92D5FF","#ffc18c","#D3C3FF","#8EFFE4"), #lwd:線寬度 #lyt:1實線 2虛線 #col邊框顏色
  filename   = NULL,                #filename輸出圖片的名稱
  fill       = c("#FFB3C0","#92D5FF","#ffc18c","#D3C3FF","#8EFFE4"),
  alpha      = 0.4, 
  cat.cex    = 1.4, #每個圓圈的標題大小
  cex        = 1.2,#每個圓圈裡數字大小
  cat.dist   = c(0.2,0.2,0.3,0.25,0.2), #圓圈標題距離圓圈
  cat.pos=c(0,-10,250,130,10),
  cat.default.pos='outer',
  margin     = 2, #空白區域(越大越多)
)

# 顯示
grid.newpage()
grid.draw(venn.plot5)
#----------------------------------------------------------------


###########################################################################################################################

######################【10 Housekeeping Genes】#######################
###【GO】
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(10Control Genes Batch Effect Correction).csv',header=T)

DEG_genes<-DEGs_df$x
DEG_genes_df<-bitr(DEG_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
#轉成ENTREZID & SYMBOL (從355DEGs -> 352DEGs)

#消失的ENSG (5個)
DEG_genes[!DEG_genes %in% DEG_genes_df$ENSEMBL]

#重複的ENSG(2個ENSG)
dup_ENSGs<-DEG_genes_df[duplicated(DEG_genes_df$ENSEMBL),]$ENSEMBL

#重複的ENTREZID(0個)
dup_ENTREZIDs<-DEG_genes_df[duplicated(DEG_genes_df$ENTREZID),]$ENTREZID
DEG_genes_df[DEG_genes_df$ENTREZID %in% dup_ENTREZIDs,]

#重複的SYMBOL(0個)
dup_SYMBOLs<-DEG_genes_df[duplicated(DEG_genes_df$SYMBOL),]$SYMBOL

DEG_genes_df[DEG_genes_df$SYMBOL %in% dup_SYMBOLs,]

###DEG ENTREZID資料集
DEG_genes<-DEG_genes_df$ENTREZID

#取和MsigDB有交集的type l interferon genes
type_I_interferons<-read.csv('RNA_DATA/Type I Interferon gene list(MsigDB).csv',header=T) #Type I Interferon gene list:7935個
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_I_interferons$Type.I.Interferon.genes) #和DEGs交集有1618
hallmark_interferons<-read.csv('RNA_DATA/Interferon_hallmark.csv',header=T) #Hallmark Interferon gene list:97個
hallmark_interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,hallmark_interferons$HALLMARK_INTERFERON_ALPHA_RESPONSE) #和DEGs交集有24

reference_interferon_genes <- c(
  "IL1RN", "IFIT5", "IFIT3", "IFIT2", "IFI6", "HERC5",
  "CMPK2", "IFIT1", "MX1", "ISG15", "OASL", "RSAD2",
  "IFI44L", "OAS3", "IFI44", "XAF1", "EIF2AK2", "LAP3",
  "LY6E", "OAS2", "OAS1", "USP18", "HERC6", "SIGLEC1",
  "IFI27") #原始3個csv交集出共同的Type I Interferon
interferon_with_reference<-intersect(interferon_in_DEGs,reference_interferon_genes)#"IL1RN"、"LY6E" 
hallmark_interferon_with_reference<-intersect(hallmark_interferon_in_DEGs,reference_interferon_genes) #"LY6E"

ego_BP<-enrichGO(gene=interferon_in_DEGs,
                 OrgDb=org.Hs.eg.db,
                 keyType = 'SYMBOL',
                 ont='BP', #也可以是CC、BP
                 pAdjustMethod = 'BH', #Benjamini-Hochberg
                 pvalueCutoff = 0.05, #1為不過濾
                 qvalueCutoff = 1, #1為不過濾
                 readable = T) #將Gene ID轉成gene Symbol(易讀)
ego_BP_df<-ego_BP@result
View(ego_BP_df)

#dotplot
dotplot(ego_BP,title='10 Housekeeping genes Type I IFN-DEGs BP(GO)',label_format = 100)
#label_format讓文字不會重疊

#barplot
barplot(ego_BP,showCategory=20,title='Type I IFN BP(GO)',label_format = 100)
#showCategory=20繪製前20個

#---【與 type I interferon 有關的 pathway】
typeI_IFN <- grepl("type I interferon|interferon-alpha|interferon-beta", ego_BP@result$Description, ignore.case = TRUE)

# 建立新的資料框
typeI_IFN_df <- ego_BP@result[typeI_IFN, ]
View(typeI_IFN_df )

#---{Type I IFN pathway genes取聯集}
type_I_IFN_pathway_genes <- strsplit(typeI_IFN_df$geneID, split = "/")
all_genes <- unlist(type_I_IFN_pathway_genes, use.names = FALSE)
unique_type_I_IFN_genes <- unique(all_genes)

IFN_ego_BP<-enrichGO(gene=unique_type_I_IFN_genes,
                 OrgDb=org.Hs.eg.db,
                 keyType = 'SYMBOL',
                 ont='BP', #也可以是CC、BP
                 pAdjustMethod = 'BH', #Benjamini-Hochberg
                 pvalueCutoff = 1, #1為不過濾
                 qvalueCutoff = 1, #1為不過濾
                 readable = T) #將Gene ID轉成gene Symbol(易讀)
IFN_ego_BP_df<-IFN_ego_BP@result
View(IFN_ego_BP_df)

#dotplot
dotplot(IFN_ego_BP,showCategory=10,title='10 Housekeeping genes Type I IFN-Signaling DEGs(GO-BP)',label_format = 100)





# 看一下結果
length(unique_type_I_IFN_genes)
head(unique_type_I_IFN_genes )

reference_interferon_genes <- c(
  "IL1RN", "IFIT5", "IFIT3", "IFIT2", "IFI6", "HERC5",
  "CMPK2", "IFIT1", "MX1", "ISG15", "OASL", "RSAD2",
  "IFI44L", "OAS3", "IFI44", "XAF1", "EIF2AK2", "LAP3",
  "LY6E", "OAS2", "OAS1", "USP18", "HERC6", "SIGLEC1",
  "IFI27") #原始3個csv交集出共同的Type I Interferon
interferon_with_reference<-intersect(unique_type_I_IFN_genes,reference_interferon_genes) #0個

# 檢視結果（這才是資料表）
View(typeI_IFN_df)


#-------------------------------------------------------------------------------------------------------------------------------

reference_interferon_genes <- c(
  "IL1RN", "IFIT5", "IFIT3", "IFIT2", "IFI6", "HERC5",
  "CMPK2", "IFIT1", "MX1", "ISG15", "OASL", "RSAD2",
  "IFI44L", "OAS3", "IFI44", "XAF1", "EIF2AK2", "LAP3",
  "LY6E", "OAS2", "OAS1", "USP18", "HERC6", "SIGLEC1",
  "IFI27") #原始3個csv交集出共同的Type I Interferon
interferon_with_reference<-intersect(interferon_in_DEGs,reference_interferon_genes)
hallmark_interferon_with_reference<-intersect(hallmark_interferon_in_DEGs,reference_interferon_genes)


#----Refernce:Other_therapy vs healthy DEGs
reference_other_therapy_IFNs<-read.csv('RNA_DATA/Other_therapy_healthy_DEGs_v2.csv',header=T)
interferon_in_other_therapy<-intersect(reference_other_therapy_IFNs$Gene,type_I_interferons$Type.I.Interferon.genes)#62
interferon_with_other_therapy<-intersect(interferon_in_DEGs,interferon_in_other_therapy)#無交集

#----Refernce:Other_therapy vs healthy DEGs
reference_untreat_IFNs<-read.csv('RNA_DATA/Untreat_healthy_DEGs_v2.csv',header=T)
interferon_in_untreat<-intersect(reference_untreat_IFNs$Gene,type_I_interferons$Type.I.Interferon.genes)#51
interferon_with_untreat<-intersect(interferon_in_DEGs,interferon_in_untreat)#無交集


#-----------{Strictly DEGs}----------------
strict_DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(11Control Genes Batch Effect Correction).csv',header=T)
strict_DEGs_genes<-strict_DEGs_df$x
strict_DEGs_genes_df<-bitr(strict_DEGs_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')


#-----------{DEGs}----------------
DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS(11Control Genes Batch Effect Correction).csv',header=T)
DEGs_genes<-DEGs_df$x
DEGs_genes_df<-bitr(DEGs_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')


#--------------【GO(MF)】------------------
#進行 enrichGO 富集分析
ego_MF<-enrichGO(gene=DEGs_genes,
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
ego_BP<-enrichGO(gene=DEGs_genes,
                 OrgDb=org.Hs.eg.db,
                 keyType = 'ENSEMBL',
                 ont='BP', #也可以是CC、BP
                 pAdjustMethod = 'BH', #Benjamini-Hochberg
                 pvalueCutoff = 0.01, #adjusted pvalue cutoff #1為不過濾
                 readable = T) #將Gene ID轉成gene Symbol(易讀)
ego_BP_df<-ego_BP@result
View(ego_BP_df)
significant_BP_result<-ego_BP@result[ego_BP@result$p.adjust<0.01,]
nrow(significant_BP_result)#253
View(significant_BP_result)

#---【與 type I interferon 有關的 pathway】
typeI_IFN <- grepl("type I interferon|interferon-alpha|interferon-beta", significant_BP_result$Description, ignore.case = TRUE)

# 建立新的資料框
typeI_IFN_df <- significant_BP_result[typeI_IFN, ]
View(typeI_IFN_df)
nrow(typeI_IFN_df)#1

#dotplot
dotplot(ego_BP,title='EnrichmentGO_BP',label_format = 100,showCategory = 30)
#label_format讓文字不會重疊

#barplot
barplot(ego_BP,showCategory=20,title='EnrichmentGO_BP_bar',label_format = 100)
#showCategory=20繪製前20個

#----{Strictly DEGs}
#進行 enrichGO 富集分析
ego_BP<-enrichGO(gene=strict_DEGs_genes,
                 OrgDb=org.Hs.eg.db,
                 keyType = 'ENSEMBL',
                 ont='BP', #也可以是CC、BP
                 pAdjustMethod = 'BH', #Benjamini-Hochberg
                 pvalueCutoff = 0.05, #adjusted pvalue cutoff #1為不過濾
                 readable = T) #將Gene ID轉成gene Symbol(易讀)
ego_BP_df<-ego_BP@result
View(ego_BP_df)
significant_BP_result<-ego_BP@result[ego_BP@result$p.adjust<0.05,]
View(significant_BP_result) #537

#---【與 type I interferon 有關的 pathway】
typeI_IFN <- grepl("type I interferon|interferon-alpha|interferon-beta", significant_BP_result$Description, ignore.case = TRUE)

# 建立新的資料框
typeI_IFN_df <- significant_BP_result[typeI_IFN, ]
View(typeI_IFN_df)

#dotplot
dotplot(ego_BP,title='EnrichmentGO_BP',label_format = 100,showCategory = 50)
#label_format讓文字不會重疊

#barplot
barplot(ego_BP,showCategory=20,title='EnrichmentGO_BP_bar',label_format = 100)
#showCategory=20繪製前20個

#--------------【KEGG pathway富集分析】------------------
kk<-enrichKEGG(gene=DEGs_genes_df$ENTREZID,#你的基因列表
               organism='hsa', #指定物種 #hsa人類
               pvalueCutoff = 0.01 #1不進行過濾
              )

View(kk@result)
nrow(kk@result[kk@result$p.adjust<0.01,])#38
View(kk@result[kk@result$p.adjust<0.01,]) 

#KEGG dotplot
dotplot(kk,title='Enrichment KEGG',font.size = 14,label_format = 100,showCategory = 30)

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

##############【Type I IFN genes pathway analysis】#################
IFN_related_genes_list<-read.csv("./RNA_DATA/Type I Interferon gene list(MsigDB).csv")
IFN_related_genes<-IFN_related_genes_list$Type.I.Interferon.genes

ego_BP<-enrichGO(gene=IFN_related_genes,
                 OrgDb=org.Hs.eg.db,
                 keyType = 'SYMBOL',
                 ont='BP', #也可以是CC、BP
                 pAdjustMethod = 'BH', #Benjamini-Hochberg
                 pvalueCutoff = 1, #1為不過濾
                 qvalueCutoff = 1, #1為不過濾
                 readable = T) #將Gene ID轉成gene Symbol(易讀)
ego_BP_df<-ego_BP@result
View(ego_BP_df)

#dotplot
dotplot(ego_BP,showCategory=10,title='Type I IFN BP(GO)',label_format = 100)
#label_format讓文字不會重疊

#barplot
barplot(ego_BP,showCategory=20,title='Type I IFN BP(GO)',label_format = 100)
#showCategory=20繪製前20個

#---【與 type I interferon 有關的 pathway】
typeI_IFN <- grepl("type I interferon|interferon-alpha|interferon-beta", ego_BP@result$Description, ignore.case = TRUE)

# 建立新的資料框
typeI_IFN_df <- ego_BP@result[typeI_IFN, ]

#---{Type I IFN pathway genes取聯集}
type_I_IFN_pathway_genes <- strsplit(typeI_IFN_df$geneID, split = "/")
all_genes <- unlist(type_I_IFN_pathway_genes, use.names = FALSE)
unique_type_I_IFN_genes <- unique(all_genes)

# 看一下結果
length(unique_type_I_IFN_genes)
head(unique_type_I_IFN_genes )

reference_interferon_genes <- c(
  "IL1RN", "IFIT5", "IFIT3", "IFIT2", "IFI6", "HERC5",
  "CMPK2", "IFIT1", "MX1", "ISG15", "OASL", "RSAD2",
  "IFI44L", "OAS3", "IFI44", "XAF1", "EIF2AK2", "LAP3",
  "LY6E", "OAS2", "OAS1", "USP18", "HERC6", "SIGLEC1",
  "IFI27") #原始3個csv交集出共同的Type I Interferon
interferon_with_reference<-intersect(unique_type_I_IFN_genes,reference_interferon_genes) #12個

# 檢視結果（這才是資料表）
View(typeI_IFN_df)

#----{所有IFN pathway}
IFN_ALL<-grepl("interferon", ego_BP@result$Description, ignore.case = TRUE)
IFN_ALL_df <- ego_BP@result[IFN_ALL, ]
View(IFN_ALL_df)



