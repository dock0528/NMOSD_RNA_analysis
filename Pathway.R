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
type_l_interferons<-read.csv('RNA_DATA/type l interferon gene list (MsigDB).csv',header=T)
interferon_in_DEGs<-intersect(DEG_genes_df$SYMBOL,type_l_interferons$Gene) #交集後的"IRF7"
#交集>>>ENSG00000185507:3665:IRF7

#-----------{Strictly DEGs}----------------
strict_DEGs_df<-read.csv('RNA_DATA/NMOSD_RNA_DEGS_strictly(Batch Effect Correction).csv',header=T)
strict_DEGs_genes<-strict_DEGs_df$x
strict_DEGs_genes_df<-bitr(strict_DEGs_genes,fromType = 'ENSEMBL',toType=c('ENTREZID','SYMBOL'),OrgDb='org.Hs.eg.db')
#轉成ENTREZID & SYMBOL (從51DEGs -> 52DEGs)
interferon_in_strict_DEGs<-intersect(strict_DEGs_genes_df$SYMBOL,type_l_interferons$Gene) #沒交集

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

