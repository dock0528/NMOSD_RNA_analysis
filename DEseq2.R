setwd("C:/Users/JANE/Desktop/Wang實驗室/NMOSD研究計畫/RNA/NMOSD_RNA_analysis")
#BiocManager::install('DESeq2')
library(DESeq2)
library(ggplot2)


Raw_count_merged<-read.csv('./RNA_DATA/Raw_count_merged_matrix.csv',row.names = 1,header=T)

#Raw count 
count_df<-as.matrix(Raw_count_merged)
rownames(count_df)<-gsub("\\.\\d+$", "",rownames(count_df)) #去除小數點以後的值

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

#建立deseq dataset
dds<-DESeqDataSetFromMatrix(countData = count_df, #順序要排對
                            colData = sample_df,
                            design=~Group)

# 保留至少在 1 個樣本裡有表達的基因
keep <- rowSums(counts(dds) > 0) >= 1 # counts(dds) 回傳一個 [基因 × 樣本] 的矩陣 
#每一行（每個基因）裡 TRUE 的個數加起來
dds  <- dds[keep, ]
#做計算
dds<-DESeq(dds)

res<-results(dds)
head(results(dds,tidy=T)) #tidy=T整齊回傳
#表達量太低>>>padj=NA


#summary of differential gene expression
summary(res)

#sort summary list by p-value
res<-res[order(res$padj),]
head(res)

#plotCounts (use plotCounts fxn to compare tumor v.s normal)
par(mfrow=c(2,3)) #創建2列3欄的圖形網格

plotCounts(dds,gene='ENSG00000105607',intgroup = 'Group')
plotCounts(dds,gene='ENSG00000115841',intgroup = 'Group')


# 先存一下原本參數，方便最後還原
op <- par(no.readonly = TRUE)
#volcano plot
#reset par
par(mar = c(5, 4, 4, 6) + 0.1)

#Make basic volcano plot
with(res,plot(log2FoldChange,-log10(pvalue),pch=20,main='NMOSD RNA DEGs',xlim=c(-3,3),ylim=c(0,40),col='grey'))

#Add color point
#blue: (log2FC)<-1 & padj<0.05 downregulation
#red: (log2FC)>1 & padj<0.05 upregulation
# down‐regulated
with(subset(res, padj < 0.05 & log2FoldChange < -1),
     points(log2FoldChange, -log10(pvalue), pch=20, col='#8FA4FF')
)

# up‐regulated
with(subset(res, padj < 0.05 & log2FoldChange >  1),
     points(log2FoldChange, -log10(pvalue), pch=20, col='#FF8080')
)


# 3. 允許畫到圖外
par(xpd = TRUE)

# 4. 在圖外加圖例，title="Change"，inset=c(-0.2,0) 往右外推 20%
legend("topright",
       inset   = c(-0.2, 0),
       legend  = c("Up", "Down"),
       title   = "Change",
       pch     = 20, #點大小
       col     = c("#FF8080", "#8FA4FF"),
       pt.cex  = 1.4, #點符號放大倍數
       bty     = "n")

# 4. 還原 par 設定
par(op)

#顯著基因
Up_DEGs<-row.names(subset(res,padj<0.05 & log2FoldChange>1)) #length(Up_DEGs):4500
Down_DEGs<-row.names(subset(res,padj<0.05 & log2FoldChange < -1)) #length(Down_DEGs):5850
DEGS<-row.names(subset(res,padj<0.05 & abs(log2FoldChange)>1)) #DEGs numbers:10350
write.csv(DEGS,file='RNA_DATA/NMOSD_RNA_DEGS.csv',row.names = F)




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

