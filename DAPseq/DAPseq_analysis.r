# =====================================================================
#  Peak Annotation in R
#  Goal: Data analysis of peaks associated with TCP15.
#
#  Input:
#    - filte_TC2_tm1_peaks.narrowPeak : bed file of peaks
#    - Ghirsutum_gene_model.gff3 : TM1 genome
#
#  Output:
#    - annotated_TC2_peaks.txt # Peak-associated gene regions
#    - annotated_TC2_gene.txt # Peak-associated genes
#    - annotated_TC2_GO_ontology.txt # Biological processes enriched by genes
#
#  Dependencies:
#    - ChIPseeker GenomicFeatures ggupset dplyr clusterProfiler ggplot2 stringr
# =====================================================================

options(ChIPseeker.downstreamDistance = 3000)

## --- Step 1: Peak annotation
peak1 <- readPeakFile("filte_TC2_tm1_peaks.narrowPeak")
cotton <- makeTxDbFromGFF("Ghirsutum_gene_model.gff3")

peakAnno <- annotatePeak(peak1, tssRegion=c(-3000, 3000), TxDb=cotton)
peakAnno_df <- as.data.frame(peakAnno)

write.table(peakAnno_df, file="annotated_TC2_peaks.txt", sep="\t", quote=FALSE)

draw_data=as.data.frame(cbind(c("Promoter:<=3kb","3'UTR","Other Intron","Distal Intergenic"),c(12.29,0.5,0.13,87.08)))
colnames(draw_data)=c("Feature","ratio")
draw_data$feature=factor(draw_data$Feature,levels=c("Promoter:<=3kb","3'UTR","Other Intron","Distal Intergenic"))
draw_data$ratio=as.numeric(draw_data$ratio)
mycols=c("#cf2c23","#fbd317","#85bad2","#d0d0d0")
ggplot(draw_data, aes(x = "", y = ratio, fill = feature)) +
  geom_bar(width = 1, stat = "identity", color = "black", size=1) +
  coord_polar("y", start = 0)+
  geom_text(aes(label = ratio), color = "black",size=8)+
  scale_fill_manual(values = mycols) +
  theme_void()
  
## --- Step 2: downstream gene analysis
genes=peakAnno_df[peakAnno_df$annotation=="Promoter (2-3kb)" | 
                   peakAnno_df$annotation=="Promoter (1-2kb)" | peakAnno_df$annotation=="Promoter (<=1kb)" | 
                   peakAnno_df$annotation=="5' UTR" | peakAnno_df$annotation=="3' UTR",]
write.table(genes, file="annotated_TC2_gene.txt", sep="\t", quote=FALSE)

## --- Step 3: GO
gene=read.table("annotated_TC2_gene.txt", sep="\t",quote = "")
genes_id2 = gene$geneId
gene_GO=read.csv("GO.csv",stringsAsFactors = FALSE)
go_BP = enricher(gene=genes_id2,
                 TERM2GENE = gene_GO[c('GO', 'GeneID')],
                 TERM2NAME = gene_GO[c('GO', 'Description')],
                 pAdjustMethod= "BH", 
                 pvalueCutoff=0.05, 
                 qvalueCutoff=0.2)
write.table(go_BP, "annotated_TC2_GO.txt", sep = '\t', row.names = FALSE, quote = FALSE)

tmp <- read.delim("annotated_TC2_GO.txt")
gene_GO <- gene_GO[!duplicated(gene_GO$GO), ]
tmp <- merge(tmp, gene_GO[c('GO', 'Ontology')], by.x = 'ID',by.y="GO")
tmp <- tmp[c(10, 1:9)]
tmp <- tmp[order(tmp$pvalue), ]
write.table(tmp, "annotated_TC2_GO_ontology.txt", sep = '\t', row.names = FALSE, quote = FALSE)