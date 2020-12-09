###  Loading required packages
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(sva)
library(IHW)
library(org.Hs.eg.db)
library(clusterProfiler)
library(fgsea)
library(EnhancedVolcano)
library(ggplot2)
library(grid)
library(gridExtra)


###  Read in the count matrix and prepare Deseq2 data object
countdata = read.table('read_counts_18_02_2020.txt', row.names = 'Gene', header = TRUE)
countdata = countdata[ ,order(names(countdata))]
countdata = countdata[complete.cases(countdata), ]
Samples = c(colnames(countdata))
colnames(countdata)
Therapy = as.factor(c('Yes', 'Yes', 'Yes', 'Yes', 'Yes', 'Yes', 'No', 'No', 'No', 'No', 'No', 'No'))
colData = data.frame(Samples, Therapy)
ddsMat = DESeqDataSetFromMatrix(countData = countdata,
                                colData = colData,
                                design = ~ Therapy)
ddsMat$Therapy = factor(ddsMat$Therapy, levels=c('No', 'Yes'))


###  Filter the count matrix
keep = rowSums(counts(ddsMat) >= 10) >= 3
ddsMat = ddsMat[keep, ]


###  Draw sample distance hatmap, MDS and PCA plots
vsd = vst(ddsMat, blind = FALSE)
sampleDists = dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$Therapy, vsd$Samples, sep = ' - ')
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = Samples)) +
  geom_point(size = 5) + coord_fixed() + ggtitle("MDS with VST data - therapy") +
  theme_grey(base_size = 17) +
  labs(x = 'MDS component 1', y = 'MDS component 2')

plotPCA(vsd, intgroup = c("Therapy"))
theme_grey(base_size = 17)
plotPCA(vsd, intgroup = c("Samples"))
theme_grey(base_size = 17)


###  Test for unaccounted batch effects using sva and add them to the design model
ddsMat1 = DESeq(ddsMat, minReplicatesForReplace = 6)
p <- ncol(attr(ddsMat1,"modelMatrix"))
m = ncol(ddsMat1)
cooksCutoff = qf(.99, p, m - p)
dat  <- counts(ddsMat1, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ Therapy, colData(ddsMat1))
mod0 <- model.matrix(~   1, colData(ddsMat1))
svseq <- svaseq(as.matrix(dat), mod, mod0)
svseq$sv
ddssva = ddsMat1
resultsNames(ddssva)
ddssva$SV1 = svseq$sv[,1] 
ddssva$SV2 = svseq$sv[,2]
ddssva$SV3 = svseq$sv[,3]
design(ddssva) = ~ SV1 + SV2 + SV3 + Therapy


###  Draw sample distance heatmap and PCA plots for the adjusted model
sampleDistMatrix_sva <- as.matrix( sampleDists_sva )
rownames(sampleDistMatrix_sva) <- paste(vsd_sva$Therapy, vsd_sva$Samples, sep = ' - ')
colnames(sampleDistMatrix_sva) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_sva,
         clustering_distance_rows = sampleDists_sva,
         clustering_distance_cols = sampleDists_sva,
         col = colors)

plotPCA(vsd_sva, intgroup = c("Therapy"))
plotPCA(vsd_sva, intgroup = c("Samples"))



###  Perform differential expression analysis with DESeq2
dds = DESeq(ddssva, minReplicatesForReplace = 6)

p <- ncol(attr(dds,"modelMatrix"))
m = ncol(dds)
cooksCutoff = qf(.99, p, m - p)

res = results(dds, alpha = 0.05, cooksCutoff = cooksCutoff,
              pAdjustMethod = 'BH', name = 'Therapy_Yes_vs_No',
              filterFun = ihw)

resOrdered = res[order(res$padj), ]
resOrderedSignif = resOrdered[!is.na(resOrdered$padj) & resOrdered$padj < 0.05, ]
resOrderedSignif2 = resOrderedSignif[abs(resOrderedSignif$log2FoldChange) > 1.5, ]


###  Perform log2FoldChange shrinkage using the 'apeglm' method
lfc_shrink = lfcShrink(dds, res = res, type = 'apeglm', coef = 'Therapy_Yes_vs_No')

res1Ordered = lfc_shrink[order(lfc_shrink$padj), ]
res1OrderedSignif = res1Ordered[!is.na(res1Ordered$padj) & res1Ordered$padj < 0.05, ]
res1OrderedSignif2 = res1OrderedSignif[abs(res1OrderedSignif$log2FoldChange) > 1.5,]


###  Gather statistics for L2FC of upregulated and downregulated DEG`s
up = res1OrderedSignif2[res1OrderedSignif2$log2FoldChange > 0, ]
down = res1OrderedSignif2[res1OrderedSignif2$log2FoldChange < 0, ]
median_up = median(up$log2FoldChange)
median_down = median(down$log2FoldChange)
iqr_up = IQR(up$log2FoldChange)
iqr_down = IQR(down$log2FoldChange)

###  Convert ENSEMBL keys to SYMBOL keys
ens.str = substr(res1OrderedSignif2[, 'Genes'], 1, 18266)
res1OrderedSignif2$Symbol = mapIds(org.Hs.eg.db,
                                   keys = ens.str,
                                   column = "SYMBOL",
                                   keytype = "ENSEMBL",
                                   multiVals = "first")


###  Draw p-value histogram
hist(lfc_shrink$pvalue, col = "lavender", main = "Terapija: Ir vs Nav", xlab = "p-values")


###  Perform GSEA analysis using the fgsea package and MSigDB Hallmark gene sets
lfc_shrink_table = as.data.frame(lfc_shrink)
lfc_shrink_table$ENSEMBL = rownames(lfc_shrink_table)
lfc_shrink_table = lfc_shrink_table[complete.cases(lfc_shrink_table), ]

translation_entrez = as.data.frame(bitr(c(lfc_shrink_table$ENSEMBL), 
                                        fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db'))
translation_symbol = as.data.frame(bitr(c(lfc_shrink_table$ENSEMBL),
                                        fromType = 'ENSEMBL', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db'))

lfc_shrink1 = merge(lfc_shrink, translation_entrez, by = 'ENSEMBL')
lfc_shrink2 = merge(lfc_shrink1, translation_symbol, by = 'ENSEMBL')

gseaDat = lfc_shrink2[which(duplicated(lfc_shrink2$ENTREZID) == F), ]

ranks = gseaDat$log2FoldChange
names(ranks) = gseaDat$ENTREZID
ranks = sort(ranks)
gmt_path = gmtPathways('msigdb.v7.1.entrez.gmt')

msig = fgsea(pathways = gmt_path, stats = ranks, minSize = 10, maxSize = 1000, eps = 0)
msig_df_sig = msig_df[msig_df$padj < 0.05, ]


###  Perform GSEA analysis using the KEGG pathway set
ranks = sort(ranks, decreasing = TRUE)
gseaRes_KEGG = gseKEGG(geneList = ranks,
                       organism = 'hsa',
                       minGSSize = 10,
                       maxGSSize = 1000,
                       pvalueCutoff = 0.05)
gseaKEGG_results <- gseaRes_KEGG@result
gseaKEGG_results = as.data.frame(gseaKEGG_results)





###  Create volcano plot of the relationship between logFoldChange and p-value
EnhancedVolcano(lfc_shrink,
                lab = lfc_shrink$Symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1.5,
                pCutoff = 7.68e-4,
                xlim = c(-6.5, 6.5),
                ylim = c(0, 11),
                pointSize = 3.5,
                labSize = 3.5,
                caption = "",
                title = "",
                subtitle = "",
                gridlines.major = FALSE,
                gridlines.minor = FALSE)


###  Create a heatmap to cluster samples and DEG`s by their expression similiarity.
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 1, hjust = 0.5, rot = 0, gp = gpar(...))
  return(res)
}

assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

heat_colors <- brewer.pal(9, "YlOrRd")
labels_row = c(res1OrderedSignif2$Symbol)
degs = c(res1OrderedSignif2$Genes)
mat  <- assay(vsd_sva)[degs,]
mat = mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd_sva)[, c("Samples","Therapy")])
pheatmap(mat, annotation_col = anno, fontsize = 10,
         fontsize_row = 13, border_color = FALSE, 
         fontsize_col = 16, labels_row = labels_row, color = heat_colors, annotation_legend = FALSE)






###  Draw boxplots for raw read count distribution for significant DEG`s
plot_sva_boxplots  = function (dati) {
  m <- estimateSizeFactors(dati)
  sva_counts = as.data.frame(counts(m, normalized=TRUE))
  head(sva_counts)
  
  for (gene in res1OrderedSignif2$Genes) {
    df = as.data.frame(t(sva_counts[gene, ]))
    
    padj = formatC(res1OrderedSignif2[gene, 5], format = 'e', digits = 2)
    log2fold = round(res1OrderedSignif2[gene, 2], 2)
    symbol = res1OrderedSignif2[gene, 7]
    
    names(df)[1] = 'Count'
    df$Therapy = Therapy
    df$Label = rownames(df)
    df[, 1] = df[, 1] + 0.5
    
    print(res1OrderedSignif2[gene, 7])
    
    par(mfrow=c(1,1))
    boxplots = ggplot(df, aes(x = Therapy, y = Count, color = Therapy)) +
      geom_boxplot(alpha = 0.8, outlier.colour = NA, coef = 500) + 
      scale_y_log10() +
      geom_point(aes(fill = Therapy), size = 5, shape  = 21, position = position_jitterdodge()) +
      theme(text = element_text(size = 18),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank()) +
      annotate('text', x = -Inf, y = Inf, hjust = -0.01, vjust = 1,
               label = paste0('P-adjusted: ', padj)) +
      annotate('text', x = -Inf, y = Inf, hjust = -0.01, vjust = 2.3,
               label = paste0('L2FC: ', log2fold))
    
    ggsave(boxplots, filename = paste0('graphs/raw_boxplots/', symbol, '_by_therapy_raw.png'), 
           height = 7, width = 7)
    
  }
}

plot_sva_boxplots(ddsMat1)


###  Draw SVA weighted/unweighted model comparison barplots for significant DEG`s
plot_mod_diff = function (no, yes) {
  m1 = estimateSizeFactors(no)
  m2 = estimateSizeFactors(yes)
  m1_counts = as.data.frame(counts(m1, normalized = TRUE))
  m2_counts = as.data.frame(counts(m2, normalized = TRUE))
  
  #print(head(m1_counts, 3))
  #print(head(m2_counts, 3))
  
  for (gene in res1OrderedSignif2$Genes) {
    
    m1_t = as.data.frame(t(m1_counts[gene,]))
    m2_t = as.data.frame(t(m2_counts[gene,]))
    
    m1_t$Sample = rownames(m1_t)
    m2_t$Sample = rownames(m2_t)
    
    names(m1_t)[1] = 'Count'
    names(m2_t)[1] = 'Count'
    
    m1_t$model = 'mod0'
    m2_t$model = 'mod1'
    
    m1_m2 = as.data.frame(rbind(m1_t, m2_t))
    
    plots = ggplot(m1_m2, aes(x=Sample, y=Count, fill=model)) +
      geom_bar(stat='identity', position='dodge')
    
    ggsave(plots, filename = paste0('graphs/publication_graphs/model_barplots/', gene, '_test.png'))
    
  }
}

plot_mod_diff(ddsMat, ddssva)


###  Create black and white boxplots for the seven most differentially 
expressed genes based on their L2FC values

list_of_genes = c('ENSG00000183742', 'ENSG00000083307', 'ENSG00000181143', 'ENSG00000181234', 'ENSG00000154997', 'ENSG00000186732', 'ENSG00000105088')

res1OrderedSignif3 = res1OrderedSignif2[which(rownames(res1OrderedSignif2) %in% list_of_genes), ]

plot_sva_boxplots_top = function (data) {
  m <- estimateSizeFactors(data)
  sva_counts = as.data.frame(counts(m, normalized=TRUE))
  sva_counts = sva_counts[which(rownames(sva_counts) %in% list_of_genes), ]
  
  p = list()
  
  for (gene in res1OrderedSignif3$Genes) {
    df = as.data.frame(t(sva_counts[gene, ]))
    
    padj = formatC(res1OrderedSignif3[gene, 5], format = 'e', digits = 2)
    log2fold = round(res1OrderedSignif3[gene, 2], 2)
    symbol = res1OrderedSignif3[gene, 7]
    
    names(df)[1] = 'Count'
    df$Therapy = Therapy
    df$Label = rownames(df)
    df[, 1] = df[, 1] + 0.5
    
    print(res1OrderedSignif2[gene, 7])
    
    p[[gene]] = ggplot(df, aes(x = Therapy, y = Count, group = Therapy)) +
      geom_boxplot(alpha = 0.8, outlier.colour = NA, coef = 500) + 
      scale_y_log10() +
      geom_point(aes(color = Therapy), size = 5, shape = 16, position = position_jitterdodge(),
                 alpha = 0.5) +
      scale_colour_manual(values=rep("black",length(Therapy))) +
      #geom_point(aes(fill = Therapy, color = '#000000'), size = 5, shape  = 21, 
      #           position = position_jitterdodge(), color = '#000000') +
      labs(title = symbol) +
      theme_bw() +
      theme(text = element_text(size = 18),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            legend.position = 'none',
            plot.title = element_text(size = 15)) +
      annotate('text', x = -Inf, y = Inf, hjust = -0.01, vjust = 1,
               label = paste0('P-adjusted: ', padj)) +
      annotate('text', x = -Inf, y = Inf, hjust = -0.01, vjust = 2.3,
               label = paste0('L2FC: ', log2fold))
  }
  tiff('graphs/publication_graphs/bw_boxplots/combined_bw_boxplots_lzw.tif',
       units = 'in', w = 15, h = 15, compression = 'lzw', res = 1200)
  do.call(grid.arrange, p)
  dev.off()
}

plot_sva_boxplots_top(ddsMat1)



###  Calculate Kendall`s correlation between Ki-67/Knosp index and read counts for the significant DEG`s
knosp_index = c(2, 1, 4, 1, 0, 4, 1, 1, 4, 2, 3, 1)
ki_index = c(1, 2, 1, 4, 4, 1, 3, 1, 2, 3, 1, 2)

index_correlation = function(x, y, z) {
  
  for (i in rownames(x)) {
    
    gene = z[i, ][, 7]
    
    counts = as.numeric(x[i, ])
    factor = c(y)
    
    data_frame = data.frame(read_counts = counts, knosp_grade = factor, Therapy = Therapy)
    frame_up = as.data.frame(data_frame[1:6, ])
    frame_down = as.data.frame(data_frame[7:12, ])
    
    corr_yes = cor.test(frame_up$knosp_grade, frame_up$read_counts, method = 'kendall')
    corr_no = cor.test(frame_down$knosp_grade, frame_down$read_counts, method = 'kendall')
    
    coef_yes = round(as.numeric(corr_yes$estimate), 2)
    pval_yes = formatC(as.numeric(corr_yes$p.value), format = 'e', digits = 2)
    coef_no = round(as.numeric(corr_no$estimate), 2)
    pval_no = formatC(as.numeric(corr_no$p.value), format = 'e', digits = 2)
    
    Therapy = as.factor(Therapy)
    
    plots = ggplot(data_frame, aes(x = knosp_grade, y = read_counts, color = Therapy, 
                                   shape = Therapy)) +
      scale_color_manual(values = c('red', 'blue')) +
      geom_point(size = 4) +
      geom_smooth(aes(group = Therapy, color = Therapy, fill = Therapy), method = lm, linetype = 'dashed', fullrange = T) +
      guides(
        shape = guide_legend(
          override.aes = list(color = rep(c('red', 'blue'), each = 1))), color = F) +
      annotate(geom = 'text', x = -Inf, y = Inf, hjust = -0.01, vjust = 1,
               label = paste0('Gene: ', gene)) +
      annotate(geom = 'text', x = -Inf, y = Inf, hjust = -0.01, vjust = 2.2,
               label = paste0('Tau Yes: ', coef_yes)) +
      annotate(geom = 'text', x = -Inf, y = Inf, hjust = -0.01, vjust = 3.6,
               label = paste0('Pval - Yes: ', pval_yes)) +
      annotate(geom = 'text', x = -Inf, y = Inf, hjust = -0.01, vjust = 55,
               label = paste0('Tau - No: ', coef_no)) +
      annotate(geom = 'text', x = -Inf, y = Inf, hjust = -0.01, vjust = 56.4,
               label = paste0('Pval - No: ', pval_no)) +
      xlab('Ki - 67 index group') + ylab('Read count')
    
    ggsave(plots, filename = paste0('correlation_graphs_6_k_index/', gene, 
                                    '_corr_k_index.png'))
  }
}

index_correlation(sva_counts_deg, ki_index, res1OrderedSignif2)
index_correlation(sva_counts_deg, knosp_index, res1OrderedSignif2)


