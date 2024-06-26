#for pieces of code
## for notes
### titles

library(CAGEr)
library(GenomicRanges)
library(GenomicFeatures) ## for function makeTxdb
library(ggplot2) ## for volcano plot
library(DESeq2)
library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)
library(heatmaps)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)

# ctss of mouse ESC were downloaded from http://fantom.gsc.riken.jp/5/sstar/FF:14345-155G7

### CAGE analysis

## start with a cage set table
### CAGE data


## load mouse tss overlapped with promoters
CAGEset <- readRDS("CAGEsetMouseESCfantom5.rds")
return.cage.anno <- function(x){
  cage.set.ESC <- tagClusters(x,sample = "CAGE_mouseMm9_ESC", returnInterquantileWidth = TRUE,
                            qLow = 0.1,qUp = 0.9)
  ## annotate peaks
  txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
  cage.range <- toGRanges(cage.set.ESC)
  anno.cage <- annotatePeak(cage.range, TxDb = txdb,  annoDb = "org.Mm.eg.db", sameStrand = TRUE, verbose = FALSE)
  cage.esc.anno <- data.frame(anno.cage@anno)
  cage.esc.anno <- subset(cage.esc.anno, cage.esc.anno$ENSEMBL != 'NA')
  return(cage.esc.anno)
}
cage.esc.anno <- return.cage.anno(CAGEset)

  ## discriminate BROAD and SHARP promoters
rownames(cage.esc.anno) <- cage.esc.anno$ensembl_gene_id
sIdx <- subset(cage.esc.anno, cage.esc.anno$interquantile_width <= 9)
bIdx <- subset(cage.esc.anno, cage.esc.anno$interquantile_width > 9)

## import RNAseq data and run DESeq2
rna.read.count <- read.csv('GSE101769_mESC_GeneCounts_1.txt', sep = '\t')
tko <- rna.read.count[,c(3,4,7,8)]
tripleKO <- rna.read.count[,c(1,2,3,4)]
PRC1.KO <- rna.read.count[,c(1,2,5,6)]

condition = c('control', 'control', 'treated', 'treated')

#  start # 
run.deseq.genes <- function(x, condition){
  coldata <- data.frame(condition = condition) #condition = type)
  row.names(coldata) <- colnames(x)
  head(coldata)
  
  dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ condition) #+ condition)
  head(dds)
  dds <- DESeq(dds)
  res <- results(dds) 
  plotMA(res, ylim = c(-10,10))

  return(res)
}
# end #

res.tko <- run.deseq.genes(tripleKO, condition)
res.prc1 <- run.deseq.genes(PRC1.KO, condition)

mcols(res.tko, use.names=TRUE) # check conditions of the experiment

# start #
subset.res <- function(x){
  broad.df <- data.frame(row.names = bIdx$ensembl_gene_id, Tag = 1:dim(bIdx)[1])
  broad.df$Tag <- 'broad'
  sharp.df <- data.frame(row.names = sIdx$ensembl_gene_id, Tag = 1:dim(sIdx)[1])
  sharp.df$Tag <- 'sharp'
  promoter_width <- rbind(broad.df, sharp.df)
  
  sign.genes <- subset(x, x$padj < 0.05)
  sign.genes <- merge(as.data.frame(sign.genes), promoter_width, by = 0)
  rownames(sign.genes) <- sign.genes$Row.names
  sign.genes$Row.names <- NULL

  sub.up <- subset(sign.genes, sign.genes$log2FoldChange > 1)
  sub.down <- subset(sign.genes, sign.genes$log2FoldChange < -1)
  list.x <- list(sub.up, sub.down, sign.genes)
  return(list.x)
}
# end #

sub.res.tko <- subset.res(res.tko)
## FOR TR.KO splits list of dataframes and assign each to a variable
test <- data.frame(sub.res.tko)
for (i in seq(sub.res.tko))
  assign(paste0("dl", i), sub.res.tko[[i]]) ## generates 3 dataframes named df1, df2, df3
## FOR PRC1 splits list of dataframes and assign each to a variable
for (i in seq(sub.res.prc1))
  assign(paste0("dl.prc1", i), sub.res.prc1[[i]]) ## generates 3 dataframes named df1, df2, df3


#### barchart for sharp and broad promoters ####
# start #
count.sharp.and.broad <- function(x, name){
  total.broad <- dim(bIdx)[1]
  total.sharp <- dim(sIdx)[1]
  count.sharp <- dim(subset(x, x$Tag == 'sharp'))[1]
  count.broad <- dim(subset(x, x$Tag == 'broad'))[1]
  normal.count.sharp <- count.sharp/total.sharp
  normal.count.broad <- count.broad/total.broad
  ## plot bar chart
  frame.count <- data.frame(row.names = c('broad', 'sharp'), name = c(normal.count.broad, normal.count.sharp))
  colnames(frame.count) <- name
  return(frame.count)
  
}
 ## count of normalised sharp and broad
frame.upreg <- count.sharp.and.broad(dl1, 'upreg')
frame.downreg <- count.sharp.and.broad(dl2, 'downreg')
frame.non.sign <- count.sharp.and.broad(dl3, 'non sign')
all.combined <- data.frame(frame.upreg, frame.downreg, frame.non.sign)

## plot barchart
plot.barchart <- function (x){
  frame.count <- data.matrix(x)
  barplot(frame.count, col=c('blue', 'orange') , border=NULL, space=0.04, font.axis=2)
}
plot.barchart(all.combined)
### end ###


## go ##
ego <- enrichGO(rownames(dl.prc11), OrgDb = org.Mm.eg.db, keyType = 'ENSEMBL', ont = 'BP',
                pvalueCutoff = 0.1, universe = rownames(rna.read.count))
dotplot(ego, showCategory = 25)

## custom MA plot for DEgenes and promoter width

## select those genes with low tau value meaning that are more likely housekeeping
tau.values <- read.csv('tau_values_MusENCODE_22PC.csv', sep = ',', row.names = 1)
housekeeping <- subset(tau.values, tau.values$Tau < 0.2)
tissue.specific <- subset(tau.values, tau.values$Tau > 0.9)

# start #
custom.volcano <- function(x, highlighted_genes, promoter_width){ # either tissue specific or housekeeping
  ## creates a dataframe of ensembl ids associated to either sharp or broad promoters
  ## remove the row below if split between sharp and broad is not required
  promoters <- data.frame(row.names = promoter_width$ensembl_gene_id, width = 1:dim(promoter_width)[1])
  width.promoter <- deparse(substitute(promoter_width))
  promoters$width <- ifelse(width.promoter == 'sIdx', 'sharp', 'broad')
  
  frame <- data.frame(row.names = rownames(x), Tag = 1:dim(x)[1])
  frame$Tag <- 'gene' # assigns a value to each element of that column
  frame[rownames(highlighted_genes),] <- 'dev' # tags genes
  ## remove the row below if split between sharp and broad is not required
  frame.width <- merge(frame, promoters, by = 0) # merges the frames with promoter widths and cage info
  
  rownames(frame.width) <- frame.width$Row.names
  frame.width$Row.names <- NULL
  # define order of plotting
  #frame$order <- ifelse(frame$Tag=="dev", 1, 2)
  
  x <- data.frame(x)
  x <- merge(x, frame.width, by = 0)
  g <- ggplot(x, aes(x=x$log2FoldChange, y=-log10(x$padj))) + 
    geom_point(aes(color=x$Tag, alpha=x$Tag, size=x$Tag, fill = x$order)) +
    scale_color_manual(values = c(dev = 'orange', gene = 'dark cyan')) + 
    theme_bw() +
    scale_size_manual(values = c(2.5, 2)) + scale_alpha_manual(values = c(dev = 1, gene = 0.02)) +
    geom_vline(xintercept=c(-1, 1), linetype="longdash", colour="black", size=0.4) +
    geom_hline(yintercept=1.3, linetype="longdash", colour="black", size=0.4)
  g
}
custom.volcano(res.prc1, housekeeping)

custom.volcano.plot <- function(results.deseq, promoter_width){ # either sIdx or bIdx
  ## create two data frames were ensembl ids are associated with a tag name indicating sharp or broad promoters
  promoters <- data.frame(row.names = promoter_width$ensembl_gene_id, Tag = 1:dim(promoter_width)[1])
  width.promoter <- deparse(substitute(promoter_width))
  promoters$Tag <- ifelse(width.promoter == 'sIdx', 'sharp', 'broad')

 #sharp.df <- data.frame(row.names = sIdx$ensembl_gene_id, Tag = 1:dim(sIdx)[1])
  #sharp.df$Tag <- 'sharp'
  res <- data.frame(results.deseq)
  res <- na.omit(res)
  merge.promoters <- merge(promoters, res, by = 0)
  #merge.broad <- merge(broad.df, res, by = 0)
  #x <- rbind(merge.sharp, merge.broad)
  x <- merge.promoters
  #x$significance[(x$padj < 0.05) & (x$log2FoldChange > 1)] <- "FDR"
  #x$significance[(x$padj < 0.05) & (x$log2FoldChange < -1)] <- "FC_FDR"
  
  g <- ggplot(x, aes(x=x$log2FoldChange, y=-log10(x$padj))) +
    geom_point(aes(color=x$Tag), alpha= 0.4, size=2) +
    scale_color_manual(values = c(sharp = 'purple', broad = 'dark cyan', dev = 'red'))
    #scale_color_manual(values = ifelse(width.promoter == 'sIdx',  "purple", "dark cyan"))
  #g
  #return(head(x,10))
  
  g1= g + theme_bw(base_size=2) +theme(axis.text.x=element_text(angle=0, size=12, vjust=1),
                                       axis.text.y=element_text(angle=0, size=12, vjust=1),
                                       axis.title=element_text(size=12),legend.position="top",  
                                       legend.key=element_blank(),legend.key.size=unit(0.5, "cm"),
                                       legend.text=element_text(size=8),title=element_text(size=8),
                                       legend.title=element_blank()) + xlab(bquote(~Log[2]~ "fold change")) +
    ylab(bquote(~-Log[10]~italic(Pvalue))) +
    #yintercept is the -log10(pvalue): if pad < 0.05 -> -log10(0.05) = 1.3
    geom_vline(xintercept=c(-2, 2), linetype="longdash", colour="black", size=0.4) +
    geom_hline(yintercept=1.3, linetype="longdash", colour="black", size=0.4)

  g1+ geom_text(aes(label = ifelse((x$padj<0.05) & (abs(x$log2FoldChange)>20), as.character(x$start), "")), 
                hjust = 0, vjust = -0.25)
  
  #return(x) prints new res.all
  g1
}
custom.volcano.plot(res.j1)


###### expression of genes based on GO  #######
###### start  #######

library(biomaRt)

mart=useMart('ENSEMBL_MART_ENSEMBL',
             dataset='drerio_gene_ensembl',
             host="sep2015.archive.ensembl.org")
#find filers
listAttributes(mart) -> list.marts
list.marts$name[grep('go', list.marts$name)]

s4 <- rownames(GermDevGenes)
G_list <- getBM(attributes= c("ensembl_gene_id", "go_id", "hgnc_symbol"), values="", mart= mart)

## retrieve ensembl ids within GO category
## find GO associate with TATA box binding protein (GO = ENSDARG00000014994)
G_list[grep('ENSDARG00000014994', G_list$ensembl_gene_id),]
# 203394 ENSDARG00000014994 GO:0003674            
# 203395 ENSDARG00000014994 GO:0006352

## subset G_list for the GO of interest
Tr.initiation.associated.genes <- subset(G_list, G_list$go_id == 'GO:0006352')
rownames(Tr.initiation.associated.genes) <- Tr.initiation.associated.genes$ensembl_gene_id
Tr.initiation.associated.genes$ensembl_gene_id <- NULL
Tr.initiation.associated.ensembl <- rownames(Tr.initiation.associated.genes)

#### start ####
## make a function to retrieve desired ensembl from GO class
mart=useMart('ENSEMBL_MART_ENSEMBL',
             dataset= 'mmusculus_gene_ensembl',
             host="may2012.archive.ensembl.org")
G_list <- getBM(attributes= c("ensembl_gene_id", "go_id", "hgnc_symbol"), values="", mart= mart)

function.biomart <- function(G_list, gene_ont1, gene_ont2){

  Tr.initiation.associated.genes <- subset(G_list, G_list$go_id == c(gene_ont1, gene_ont2))
  rownames(Tr.initiation.associated.genes) <- Tr.initiation.associated.genes$ensembl_gene_id
  Tr.initiation.associated.genes$ensembl_gene_id <- NULL
  return(Tr.initiation.associated.genes)
}
# GO:0007275 mult organism dev; GO:0009790 embryo development; GO:0009653 tissue morphogenesis; 0007420 brain development
# GO:0009987 cellular processes; GO:0044764 mult organism cellular processes; ->ENSMUSG00000059070<- translation
ensembl.dev <- function.biomart(G_list, 'GO:0007420', 'GO:0009790') 

#### end ####

#### start ####
## plot heatmap
CAGEset <- readRDS("pooled_genes_with_tcs_perGene.rds")
#make txdb to assign promoters


plot.heatmap.promoters <- function(x){
  ## download the CGI coordinates from UCSC table browser
  url <- "http://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/cpgIslandExt.txt.gz"
  con <- gzcon(url(url, open = "r"))
  cgi.gr <- read.table(textConnection(readLines(con)), sep = "\t", col.names = c("bin",
                                                                                 "chrom", "chromStart", "chromEnd", "name", "length", "cpgNum", "gcNum", "perCpg",
                                                                                 "perGc", "obsExp"))
  ## convert to GRanges
  cgi.gr <- GRanges(cgi.gr$chrom, IRanges(start = cgi.gr$chromStart + 1, end = cgi.gr$chromEnd,
                                          names = cgi.gr$name), cpgNum = cgi.gr$cpgNum, gcNum = cgi.gr$gcNum, perCpg = cgi.gr$perCpg,
                    perGc = cgi.gr$perGc, obsExp = cgi.gr$obsExp, seqlengths = seqlengths(BSgenome.Mmusculus.UCSC.mm9))
  
  dominantCTSS_flankSeq <- getSeq(Mmusculus, promotersflank)
  mat_count <- matrix(0, nrow = length(dominantCTSS_flankSeq), ncol = width(dominantCTSS_flankSeq)) 
  # pattern finding:
  patternsOccurrence <- getPatternOccurrenceList(regionsSeq = dominantCTSS_flankSeq,
                                                 patterns = "CG")
  # fill in matrix with 1s at positions where the motif occurs
  mat_count[cbind(patternsOccurrence[[1]]$sequence, patternsOccurrence[[1]]$position)] <- patternsOccurrence[[1]]$value
  
  hm.l = CoverageHeatmap(
    windows=promotersflank,
    track=prom.cage,
    weight=weight,
    coords = c(-500,500),      # change this depending on your window size 
    label = "")
  # within heatmaps package
  
  x <- na.omit(x)
  prom.cage <- GRanges(seqnames = x$chr,
                               ranges=IRanges(start = x$start, end = x$end),
                               strand = x$strand,
                               IQwidth = x$interquantile_width,
                              tpm = x$tpm)
  seqlengths(BSgenome.Mmusculus.UCSC.mm9)
  seqinfo(prom.cage) <- seqinfo(Mmusculus)[seqlevels(prom.cage)]
  names(prom.cage) <- seqnames(prom.cage)
  

  weight <- prom.cage$tpm
  # generate promoters object
  mm9 <- TxDb.Mmusculus.UCSC.mm9.knownGene
  promotersflank <- promoters(mm9, upstream=500, downstream=500)
  promotersflank <- promoters[width(trim(promoters)) == 1000]
  seqinfo(promotersflank) <- seqinfo(Mmusculus)[seqlevels((promotersflank))]
  promflank.ordered <- dropSeqlevels(promotersflank,"chrY_random ",pruning.mode = "coarse")
  names(promflank.ordered) <- seqnames(promflank.ordered)  
  
  hm.l = CoverageHeatmap(
    windows = promflank.ordered,
    track = prom.cage,
    weight=weight,
    coords = c(-1000,500),
    label = "tpm.dominant_ctss")
  
  #smooth heatmap and scale 
  hm_smoothed.l <- smoothHeatmap(hm.l, sigma = c(2, 3), output.size=c(5000, 500))
  scale(hm_smoothed.l) <- quantile(weight, c(0.1, 0.9))
  
  # plot heatmaps
  samplename <- deparse(substitute(x))
  pdf(paste('SmoothHeatmap', samplename, '.pdf'), height = 8, width = 6)
  plotHeatmapList(hm_smoothed.l,
                  cex.label = 1.5,
                  legend = TRUE,
                  legend.pos = "l",
                  legend.width = 0.15,
                  color = c("white", "blue"))
  
  dev.off()
  
}


function.density.heatmap <- function(x, subset){
  rownames(x) <- x$ensembl_gene_id
  subset.names <- rownames(subset)
  x <- x[subset.names,]
  x <- na.omit(x)
  PromotersTSS<-GRanges(seqnames = x$chr,
                                 ranges=IRanges(start = x$dominant_ctss, end = x$dominant_ctss),
                                 strand = x$strand,
                                 interquantileWidth = x$interquantile_width,
                                 seqlengths = seqlengths(Mmusculus))
  ## centers the cage dominant to 0 if up and downstream are equally set
  PromotersTSSflank <- promoters(PromotersTSS, upstream = 600, downstream = 600)
  PromotersTSSflank <- PromotersTSSflank[order(PromotersTSSflank$interquantileWidth),]
  
  PromotersTSSflankSeq <- getSeq(Mmusculus, PromotersTSSflank)
  #pdf(paste('SmoothHeatmap', samplename, '.pdf'), height = 8, width = 6)
  
  plotPatternDensityMap(PromotersTSSflankSeq, c("TATAA"))
                       
  #dev.off()
  #return(PromotersTSSflank)
  
  sIdx <- PromotersTSSflank$interquantileWidth <= 9
  bIdx <- PromotersTSSflank$interquantileWidth > 9
  # plot average dinucleotide profile for sharp promoters
  par(mfrow = c(1,2), mar = c(4.5,4,1,1))

  plotPatternOccurrenceAverage(regionsSeq = PromotersTSSflankSeq[sIdx],
                               patterns = c("WW", "SS"), flankUp = 400, flankDown = 400,
                               smoothingWindow = 3, color = c("red3", "blue3"), cex.axis = 0.9)
  # plot average dinucleotide profile for broad promoters
  plotPatternOccurrenceAverage(regionsSeq = PromotersTSSflankSeq[bIdx],
                               patterns = c("WW", "SS"), flankUp = 400, flankDown = 400,
                               smoothingWindow = 3, color = c("red3", "blue3"), cex.axis = 0.9)

  
}

density.heatmap.pgc.shift <- function.density.heatmap(CAGEset, dl1) #dl1 = upregulated
#### end ####

## import RNAseq dataset
soma.MO <- totalReadsCount[,c(22,24,18,20)]
condition <- c('MO', 'MO', 'wt', 'wt')
soma.MO.res <- as.data.frame(run.deseq.genes(soma.MO, condition))
## subset 
###### end  #######
