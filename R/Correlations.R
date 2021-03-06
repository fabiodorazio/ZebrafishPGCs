## plot correlation of tpm between PGCs and somatic cells

rna.pgc.256.high <- read.csv('UpregInPGCOnlyAtHigh.txt', sep = '\t')
rna.soma.256.high <- read.csv('UpregInSomaOnlyAtHigh.txt', sep = '\t')

setwd("~/Desktop/PhD-March-2019-backup/FabioRNAseq (2)/2108/RNAseqtables")
tpm <- read.csv('tpmStagesFiltered.txt', sep = '\t')
tpm$colour <- 'black'

setwd("~/Desktop/PhD-March-2019-backup/FabioRNAseq (2)/2108/RNAseqtables/")
germ.plasm <- read.csv('UpregInPGCvsSomaAt256cell.txt', sep = ',')
germ.plasm.id <- unique(germ.plasm$Gene.stable.ID)

for(i in 1:nrow(tpm)){
  if(rownames(tpm)[i] %in% germ.plasm.id){
    tpm$colour[i]="red"
  }
}

tpm.pgc <- tpm[rownames(rna.pgc.256.high),]
tpm.pgc <- na.omit(tpm.pgc)
tpm.soma <- tpm[rownames(rna.soma.256.high),]
tpm.soma <- na.omit(tpm.soma)

plot(log(tpm.pgc$sHighPGC1), log(tpm.pgc$sHighSoma1), pch = 16, col = tpm.pgc$colour,
     xlab = bquote(~Log[2]~ "(tpm) PGC"), ylab = bquote(~Log[2]~ "(tpm) Soma"))


