library(DiffBind)

idr.T1 <- read.csv("Adults-dripseq-samples-idr-dsA.csv")
idr.T1 <- dba(sampleSheet=idr.T1)
idr.T1 <- dba.count(idr.T1, bUseSummarizeOverlaps=TRUE,summits=FALSE)
dba.plotPCA(idr.T1,  attributes=DBA_CONDITION, label=DBA_ID)

idr.T1 <- dba.normalize(idr.T1, normalize=DBA_NORM_LIB,library=DBA_LIBSIZE_FULL)
idr.T1.contrast <- dba.contrast(idr.T1, categories=DBA_CONDITION, minMembers=2)
idr.T1.contrast
idr.T1.contrast <- dba.analyze(idr.T1.contrast, method=DBA_ALL_METHODS)
dba.show(idr.T1.contrast, bContrasts=TRUE, th=0.1)
idr.T1.out <- dba.report(idr.T1.contrast, method=DBA_DESEQ2, th=0.1)
idr.T1.df <- as.data.frame(idr.T1.out)
write.table(idr.T1.df, file="results/Adults_dsA_DRIPseq_idr_bysex_deseq3.txt", sep="\t", quote=F, row.names=F)
idr.T1.out <- dba.report(idr.T1.contrast, method=DBA_DESEQ2, th=1)
idr.T1.df <- as.data.frame(idr.T1.out)
write.table(idr.T1.df, file="results/Adults_dsA_DRIPseq_idr_all_deseq3.txt", sep="\t", quote=F, row.names=F)