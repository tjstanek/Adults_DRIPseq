directory <- "htseq_ct/"

sampleFiles <- grep(".txt",list.files(directory),value=TRUE)
sampleCondition <- c("F379.1", "F379.2","F732.1","F732.2","M379.1", "M379.2","M732.1","M732.2")
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)

library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq

keep <- rowSums(counts(ddsHTSeq)) > 0
ddsHTSeq0 <- ddsHTSeq[keep,]

rld <- rlog(ddsHTSeq0, blind=TRUE)
head(assay(rld), 3)
write.table(assay(rld), "RNAseq_htseq_rlog_normalized.txt", sep = "\t", row.names = TRUE, quote = FALSE, col.names = FALSE)