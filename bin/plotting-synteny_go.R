#!/opt/conda/bin/Rscript --vanilla

library(pheatmap)

erefd<-read.table("Go_summary_topSynteny.tsv", h=T, sep="\t")
newdata <- erefd[order(erefd$Count_significant, decreasing = T),]
rownames(newdata) <- paste(newdata$GO_ID,newdata$GO_term)
df = subset(newdata, select = -c(GO_term,GO_ID,Count_significant) )
df[is.na(df)] <- 1
df2<-as.matrix(df)
my_palette <- colorRampPalette(c("purple","red", "orange", "yellow", "white")) (n=20)
pdf("Go_summary_topSynteny.pdf", width=6, height=5)
pheatmap::pheatmap(log10(head(df2, n=30)), col=my_palette, cluster_rows = F, treeheight_row = 0, treeheight_col = 0, legend=T)
dev.off()

erefd<-read.table("Go_summary_botSynteny.tsv", h=T, sep="\t")
newdata <- erefd[order(erefd$Count_significant, decreasing = T),]
rownames(newdata) <- paste(newdata$GO_ID,newdata$GO_term)
df = subset(newdata, select = -c(GO_term,GO_ID,Count_significant) )
df[is.na(df)] <- 1
df2<-as.matrix(df)
my_palette <- colorRampPalette(c("purple","red", "orange", "yellow", "white")) (n=20)
pdf("Go_summary_botSynteny.pdf", width=6, height=5)
pheatmap::pheatmap(log10(head(df2, n=30)), col=my_palette, cluster_rows = F, treeheight_row = 0, treeheight_col = 0, legend=T)
dev.off()


erefd<-read.table("Go_summary_averhigh.tsv", h=T, sep="\t")
newdata <- erefd[order(erefd$Count_significant, decreasing = T),]
rownames(newdata) <- paste(newdata$GO_ID,newdata$GO_term)
df = subset(newdata, select = -c(GO_term,GO_ID,Count_significant) )
df[is.na(df)] <- 1
df2<-as.matrix(df)
my_palette <- colorRampPalette(c("purple","red", "orange", "yellow", "white")) (n=20)
pdf("Go_summary_averhigh.pdf", width=6, height=5)
pheatmap::pheatmap(log10(head(df2, n=30)), col=my_palette, cluster_rows = F, treeheight_row = 0, treeheight_col = 0, legend=T)
dev.off()

erefd<-read.table("Go_summary_averlow.tsv", h=T, sep="\t")
newdata <- erefd[order(erefd$Count_significant, decreasing = T),]
rownames(newdata) <- paste(newdata$GO_ID,newdata$GO_term)
df = subset(newdata, select = -c(GO_term,GO_ID,Count_significant) )
df[is.na(df)] <- 1
df2<-as.matrix(df)
my_palette <- colorRampPalette(c("purple","red", "orange", "yellow", "white")) (n=20)
pdf("Go_summary_averlow.pdf", width=6, height=5)
pheatmap::pheatmap(log10(head(df2, n=30)), col=my_palette, cluster_rows = F, treeheight_row = 0, treeheight_col = 0, legend=T)
dev.off()

erefd<-read.table("Go_summary_highScore.tsv", h=T, sep="\t")
newdata <- erefd[order(erefd$Count_significant, decreasing = T),]
rownames(newdata) <- paste(newdata$GO_ID,newdata$GO_term)
df = subset(newdata, select = -c(GO_term,GO_ID,Count_significant) )
df[is.na(df)] <- 1
df2<-as.matrix(df)
my_palette <- colorRampPalette(c("purple","red", "orange", "yellow", "white")) (n=20)
pdf("Go_summary_highScore.pdf", width=6, height=5)
pheatmap::pheatmap(log10(head(df2, n=30)), col=my_palette, cluster_rows = F, treeheight_row = 0, treeheight_col = 0, legend=T)
dev.off()

erefd<-read.table("Go_summary_lowScore.tsv", h=T, sep="\t")
newdata <- erefd[order(erefd$Count_significant, decreasing = T),]
rownames(newdata) <- paste(newdata$GO_ID,newdata$GO_term)
df = subset(newdata, select = -c(GO_term,GO_ID,Count_significant) )
df[is.na(df)] <- 1
df2<-as.matrix(df)
my_palette <- colorRampPalette(c("purple","red", "orange", "yellow", "white")) (n=20)
pdf("Go_summary_lowScore.pdf", width=6, height=5)
pheatmap::pheatmap(log10(head(df2, n=30)), col=my_palette, cluster_rows = F, treeheight_row = 0, treeheight_col = 0, legend=T)
dev.off()


