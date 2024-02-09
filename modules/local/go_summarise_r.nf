process TEMP_R_PLOTS {

   label 'process_single'
   tag "$sample_id"
   container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'
   publishDir "$params.outdir/GO_results" , mode: "copy"

   input:
   path("*")

   output:
   path( "*.pdf" ), emit: go_summary_pdf
   '''
   #!/usr/bin/Rscript
   library(pheatmap)
   
   erefd<-read.table("Go_summary_topSynteny.tsv", h=T, sep="\t")
   newdata <- erefd[order(erefd$Count_significant, decreasing = T),]
   rownames(newdata) <- paste(newdata$GO_ID,newdata$GO_term)
   df = subset(newdata, select = -c(GO_term,GO_ID,Count_significant) )
   df[is.na(df)] <- 1
   df2<-as.matrix(df)
   df3<-as.data.frame(df2)
   my_mean<-mean(df2)
   if (my_mean == 1){
   	#Do nothing, no results were significant, so plotting will fail
   } else {
   	my_palette <- colorRampPalette(c("purple","red", "orange", "yellow", "white")) (n=20)
   	pdf("Go_summary_topSynteny.pdf", width=6, height=5)
   	pheatmap::pheatmap(log10(head(df2, n=30)), col=my_palette, cluster_rows = F, treeheight_row = 0, treeheight_col = 0, legend=T)
   	dev.off()
   }
   
   
   erefd<-read.table("Go_summary_botSynteny.tsv", h=T, sep="\t")
   newdata <- erefd[order(erefd$Count_significant, decreasing = T),]
   rownames(newdata) <- paste(newdata$GO_ID,newdata$GO_term)
   df = subset(newdata, select = -c(GO_term,GO_ID,Count_significant) )
   df[is.na(df)] <- 1
   df2<-as.matrix(df)
   my_mean<-mean(df2)
   if (my_mean == 1){
   	#Do nothing, no results were significant, so plotting will fail
   } else {
   	my_palette <- colorRampPalette(c("purple","red", "orange", "yellow", "white")) (n=20)
   	pdf("Go_summary_botSynteny.pdf", width=6, height=5)
   	pheatmap::pheatmap(log10(head(df2, n=30)), col=my_palette, cluster_rows = F, treeheight_row = 0, treeheight_col = 0, legend=T)
   	dev.off()
   }
   
   
   erefd<-read.table("Go_summary_averhigh.tsv", h=T, sep="\t")
   newdata <- erefd[order(erefd$Count_significant, decreasing = T),]
   rownames(newdata) <- paste(newdata$GO_ID,newdata$GO_term)
   df = subset(newdata, select = -c(GO_term,GO_ID,Count_significant) )
   df[is.na(df)] <- 1
   df2<-as.matrix(df)
   my_mean<-mean(df2)
   if (my_mean == 1){
   	#Do nothing, no results were significant, so plotting will fail
   } else {
   	my_palette <- colorRampPalette(c("purple","red", "orange", "yellow", "white")) (n=20)
   	pdf("Go_summary_averhigh.pdf", width=6, height=5)
   	pheatmap::pheatmap(log10(head(df2, n=30)), col=my_palette, cluster_rows = F, treeheight_row = 0, treeheight_col = 0, legend=T)
   	dev.off()
   }
   
   erefd<-read.table("Go_summary_averlow.tsv", h=T, sep="\t")
   newdata <- erefd[order(erefd$Count_significant, decreasing = T),]
   rownames(newdata) <- paste(newdata$GO_ID,newdata$GO_term)
   df = subset(newdata, select = -c(GO_term,GO_ID,Count_significant) )
   df[is.na(df)] <- 1
   df2<-as.matrix(df)
   my_mean<-mean(df2)
   if (my_mean == 1){
   	#Do nothing, no results were significant, so plotting will fail
   } else {
   	my_palette <- colorRampPalette(c("purple","red", "orange", "yellow", "white")) (n=20)
   	pdf("Go_summary_averlow.pdf", width=6, height=5)
   	pheatmap::pheatmap(log10(head(df2, n=30)), col=my_palette, cluster_rows = F, treeheight_row = 0, treeheight_col = 0, legend=T)
   	dev.off()
   }
   
   erefd<-read.table("Go_summary_highScore.tsv", h=T, sep="\t")
   newdata <- erefd[order(erefd$Count_significant, decreasing = T),]
   rownames(newdata) <- paste(newdata$GO_ID,newdata$GO_term)
   df = subset(newdata, select = -c(GO_term,GO_ID,Count_significant) )
   df[is.na(df)] <- 1
   df2<-as.matrix(df)
   my_mean<-mean(df2)
   if (my_mean == 1){
   	#Do nothing, no results were significant, so plotting will fail
   } else {
   	my_palette <- colorRampPalette(c("purple","red", "orange", "yellow", "white")) (n=20)
   	pdf("Go_summary_highScore.pdf", width=6, height=5)
   	pheatmap::pheatmap(log10(head(df2, n=30)), col=my_palette, cluster_rows = F, treeheight_row = 0, treeheight_col = 0, legend=T)
   	dev.off()
   }
   
   erefd<-read.table("Go_summary_lowScore.tsv", h=T, sep="\t")
   newdata <- erefd[order(erefd$Count_significant, decreasing = T),]
   rownames(newdata) <- paste(newdata$GO_ID,newdata$GO_term)
   df = subset(newdata, select = -c(GO_term,GO_ID,Count_significant) )
   df[is.na(df)] <- 1
   df2<-as.matrix(df)
   my_mean<-mean(df2)
   if (my_mean == 1){
   	#Do nothing, no results were significant, so plotting will fail
   } else {
   	my_palette <- colorRampPalette(c("purple","red", "orange", "yellow", "white")) (n=20)
   	pdf("Go_summary_lowScore.pdf", width=6, height=5)
   	pheatmap::pheatmap(log10(head(df2, n=30)), col=my_palette, cluster_rows = F, treeheight_row = 0, treeheight_col = 0, legend=T)
   	dev.off()
   }
   '''
}
