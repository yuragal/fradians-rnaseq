dir="~/matrosov/axSA_assemblies/drap/stranded/trinity_norm_altogether/abundance_estimation"

#ExN50 plot
#See https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats#contig-ex90n50-statistic-and-ex90-transcript-count
exn50=read.table(paste0(dir,'/salmon.isoform.TMM.EXPR.matrix.ExN50.stats'),sep='\t',header=TRUE)
gl=seq(0,25,5)
svg('ExN50.svg',width=6,height=5)
par(mar=c(5,4,4,5)+.1)
plot(exn50$E,exn50$ExN50,type='l',col='red',xlab="Expression percentile",ylab="N50, bp")
par(new=TRUE)
plot(exn50$E,exn50$num_transcripts,type='l',col='blue',xaxt="n",yaxt="n",xlab="",ylab="")
axis(4,labels=paste0(gl,'K'),at=gl*1000)
mtext("Number of transcripts",side=4,line=3)
legend("topleft",col=c("red","blue"),lty=1,legend=c("ExN50","Transcripts"))
abline(v=96, lty=3)
dev.off()
#axis(1,labels=c(96),at=c(96))

#Nx plots for all assemblies
#cd /home/ygalachyants/axSA_assemblies/drap/stranded/trinity_norm_altogether/trinotate
#cat transcripts.fa | awk 'NR%2==1{printf("%s\t",$1);next}NR%2==0{printf("%d\n",length($1))}' | cut -f2 | sort -k1nr > transcripts.len
#ln -s transcripts.len trinity_norm_altogether.len
#cat ../../metaRF_trinity/transcripts_fpkm_1.fa | awk 'NR%2==1{printf("%s\t",$1);next}NR%2==0{printf("%d\n",length($1))}' | cut -f2 | sort -k1nr > metaRF_trinity.len
#cat ../../metaRF_oases/transcripts_fpkm_1.fa | awk 'NR%2==1{printf("%s\t",$1);next}NR%2==0{printf("%d\n",length($1))}' | cut -f2 | sort -k1nr > metaRF_oases.len
#cat ../../oases_norm_altogether/transcripts_fpkm_1.fa | awk 'NR%2==1{printf("%s\t",$1);next}NR%2==0{printf("%d\n",length($1))}' | cut -f2 | sort -k1nr > oases_norm_altogether.len
tl=list(
	meta_oases=read.table(paste0(dir,'/../trinotate/metaRF_oases.len')),
	merge_oases=read.table(paste0(dir,'/../trinotate/oases_norm_altogether.len')),
	meta_tr=read.table(paste0(dir,'/../trinotate/metaRF_trinity.len')),
	merge_tr=read.table(paste0(dir,'/../trinotate/trinity_norm_altogether.len'))
)

for(k in c(1:4)){
	tl[[k]]$cum_len=sapply(rownames(tl[[k]]), function(i) sum(tl[[k]][1:i,1]))
	tl[[k]]$Nx=tl[[k]]$cum_len/tl[[k]]$cum_len[dim(tl[[k]])[1]]
}

lapply(tl, function(i) max(i[2])) %>% unlist /10^6

colv=c('blue','green','orange','red')
line_names=c(paste(c('AM','MA'),c(rep('oases',2),rep('trinity',2)),sep='-'))

svg('Nx.svg',width=6,height=5)
plot(tl[[1]]$Nx,tl[[1]]$V1,type='l',col=colv[1], xaxs='i',yaxs='i',log='y',xlab='Normalized assembly length',ylab='Contig lenght, bp')
lines(tl[[2]]$Nx,tl[[2]]$V1,type='l',col=colv[2])
lines(tl[[3]]$Nx,tl[[3]]$V1,type='l',col=colv[3])
lines(tl[[4]]$Nx,tl[[4]]$V1,type='l',col=colv[4])
legend("topright",col=colv,lty=1,legend=line_names)
abline(v=0.5, lty=3)
dev.off()

#Abundance estimation from Salmon data
# See for more details:
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html


library(DESeq2)
library(magrittr)
library(tximport)
#library(tximportData)

samples <- read.table(file.path(dir,"samples.1"), header=FALSE)
colnames(samples)[2]='file'

samples$rep=gsub('.*_rep','',samples$file)
#Culture type -- "normal exponential" -- NExp -- or "dark-synchronized, light-treated" -- DSLT
samples$ctype=gsub('_rep[0-9]+','',samples$file) %>% gsub('-[0-9]+','',.) %>% gsub('lag','DSLT',.) %>% gsub('exp','NExp',.)
#TOT -- time of treatment
samples$TOT=c('-','-',0,0,20,20,40) %>% factor
samples$id=paste0(samples$ctype,samples$rep,'-',samples$TOT) %>% gsub('--','',.)
s=samples[,c(8,6,5,7)]
files <- file.path(dir,paste0("salmon_",samples$file), "quant.sf")
names(files) <- samples$id
#tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))

#Set txOut=TRUE so thet to output transcpirt abundances only
txi=tximport(files,type = 'salmon',txOut = TRUE)
#View(txi$counts)

ddsTxi=DESeqDataSetFromTximport(txi, colData = s, design = ~ ctype)
dds=DESeq(ddsTxi)


library(dplyr)
library(ggplot2)
library(vsn)

ntd=normTransform(dds)
vsd=vst(dds, blind = FALSE)
rld=list(NExp_DSLT=rlog(dds, blind = FALSE))

df <- bind_rows(
	as_data_frame(assay(ntd)[, 1:2]) %>% mutate(transformation = "log2(x + 1)"),
	as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
	as_data_frame(assay(rld[[1]])[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 50) + coord_fixed() + facet_grid( . ~ transformation)

source('~/Desktop/LIN/axSA/RNAseq/26042018/functions.R')
p1=meanSdPlot(assay(ntd),plot=FALSE ,show.legend=FALSE)$gg + ggtitle('ntd')
p2=meanSdPlot(assay(vsd),plot=FALSE ,show.legend=TRUE)$gg + ggtitle('vsd')
p3=meanSdPlot(assay(rld[[1]]),plot=FALSE ,show.legend=TRUE)$gg + ggtitle('rld')
p=multiplot(p1,p2,p3,cols=3)

#Exploratory analysis plots
library(pheatmap)
library(RColorBrewer)
library(viridis)

sampleDists <- dist(t(assay(rld[[1]])))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- rld[[1]]$id
colnames(sampleDistMatrix) <- NULL
#colors <- colorRampPalette( rev(brewer.pal(7, "Blues")) )(255)
p1=pheatmap(sampleDistMatrix,
	clustering_distance_rows = sampleDists,
	clustering_distance_cols = sampleDists,
	col = viridis(255),
	cellwidth = 30,
	cellheight = 30
	)

pcaData=plotPCA(rld$NExp_DSLT, intgroup = c('ctype','TOT'),returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p2=ggplot(pcaData, aes(PC1, PC2, shape=ctype,color=TOT)) +
	geom_point(color = "black", size = 10) +
	geom_point(size=9) +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
	coord_fixed(ratio=3.5) +
	theme_linedraw() + 
	scale_shape_discrete(name="Culture\ntype") +
	scale_color_viridis(discrete = TRUE,name="Time\nof\ntreatment",alpha = 0.9) +
	xlim(-100,50) + ylim(-15,25)

pcaData=plotPCA(rld$Dark_Light, intgroup = c('TOT','Light'),returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p3=ggplot(pcaData, aes(PC1, PC2, shape=Light,color=TOT)) +
	geom_point(color = "black", size = 10) +
	geom_point(size=9) +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
	coord_fixed(ratio=1.2) +
	theme_linedraw() + 
	scale_shape_discrete(name="Light\nexposure") +
	scale_color_viridis(discrete = TRUE,name="Time\nof\ntreatment",alpha = 0.9) +
	xlim(-20,15) + ylim(-15,15)
ggsave(p1,filename = 'Exploratory_sample2sample_rld_clustering.svg',width=5,height = 4)
ggsave(p2,filename = 'Exploratory_samples_rld_PCA.svg',width=5,height = 4)
ggsave(p3,filename = 'Exploratory_samples_rld_PCA_DSLT.svg.svg',width=5,height = 4)
# mds <- as.data.frame(colData(rld))  %>%
# 	cbind(cmdscale(sampleDistMatrix))
# ggplot(mds, aes(x = `1`, y = `2`, color = V1)) +
# 	geom_point(size = 3) + coord_fixed()

rm(p1)
rm(p2)
rm(p3)
rm(ntd)
rm(vsd)
rm(rld)


#Set lfcThreshold to 1 and alpha to 1e-2
#https://support.bioconductor.org/p/101504/#101507
#This seems quite conservative. However, if one need more, set lfc to 2 and a


#!!! UPDATE
#Filtering should be performed in such a way that uses LFC threshold
#to compute significance upon testing procedure and adjusted p-values.
#
#See https://support.bioconductor.org/p/101504/ for more info.
#So, 
# 1) get DESeq object (optionally, with some initial filtering)
# 2) perform testing of genes according to experimental design (functions 'results' or 'lfcShrink')
# 3) then choose genes according to adjusted p-value, q-value, s-value etc., depending on testing procedure used
# 
library(apeglm)
#dds=dds[ rowSums(counts(dds)) > 1, ]
dds=ddsTxi[ rowSums(counts(ddsTxi)) > 10 & !names(ddsTxi) %in% cont_ids$tid, ]
dds=DESeq(dds)
#!!! Important to have LFC-values consistent with sequecne of contrast levels
#NExp vs DSLT
dds@colData$ctype %<>% as.character() %>% factor(.,levels=c('NExp','DSLT'))
# LFC_treshold=2
# resShrink %>% summary
# out of 25525 with nonzero total read count
# s-value < 0.005
# LFC > 2.00 (up)    : 1129, 4.4%
# LFC < -2.00 (down) : 1549, 6.1%
# resShrink %>% summary
# out of 25525 with nonzero total read count
# s-value < 0.005
# LFC > 1.00 (up)    : 2823, 11%
# LFC < -1.00 (down) : 3133, 12%
LFC_treshold=1
#res=results(dds,lfcThreshold = LFC_treshold, alpha=0.05)
res <- lfcShrink(dds, coef=2, type="apeglm", lfcThreshold = LFC_treshold)
#resShrinkAshr <- lfcShrink(dds, coef="ctype_NExp_vs_DSLT", type="ashr",lfcThreshold = LFC_treshold)
# n1=subset(resShrink,svalue<=0.005) %>% rownames()
# n2=subset(res,padj<=0.005) %>% rownames()
# setdiff(n1,n2) %>% length()
# setdiff(n2,n1) %>% length()
# 
# Nexp_vs_DSLT=subset(resShrink,svalue<=0.005)
deseq_results=list()
contrast='NExp_DSLT'
deseq_results[[contrast]]=list(ddstxi=ddsTxi,dds=dds,res=res)


#It seems apeglm reveals twice more genes than normal test while significance is comparable
#DESeq2::plotMA(res, xlim=c(0.5,100000), cex=.3)
DESeq2::plotMA(resShrink, ylim=c(-1,1)*8,xlim=c(0.5,100000), cex=.5)
abline(v=c(1,10,100,1000,10000), col="dodgerblue", lwd=1,lty=3)
abline(h=c(-1,1), col="dodgerblue",  lwd=1,lty=3)


#Compare subgroup of DSLT-samples
contrast='Dark_Light'
ddsTxiLight=ddsTxi
#ddsTxiLight@colData$Light=c('always','always','non_treated','non_treated','treated','treated','treated') %>% factor
ddsTxiLight@colData$Light=c(rep('-',2),rep('Dark',2), rep('Light',3)) %>% factor
ddsLight=ddsTxiLight[
	rowSums(counts(ddsTxiLight)) > 10 & 
		!names(ddsTxiLight) %in% cont_ids$tid, 
	ddsTxiLight$ctype=='DSLT']
ddsLight$Light=droplevels(ddsLight$Light)
design(ddsLight)=formula(~ Light)
ddsLight=DESeq(ddsLight)
rldLight=rlog(ddsLight,blind=FALSE)
pcaData=plotPCA(rldLight, intgroup = c('TOT','Light'),returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, shape=Light,color=TOT)) +
	geom_point(color = "black", size = 4) +
	geom_point(size=3) +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
	coord_fixed(ratio=1.2) + theme_linedraw() + scale_shape_discrete(name="Light\nexposure") + 
	scale_color_viridis(discrete = TRUE,name="Time\nof\ntreatment")

#Same result for MDS
#mds=as.data.frame(colData(rld)) %>% cbind(cmdscale(dist(t(assay(rld))) %>% as.matrix()))
#ggplot(mds, aes(x = `1`, y = `2`, color = TOT)) +	geom_point(size = 3) + coord_fixed()

resLight=lfcShrink(ddsLight, coef=2, type="apeglm", lfcThreshold = 0.585)

deseq_results[[contrast]]=list(ddstxi=ddsTxiLight,dds=ddsLight,res=resLight)

#Draw MAplot for both contrasts
svg(filename = 'Exploratory_MAplot.svg',width=8,height=4)
par(mfrow=c(1,2))
DESeq2::plotMA(
	resShrink, 
	ylim=c(-1,1)*8,
	xlim=c(0.5,100000),
	cex=.3,
	main='NExp/DSLT',
	sub=expression(paste('logFC > 1, ',italic(s),' < 0.005')))
abline(v=c(1,10,100,1000,10000), col="dodgerblue", lwd=1,lty=3)
abline(h=c(-1,1), col="dodgerblue",  lwd=1,lty=3)
DESeq2::plotMA(
	resLight, 
	xlim=c(0.5,100000),
	ylim=c(-3,3),
	alpha=0.05,
	cex=.3,
	main='DSLT Dark/Light',
	sub=expression(paste('logFC > 0.585, ',italic(s),' < 0.05')))
abline(v=c(1,10,100,1000,10000), col="dodgerblue", lwd=1,lty=3)
abline(h=c(-.585,.585), col="dodgerblue",  lwd=1,lty=3)
dev.off()

library(mixOmics)
library(RColorBrewer)

de = (res$padj < 0.01)
de[is.na(de)] = FALSE
which(de)
rld = rlog(dds, blind = FALSE)
#vsd = vst(dds, blind = FALSE)
plotPCA(rld,intgroup=c("V1"))


library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:50]
df <- as.data.frame(colData(dds)[,c("V1")])
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)


m=assay(rld)[de, ]
#rownames(m)=gsub("MSTRG.","",rownames(m))

cim(m, ylab = "Genes", xlab = "Samples",color = colorRampPalette(brewer.pal(9, "Spectral"))(255), symkey = FALSE, lhei = c(1, 3))


pcols=c("green2","green4","green5","red")
sl=samples$V1 %>% substring(.,1,6) %>% as.factor

pcol=pcols[sl %>% as.numeric]

cim(m, xlab = "Samples",color = colorRampPalette(brewer.pal(9, "Spectral"))(255), symkey = FALSE,row.names=FALSE,cut.tree=c(0.15,0.1),col.sideColors = pcol,legend=list(legend=levels(sl),col=pcols,title="Growth phase"))


cimr=cim(m, xlab = "Samples",color = colorRampPalette(brewer.pal(9, "Spectral"))(255), symkey = FALSE,row.names=FALSE,cut.tree=c(0.15,0.1),col.sideColors = pcol,legend=list(legend=levels(sl),col=pcols,title="Growth phase"))

cim(t(m), ylab = "Samples", xlab="Genes",color = colorRampPalette(brewer.pal(9, "Spectral"))(255), symkey = FALSE,col.names=FALSE,cut.tree=c(0.3,0.20),row.sideColors = pcol,legend=list(legend=c("exponential","lag"),col=c("red3","green4"),title="Growth phase"),save = "DE_1091_Genes_vs_Samples.svg")