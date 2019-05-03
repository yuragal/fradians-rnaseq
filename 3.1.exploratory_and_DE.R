#matrosov is mount point of HPC cluster $HOME
dir="~/matrosov/axSA_assemblies/drap/stranded/trinity_norm_altogether/abundance_estimation"

#Abundance estimation from Salmon data
# See for more details:
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html


library(dplyr)
library(ggplot2)
library(vsn)
library(DESeq2)
library(magrittr)
library(tximport)
source('./accessory_scripts/functions.R')
#library(tximportData)

#samples.1 is a table
# exp     exp_rep1        /store/ygalachyants/reads/20160805.A-13-1_R1.fq /store/ygalachyants/reads/20160805.A-13-1_R2.fq
# exp     exp_rep2        /store/ygalachyants/reads/20160805.A-13-2_R1.fq /store/ygalachyants/reads/20160805.A-13-2_R2.fq
# lag-0   lag-0_rep1      /store/ygalachyants/reads/20160805.A-3-1_R1.fq  /store/ygalachyants/reads/20160805.A-3-1_R2.fq
# lag-0   lag-0_rep2      /store/ygalachyants/reads/20160805.A-3-2_R1.fq  /store/ygalachyants/reads/20160805.A-3-2_R2.fq
# lag-20  lag-20_rep1     /store/ygalachyants/reads/20160805.A-3-3_R1.fq  /store/ygalachyants/reads/20160805.A-3-3_R2.fq
# lag-20  lag-20_rep2     /store/ygalachyants/reads/20160805.A-3-4_R1.fq  /store/ygalachyants/reads/20160805.A-3-4_R2.fq
# lag-40  lag-40_rep1     /store/ygalachyants/reads/20160805.A-3-6_R1.fq  /store/ygalachyants/reads/20160805.A-3-6_R2.fq

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

ntd=normTransform(dds)
vsd=vst(dds, blind = FALSE)
rld=list(NExp_DSLT=rlog(dds, blind = FALSE))

df <- bind_rows(
    as_data_frame(assay(ntd)[, 1:2]) %>% mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
    as_data_frame(assay(rld[[1]])[, 1:2]) %>% mutate(transformation = "rlog"))
colnames(df)[1:2] <- c("x", "y")
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 50) + coord_fixed() + facet_grid( . ~ transformation)

p1=meanSdPlot(assay(ntd),plot=FALSE ,show.legend=FALSE)$gg + ggtitle('ntd')
p2=meanSdPlot(assay(vsd),plot=FALSE ,show.legend=TRUE)$gg + ggtitle('vsd')
p3=meanSdPlot(assay(rld[[1]]),plot=FALSE ,show.legend=TRUE)$gg + ggtitle('rld')
p=multiplot(p1,p2,p3,cols=3)

#Exploratory analysis plots
#############################################
#Figure 3a
#############################################
library(pheatmap)
#library(RColorBrewer)
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

#############################################
#Figure 3b
#############################################
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
#############################################
#Figure 3c
#############################################
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

rm(p1)
rm(p2)
rm(p3)
rm(ntd)
rm(vsd)
rm(rld)

#############################################
#DE-analysis
#############################################
#See https://support.bioconductor.org/p/101504/ for more info.
#So, 
# 1) get DESeq object (optionally, with some initial filtering)
# 2) perform testing of genes according to experimental design (functions 'results' or 'lfcShrink')
# 3) then choose genes according to adjusted p-value, q-value, s-value etc., depending on testing procedure used
# 
library(apeglm)
contrast='NExp_DSLT'

dds=ddsTxi[ rowSums(counts(ddsTxi)) > 10 & !names(ddsTxi) %in% cont_ids$tid, ]
dds=DESeq(dds)
#!!! Important to have LFC-values consistent with sequecne of contrast levels
#NExp vs DSLT
dds@colData$ctype %<>% as.character() %>% factor(.,levels=c('NExp','DSLT'))
# resShrink %>% summary
# out of 25525 with nonzero total read count
# s-value < 0.005
# LFC > 1.00 (up)    : 2823, 11%
# LFC < -1.00 (down) : 3133, 12%
LFC_treshold=1
res <- lfcShrink(dds, coef=2, type="apeglm", lfcThreshold = LFC_treshold)
deseq_results=list()
deseq_results[[contrast]]=list(ddstxi=ddsTxi,dds=dds,res=res)


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
resLight=lfcShrink(ddsLight, coef=2, type="apeglm", lfcThreshold = 0.585)
deseq_results[[contrast]]=list(ddstxi=ddsTxiLight,dds=ddsLight,res=resLight)

#############################################
#Figure S4
#############################################
#Draw MAplot for both contrasts
svg(filename = 'Exploratory_MAplot.svg',width=8,height=4)
par(mfrow=c(1,2))
DESeq2::plotMA(
    deseq_results[['NExp_DSLT']], 
    ylim=c(-1,1)*8,
    xlim=c(0.5,100000),
    cex=.3,
    main='NExp/DSLT',
    sub=expression(paste('logFC > 1, ',italic(s),' < 0.005')))
abline(v=c(1,10,100,1000,10000), col="dodgerblue", lwd=1,lty=3)
abline(h=c(-1,1), col="dodgerblue",  lwd=1,lty=3)
DESeq2::plotMA(
    deseq_results[['Dark_Light']], 
    xlim=c(0.5,100000),
    ylim=c(-3,3),
    alpha=0.05,
    cex=.3,
    main='DSLT Dark/Light',
    sub=expression(paste('logFC > 0.585, ',italic(s),' < 0.05')))
abline(v=c(1,10,100,1000,10000), col="dodgerblue", lwd=1,lty=3)
abline(h=c(-.585,.585), col="dodgerblue",  lwd=1,lty=3)
dev.off()
