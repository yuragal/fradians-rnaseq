#matrosov is mount point of HPC cluster $HOME
dir="~/matrosov/axSA_assemblies/drap/stranded/trinity_norm_altogether/abundance_estimation"

#############################################
#Figure 2a
#############################################
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

#############################################
#Figure 2b
#############################################
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



library(trinotateR)

ann_dir='~/matrosov/axSA_assemblies/drap/stranded/trinity_norm_altogether/trinotate/'
ta=read_trinotate(paste0(ann_dir,'trinotate.annotation_report.nr.csv'))
cont_ids=read.table(paste0(ann_dir,'filter_primates/contaminant.contig.ids'),stringsAsFactors = FALSE)
cont_ids$V2=cont_ids$V1 %>% gsub(pattern = "\\.p[0-9]",replacement = "",.)
names(cont_ids)=c('pid','tid')
ta_clean=ta[!ta$gene_id %in% cont_ids$tid,]
summary_trinotate(ta_clean)

boxplot(ta_clean$transcript %>% nchar,ta_clean$peptide %>% nchar,log='y')
ta_clean$peptide %>% nchar %>% summary(na.rm = TRUE)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   87     200     332     409     511    5205    5450 
ta_clean$peptide %>% nchar %>% sd(na.rm = TRUE)
# 317.8367

#############################################
#Figure S3a
#############################################
#Trinotate results by COG categories
#bash
#zcat ~/matrosov/tools/annotation/Trinotate-v3.1.1/NOG.annotations.tsv.gz > \n
#  ~/matrosov/axSA_assemblies/drap/stranded/trinity_norm_altogether/trinotate/NOG.annotations.tsv

source('./functions.R')
source('./annotation_functions.R')
ids=list(
    Assembly=ta_clean$gene_id, 
    NExp_DSLT=rownames(deseq_results$NExp_DSLT$res)[which(deseq_results$NExp_DSLT$res$svalue<=0.005)],
    Dark_Light=rownames(deseq_results$Dark_Light$res)[which(deseq_results$Dark_Light$res$svalue<=0.05)]
)
p=draw_COG_histogram(triR = ta_clean,ids)
#p+geom_text(data=p$data[p$data$cog=='S',],aes(x=label,y=15,label=freq))

#############################################
#Figure S2a
#############################################
#Additional annotations from TRAPID and Mercator
library(VennDiagram)
library(limma)
#Venn plot
cnts=read.table(paste0(ann_dir,'more_annotations/compare_merged.Mercator_TRAPID_Both.csv'), header=FALSE,stringsAsFactors = FALSE)
rownames(cnts)=cnts$V1
cnts=cnts[,2:4]
colnames(cnts)=c('Mercator','TRAPID','Both')
vp<-venn.diagram(x=list(
    Mercator=rownames(cnts)[which(cnts$Mercator==1)],
    TRAPID=rownames(cnts)[which(cnts$TRAPID==1)],
    Both=rownames(cnts)[which(cnts$Both==1)]),
    cat.fontfamily = "arial",
    fontfamily = "arial",
    cex=0.5, cat.cex=0, lwd=0.5,
    fill=viridis(3,alpha = 0.7),
    #fill = 2:4, alpha = 0.3,
    filename = "Annotation_more_annotations.png",
    units = 'px',width=400, height = 400
    )
#Left -- mercator, Right -- TRAPID, bottom -- Both

#############################################
#Figure S2b
#############################################
#Historgrams of similarities of trinotate and new TRAPID/Mercator annotations
annots=read.table(paste0(ann_dir,'more_annotations/merged.goids.txt'), header=FALSE,stringsAsFactors = FALSE, sep='\t')
colnames(annots)=strsplit('id,blast,pfam,TM',',') %>% unlist
simids=which((annots$blast !='' | annots$pfam != '') & annots$TM != '')
# for 6314 transcripts annotations were updated
annots$triGO=strsplit(paste0(annots$blast,annots$pfam),',')
annots$tmGO=strsplit(annots$TM,',')
library(reshape2)
library(GOSemSim)
source('functions.R')
simDataHs=get_simData(OrgDB = 'org.Hs.eg.db') #Get GO DAG and IC for A.thaliana
onts=strsplit('BP,MF,CC',',') %>% unlist

#it takes about 5-10 minutes each
sims=list(Jiang=sim_by_ontology(annots[simids,c('triGO','tmGO')],simData = simData,measure='Jiang'))
sims$Wang=sim_by_ontology(annots[simids,c('triGO','tmGO')],simData = simData,measure='Wang')
sims$Lin=sim_by_ontology(annots[simids,c('triGO','tmGO')],simData = simData,measure='Lin')
simsdf=melt(sims)[,-3]
simsdf=simsdf[!is.na(simsdf$distance),]
names(simsdf)=c('id','ont','distance','measure')
simsdf$ont %<>% factor(.,levels=onts)
simsdf$ont_measure=paste(simsdf$ont,simsdf$measure,sep=':')
simsdf$label=AddNameStat(simsdf, "ont_measure", "distance", stat = "count")


p1=ggplot(data=simsdf,aes(x=distance,fill=factor(ont,levels=onts))) +
    #geom_freqpoly(binwidth=0.1) +
    #geom_bar(breaks=seq(0,1,by=0.05),stat) +
    geom_histogram(breaks=seq(0,1,by=0.05)) +
    scale_fill_viridis(discrete = TRUE, option = 'plasma', direction = -1,begin=.15,end=0.8,name='Ontology') +
    theme_minimal() +
    facet_wrap(~ label,ncol=3) + 
    #coord_cartesian(ylim=c(0,15)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position =  'NONE') + 
    ylab('Count') + xlab('Similarity bin')
ggsave(p1,filename='Annotation_Semantic_similarity_Trinotate_vs_TRAPID_Mercator.svg',width=6,height=6)


#############################################
#Figure S2c
#############################################
#Histogram of most abundant GO-terms by Trinotate annotation and those added from TRAPID/Mercator
a=annots$triGO %>% unlist
b=annots$tmGO %>% unlist
slimres_tri=get_slim_GOs(a)
slimres_tm=get_slim_GOs(b)
slimres_top[['tri']]=get_slim_topGOs(slimres_tri,min=1)
slimres_top[['tm']]=get_slim_topGOs(slimres_tm,min=1)
onts=strsplit('BP,MF,CC',',') %>% unlist

slimdf=melt(slimres_top)
slimdf=slimdf[slimdf$variable=='Count',][,-2]
colnames(slimdf)[2:4]=c('Count','Ont','Annot')

p1=ggplot(data=slimdf,aes(x=reorder(Term,-Count),y=Count,fill=Annot)) +
    geom_bar(stat='identity') +
    scale_fill_viridis(
        discrete = TRUE,
        option = 'plasma',
        direction = -1,
        begin=.15,
        end=0.8,
        name='Annotation'
    ) +
    facet_wrap(~ Ont,scales='free_x') + 
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        #legend.position =  'NONE',
        plot.margin = margin(10, 10, 10, 90)
    ) + 
    ylab('Count') + xlab('Term')
p1
ggsave(p1,filename='Annotation_GOslim_topGOs.svg',width=15,height=5)

#############################################
#Figure S3b
#############################################
ids=list(
    Assembly=ta_clean$gene_id, 
    NExp_DSLT=rownames(deseq_results$NExp_DSLT$res)[which(deseq_results$NExp_DSLT$res$svalue<=0.005)],
    Dark_Light=rownames(deseq_results$Dark_Light$res)[which(deseq_results$Dark_Light$res$svalue<=0.05)]
)
slimres_top_DE=list()
i=which(annots$gid %in% ids$NExp_DSLT)
a=annots[i,c('triGO','tmGO')] %>% unlist %>% unname
i=which(annots$gid %in% ids$Dark_Light)
b=annots[i,c('triGO','tmGO')] %>% unlist %>% unname
slimres_top_DE[['NExp_DSLT']]=get_slim_GOs(a)
slimres_top_DE[['Dark_Light']]=get_slim_GOs(b)
slimres_top_DE[['NExp_DSLT']]=get_slim_topGOs(slimres_top_DE[['NExp_DSLT']],1)
slimres_top_DE[['Dark_Light']]=get_slim_topGOs(slimres_top_DE[['Dark_Light']],1)

slimdf=melt(slimres_top_DE)
slimdf=slimdf[slimdf$variable=='Count',][,-2]
colnames(slimdf)[2:4]=c('Count','Ont','Contrast')

p1=ggplot(data=slimdf,aes(x=reorder(Term,-Count),y=Count,fill=factor(Contrast,levels=c('NExp_DSLT','Dark_Light')))) +
    geom_bar(stat='identity',position=position_dodge(preserve = 'single')) +
    scale_fill_viridis(discrete = TRUE, option = 'plasma', begin=.15, end=0.8, name='Contrast') +
    facet_wrap(~ Ont,scales='free_x') +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                #legend.position =  'NONE',
                plot.margin = margin(10, 10, 10, 90)) +
    ylab('Count') +
    xlab('Term') +
    scale_y_log10()
p1

ggsave(p1,filename='Annotation_GOslim_topGOs_DE.svg',width=15,height=5)


###########################################################
# Transfer annotations and gene names 
# from phtri, thapse and fracy using BBpHits
# and KEGGREST
# See ~matrosov/axSA_assemblies/drap/stranded/trinity_norm_altogether/trinotate/TD2SAgenome/cmds
# for details on BLAST-analysis pipeline
###########################################################
dann=check_diatom_annots(ta_clean=ta_clean)
View(dann)

#Check correspondence protein-coding transcripts annotated with Trinotate KOs
# and those obtained from blast vs diatoms analysis
idx=which(!is.na(ta_clean$prot_id) & grepl('KO:',ta_clean$Kegg))
ta_kos=ta_clean$Kegg[idx] %>% sapply(.,function(i) unlist(strsplit(i,'KO:'))[2])
names(ta_kos)=ta_clean$prot_id[idx]

idx=which(!is.na(dann$ko))
dia_kos=dann$ko[idx]
names(dia_kos)=dann$prot_id[idx] %>% sub('\\.p\\d$','',.)
names(ta_kos) %in% names(dia_kos) %>% which %>% length

################################
# Trinotate kos -- 6497
# Diatom kos        -- 4283
# !!!Shared kos -- 3535
################################

#Compare KOs for shared transcripts
a=ta_kos[which(names(ta_kos) %in% names(dia_kos))] %>% .[order(names(.),.)]
b=dia_kos[which(names(dia_kos) %in% names(ta_kos))] %>% .[order(names(.),.)]
ko_df=cbind(a,b) %>% as.data.frame(.,stringsAsFactors=FALSE)
ko_df[ko_df$a !=ko_df$b,] %>% dim
rm(a,b,ko_df)
# From 3535 shared KOs, 280 do not match between annoations
# I tend to believe diatom KOs

#Let's try to add 2962 KO-annotations missing in Diatoms and present in Trinotate
#First, check if Trinotate KOs are reasonably assigned
#For this, get NR BBpHs from ta_clean and look if these are mostly belong to protists by taxonomy

a=names(ta_kos)
b=names(dia_kos)
#(!a %in% b) %>% which %>% length
a=a[(!a %in% b) %>% which]
names(a)=ta_kos[(!a %in% b) %>% which]
idx=which(ta_clean$gene_id %in% a)

uann=check_uniprot_annot(ta_clean = ta_clean, idx=idx, ta_kos = ta_kos)
# !!! Turned out 916 transcripts were in both Diatom and SwissProt sets, but have no KO there in dann.
#  In SwissProt set annotations seem more complete for duplicated transcripts
#  Another 80 diatom-annotated transcripts have no KO assigned
#  but KO are present for them in ta_clean
dann=diatom2uniprot_annots(ta_clean=ta_clean,dann=dann,uann=uann)

# !!! 456 transcripts from ta_clean have kegg annotation but their sequences 
#   neither match diatom proteomes 
#   nor they have BBp hits belonging to diatoms or protist against nr

idx=(!is.na(ta_clean$prot_id) & grepl('KO:',ta_clean$Kegg))
idx=which(!ta_clean$prot_id %in% dann$prot_id & idx)
#idx=which(!names(ta_kos) %in% sub('\\.p\\d+','',dann$prot_id))
ann_merged=rbind(
    dann,
    convert_trinotate_records(ta_clean=ta_clean,idx=idx,nrBBpHs = nrBBpHs))


########################################################################
# Map OrthoMCL IDs to transctipts
########################################################################
omclf='~/matrosov/axSA_assemblies/drap/stranded/trinity_norm_altogether/trinotate/TD.orthoMCL.DBv5/orthologGroups'
omcl=read.table(omclf,header = FALSE,stringsAsFactors = FALSE)
t=omcl$V2
names(t)=omcl$V1
omcl=t
omcl[which(omcl=='NO_GROUP')]=NA
omcl=omcl[which(names(omcl) %in% ta_clean$prot_id)]
ta_clean$omcl=NA
idx=ta_clean$prot_id %in% names(omcl)
ta_clean$omcl[idx]=omcl[ta_clean$prot_id[idx]]


save(ta_clean,kos,fkegg_annots,kegg_annots,fkegg_annots2merge,uann,BBpHs,nrBBpHs,nrids,nrids2tx,defs,file,ann_merged,koe,file = 'Trinotate_merge_KEGG.RData')
rm(list=strsplit('idx,idx1,idx2,df,z,dupl,cl,p,p1,s,sheet,go_cl,go_counts,genes,circ,chord,cogs,cont_ids,gos,ns',',') %>% unlist())




