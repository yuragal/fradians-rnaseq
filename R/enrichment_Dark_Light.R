###################################################################################
# README FIRST
###################################################################################
# There is a set of variables in global environment, that stores the objects computed in various DE-analyses.
# Just assign a specific object to desired 'local' variable to start pipeline that uses standard variable names.
# 
# For contrasts NExp/DSLT and DSLT-Dark/Light:
# 1) topGO_results contains topGO data and results
# 2) deseq_results contains ddsTxi, dds, and result objects
# 3) phm_data contains go_cl, go_counts, sz, contrast, ann_row and ann_colors slots
# 4) goplot_data contains sigGOs, sigGenes, gos auxillary objects and circ and chord slots used in GOCluster2 function

#these are stored in DE.Rdata file
load('DE.Rdata')

source('functions.R')
onts=c('BP','MF','CC')
contrast='Dark_Light'

#Construct 'Universe' as a set of annotationts for assembled transcripts
#See https://gist.github.com/slavailn/dcff753cc32fd9e053590308c831b057
uf='~/matrosov/axSA_assemblies/drap/stranded/trinity_norm_altogether/trinotate/trinotate.go_annotations.txt'
uni=readMappings(uf)
geneNames=names(uni)
#map GOs to genes in universe
# l=uni %>% unlist()
# lu=l %>% unique()
# uni_go2gene=lapply(lu,function(i) l[which(l %in% i)] %>% names %>% gsub(perl=TRUE,"(.*)_(\\d)\\d+","\\1_\\2",.))
# names(uni_go2gene)=lu
# rm(list = c('l','lu'))

#Load dds and res DESeq2 objects, topGO slots
dds=deseq_results[[contrast]]$dds
res=deseq_results[[contrast]]$res
topGO_data=topGO_results[[contrast]]$data
topGO_res=topGO_results[[contrast]]$res

###################################################################################
# Run GSEA with topGO
###################################################################################
#Set significance threshold for s-value to 0.05
#!!! Notice !!! Number of genes that are used in GSEA is 406 while there are 858 resLight DE-genes
#i.e. 452 genes do not have functional annotation !!!
geneList=geneNames %in% rownames(subset(res,svalue<=0.05)) %>% as.integer() %>% factor()
names(geneList) <- geneNames
#ddsLight@colData$Labels=factor(c(rep('Dark',2),rep('Light',3)))

topGO_data=get_topGO_data(gene2GO = uni, genelist = geneList)
topGO_res=get_topGO_res(
	topGO_data,
	alpha=0.05,
	#printDAG = sprintf('./topGO/topGO_%s_',contrast)
	)
topGO_results[[contrast]]=list(data=topGO_data,res=topGO_res)

###################################################################################
# Draw GO-to-transcript infographics
###################################################################################
#sigGOs -- all enriched terms with p-val <=0.05

sigGOs=plot_counts_by_GO_ggplot(
	tgo=topGO_results[[contrast]],
	dds=dds,
	res=subset(res,svalue<=0.05),
	ret_ggplot = FALSE)

sigGenes=plot_counts_by_GO_ggplot(
	tgo=topGO_results[[contrast]],
	dds=dds,
	subset(res,svalue<=0.05),
	topN=NULL)
#Check if sigGOs and sigGenes objects are consistent
print(
	(sigGenes$tgs %>% unique() %>% length) == 
		(sigGOs$genes %>% unlist %>% unique %>% length()))

sigGOs=sigGOs[,c(1,2,4,3,5)]
names(sigGOs)=strsplit('category,ID,term,adj_pval,genes',',') %>% unlist
sigGOs$category %<>% factor()
sigGOs$ID %<>% factor()
sigGOs$term %<>% factor()
sigGOs$genes=apply(sigGOs,1,function(r) paste(r$genes,collapse=', ')) %>% factor

################################################################################
# Cluster GO-terms
################################################################################
#For clustering transcripts by function, select GO-terms that have more than 5 transcripts
#p-value threshold was set to 0.05 for Dark/Light contrast
#There are 31 out of 183 GO-terms passing filter
#
#Another variants are 
#		p<=0.05 & transcripts per GO-term >=3 (90 GOs)
#		p<=0.01 & transcripts per GO-term >=3 (43 GOs)
#		
#		In both variants, manual group assignment was performed. However, clustering terms on similarity is weak.
#		Manually assigned groups do not bundle into consecutive colored blocks in the annotation track.
#		That is why I tend to leave the first variant of analysis where 31 GO-terms are clustered
idx=which(
	sigGOs$genes
		%>% as.character
		%>% strsplit(.,', ')
		%>% lapply(.,function(l)length(l))
		%>% unlist >= 6)
		#%>% unlist >= 3 & sigGOs$adj_pval <=0.01)
gos=sigGOs[idx,c(1:3)]
gos$category %<>% as.character()
gos$ID %<>% as.character()
gos$term %<>% as.character()
rownames(gos)=gos$ID

#Jiang distance used here to cluster GO-terms.
dist_measure = 'Jiang'
go_cl=cluster_GOterms(gos = gos, simData = simData, dist_measure = dist_measure)
go_cl$dist_measure=dist_measure
idx=match(go_cl$dist %>% rownames(),sigGOs$ID)
go_counts=get_counts_up_down_by_GO(
	sigGOs[idx,],
	sigGenes,
	ret_term_names = TRUE)
go_counts$label=
	substr(go_counts$term,1,60) %>% 
	sapply(.,function(i) 
		ifelse(nchar(i)<60,i,paste0(i,'...')))
colnames(go_counts)[(go_counts %>% colnames() %>% grep('up',.))]='Dark Up'
colnames(go_counts)[(go_counts %>% colnames() %>% grep('down',.))]='Light Up'

# This is to for classifying larger sets of GO-terms (GO43 and GO90)
#
# df=data.frame(group=phm_data$Dark_Light$ann_row$Group %>% as.character(),ID=phm_data$Dark_Light$go_counts$id)
# df1=data.frame(sigGOs[idx,c(1:3)], go_counts[,3:4])
# df1$ID %<>% as.character()
# df1$group=df$group[match(df1$ID,df$ID)]

# group90=read.xlsx(
# 	file='Suppl_Tables.xlsx',
# 	sheetName=sprintf("%s_GOCluster_90",contrast))[,c(2,6)]
# 
# group43=read.xlsx(
# 	file='Suppl_Tables.xlsx',
# 	sheetName=sprintf("%s_GOCluster_43",contrast))

write.xlsx(
	#df1,
	sigGOs[idx,c(1:3)] %>% as.data.frame() %>% cbind(.,go_counts[,3:4]),
	file='Suppl_Tables.xlsx',
	append=TRUE,
	sheetName=sprintf("%s_GOCluster",contrast),
	row.names = FALSE)
# Manually classify GO-terms to some feasible groups
# The main idea is to decrease the number of categories to some number < 20.
# Then name column as 'group', save table and read it back to R.

#Draw heatmap, save annotations of heatmap to list phm_data
#Export image size is 1500x800px
draw_GO_term_clusters_heatmap(
	go_cl = go_cl,
	go_counts = go_counts,
	sz=10,
	contrast = contrast,
	#sheetName = sprintf("%s_GOCluster_43",contrast)
	)
phm_data[[contrast]]=draw_GO_term_clusters_heatmap(
	go_cl = go_cl,
	go_counts = go_counts,
	sz=10,
	contrast = contrast,
	silent=TRUE)


###################################################################################
# Draw results of clustering transcripts by GO-terms and LFC (GOCluster)
###################################################################################
go_colors=
	phm_data[[contrast]]$ann_colors$Group[
		phm_data[[contrast]]$ann_row$Group %>% 
			as.character()] %>% unname

p=draw_GO_term_clusters_goplot(
	sigGOs=sigGOs,
	tgos=gos[rownames(phm_data[[contrast]]$ann_row),],
	res=res,
	term.col=go_colors,
	lfc.col=viridis::inferno(3,begin = 0.3,end=0.9,direction = -1),
	lfc=1.5,
	ret='graph')
ggsave(p,filename=sprintf("%s_31GOs_GOCluster2.svg",contrast),width=15,height=10)

goplot_data[[contrast]]$sigGOs=sigGOs
goplot_data[[contrast]]$sigGenes=sigGenes
goplot_data[[contrast]]$gos=gos

goplot_data[[contrast]][[c('circ','chord')]]=draw_GO_term_clusters_goplot(
	sigGOs=sigGOs,
	tgos=sigGOs,
	res=res,
	term.col = NULL,
	ret='data')

################################################################################
# GO-chord for a subset of enriched GO-terms
################################################################################
# Now select a subset of higly-signigicant GO-terms for which chord-graph is plotted
idx=which(sigGOs$adj_pval<=0.0005)
# save them to goshs* to match them later with a set of GO-terms
# of lower significance (gos* were plotted in GOCluster graph)
gos=sigGOs[idx,c(1:3)]
gos$category %<>% as.character()
gos$ID %<>% as.character()
gos$term %<>% as.character()
rownames(gos)=gos$ID

#cluster GO-terms
dist_measure = 'Wang'
go_cl=cluster_GOterms(gos = gos, simData = simData,dist_measure=dist_measure)
# Warning message:
# 	GO-terms GO:0102985 not present in database!
# 	Assign zeroes for these terms in the resulted distance matrix.

idx=match(go_cl$dist %>% rownames(),sigGOs$ID)
go_counts=get_counts_up_down_by_GO(
	sigGOs[idx,],
	sigGenes,
	ret_term_names = TRUE)
go_counts$label=
	substr(go_counts$term,1,60) %>% 
	sapply(.,function(i) 
		ifelse(nchar(i)<60,i,paste0(i,'...')))
colnames(go_counts)[(go_counts %>% colnames() %>% grep('up',.))]='Dark Up'
colnames(go_counts)[(go_counts %>% colnames() %>% grep('down',.))]='Light Up'

#Draw heatmap, save annotations of heatmap to list phm_data
#Export image size is 1500x800px
draw_GO_term_clusters_heatmap(
	go_cl = go_cl,
	go_counts = go_counts,
	sz=10,
	contrast = contrast,
	sheetName = sprintf("%s_GOChord",contrast)
)
phm_data[[paste0(contrast,'_GOChord')]]=draw_GO_term_clusters_heatmap(
	go_cl = go_cl,
	go_counts = go_counts,
	sz=10,
	contrast = contrast,
	sheetName = sprintf("%s_GOChord",contrast),
	silent=TRUE)

#chg -- chord groups (14 GO-terms)
#chl -- chord levels (14 GO-terms)
#clg -- cluster groups (31 GO-term)
#cll -- cluster levels (31 GO-terms)
chg=phm_data$Dark_Light_GOChord$ann_row$Group %>% as.character
chl=phm_data$Dark_Light_GOChord$ann_row$Group %>% levels
clg=phm_data$Dark_Light$ann_row$Group %>% as.character
cll=phm_data$Dark_Light$ann_row$Group %>% levels
idx=match(chl,cll) %>% na.omit()
#matched colors
mc=phm_data$Dark_Light$ann_colors$Group[cll[idx]]
st=length(phm_data$Dark_Light$ann_colors$Group) + 4
colvec=pals::kelly()[st:(st + length(chl) - length(mc) - 1)] %>% rev()
names(colvec)=chl[attr(idx,'na.action')]
colvec=c(colvec,mc)
phm_data$Dark_Light_GOChord$ann_colors$Group=colvec[chg %>% unique]
rm(list = strsplit('chg,chl,clg,cll,idx,mc,st,colvec',',') %>% unlist)

pheatmap(
	go_cl$dist,
	cluster_rows = FALSE,
	cluster_cols = FALSE,
	color=viridis::viridis(255,begin = 0.25),
	labels_row = phm_data$Dark_Light_GOChord$go_counts$label,
	annotation_row = phm_data$Dark_Light_GOChord$ann_row,
	annotation_colors=phm_data$Dark_Light_GOChord$ann_colors,
	cellwidth = 10,cellheight = 10
)

# Now make GOChord graph
circ=goplot_data$Dark_Light$circ
#geneids -- subset of DE-genes associated with enriched GO-terms
gos=gos[go_cl$dist %>% rownames,]
geneids=sigGOs[match(gos$ID,sigGOs$ID),'genes'] %>%
	as.character() %>%
	strsplit(.,', ') %>%
	unlist %>%
	unique %>%
	toupper()
genes=circ[circ$genes %in% geneids,c(5,6)] %>% unique()
rownames(genes)=genes$ID
colnames(genes)=c('ID','logFC')

#Create chord object taking order of GO-terms from clustering results
chord=chord_dat(circ,genes,go_cl$dist %>% rownames)
#Fix transcript order (group one transcript with set of similar ones)
chord=chord[c(1:21,45,22:44),]

p=GOChord(
	chord, 
	gene.space=0.2,
	gene.size=3,
	process.label=5,
	lfc.col=viridis::inferno(3,begin = 0.3,end=0.9,direction = -1),
	border.size=0.2,space=0.01,
	lfc.min=-1.5,lfc.max = 1.5
	#limit = c(1,3)
)
#save for annotation of bubble plot 
chord_gos=p$layers[[3]]$data$labels

#Transcript names
p$layers[[4]]$data$labels %<>% 
	gsub('TRINITY_','',.) %>%
	gsub('CONTIG.*','',.) %>%
	gsub('_.*','',.)

f=p$layers[[5]]$aes_params$fill
m=match(f,f %>% unique())
n=f %>% unique() %>% length()
#c=viridis(n,alpha = 0.7)
#Brewer palettes are not suitable when number of GOs more than 12
#c=brewer.pal(n,name='Accent')
#c=brewer.pal(n,name='Set1')
c=add.alpha(phm_data$Dark_Light_GOChord$ann_colors$Group[phm_data$Dark_Light_GOChord$ann_row$Group %>% as.character],0.9)

#color of GO terms guides
p$guides[[1]]$override.aes$fill=c
#color of GO terms circle segments
r=which(p$layers[[2]]$data$id==1) %>% length()
p$layers[[2]]$aes_params$fill=sapply(c,function(i) rep(i,r)) %>% as.vector()
#color of ribbons connecting GO terms and genes
p$layers[[5]]$aes_params$fill=sapply(m, function(i) c[i])

df=p$layer[[2]]$data
df=sapply(
	1:length(c),
	function(i) 
		c(
			mx=mean(df$x[df$id==i]),
			my=mean(df$y[df$id==i]))
	) %>% 
	t %>% 
	data.frame(.)
df$id=rownames(df)
df=cbind(df,group=phm_data$Dark_Light_GOChord$ann_row$Group, label=go_counts$label)

p1=p+geom_text(data=df,aes(x=mx,y=my,label=id))
require(ggrepel)
p1=p1+geom_text_repel(
	data=df[df$my>=0,],aes(x=mx,y=my,label=group),
	nudge_x       = 10 - df$my[df$my>=0],
	segment.size  = 0.3,
	segment.color = "gray30",
	direction     = "x",
	size          = 2
	#vjust         = 0
	)
p1=p1+geom_text_repel(
	data=df[df$my<0,],aes(x=mx,y=my,label=group),
	nudge_x       = 10 - df$my[df$my<0],
	segment.size  = 0.3,
	segment.color = "gray30",
	direction     = "x",
	size          = 2
	#vjust         = 0
)
#draw
p1
#save as 800x1000 svg

###################################################################################
# Draw GOBubble plot
###################################################################################
#GOBubble(circ, display = 'multiple', colour=viridis(3), bg.col = T, labels = 3)
#reduced_circ <- reduce_overlap(circ, overlap = 0.75)
p=GOBubble(circ, colour=colvec.onts,display = 'multiple', labels = 3)
p1=ggplot(p$data,aes(x=zscore,y=adj_pval,label=id,fill=category,size=count))
p1=p1 + geom_point(shape=21,stroke=0.2)+theme_light()
p1=p1 + scale_size(trans='log10',range=c(1,12),breaks=c(10,100),name='Abundance')
p1=p1 + scale_fill_manual(values=add.alpha(colvec.onts,0.7)) #scale_fill_viridis(discrete=TRUE,alpha = 0.4)
p1=p1 + guides(fill=guide_legend(title='Ontology'))
#,labels=c('Biological process','Molecular function','Cellular component'))

anndf=p$data[which(p1$data$adj_pval>3),]
#!is.na(match(p1$data$term,chord_gos))),]

anndf$id[match(chord_gos,anndf$id)]=1:length(chord_gos)
p1=p1+geom_text(data=anndf[grep('^\\d',anndf$id),],aes(x=zscore,y=adj_pval,label=id,size=1.5)) 
#p1=p1+geom_text_repel(data=anndf[grep('^GO',anndf$id),],aes(x=zscore,y=adj_pval,label=id,size=1.2))

#anndf=p$data[which(!is.na(match(p1$data$term,chord_gos))),]
#anndf$term=1:length(chord_gos)
p1=p1+labs(x='Z-score',y=expression(paste('-log'[10],'(adjusted ',italic('p'),'-value)')))
#p1+facet_wrap(~category)

p1
ggsave(p1,filename = 'DSLT_GO_Bubbleplot_plasma1.svg',width=7,height = 6)


###################################################################################
# Export data to supplementary tables
###################################################################################
df=subset(res,svalue<=ifelse(contrast=='Dark_Light',0.05,0.005)) %>% as.data.frame()
df$baseMean %<>% sprintf("%.3f",.) %>% as.numeric()
df$log2FoldChange %<>% sprintf("%.3f",.) %>% as.numeric()
df$lfcSE %<>% sprintf("%.3f",.) %>% as.numeric()
df$svalue %<>% sprintf("%.3e",.) %>% as.numeric()
write.xlsx(
	df,
	file='Suppl_Tables.xlsx',
	append=TRUE,
	sheetName=sprintf("%s_DE_transcripts",contrast),
	row.names = TRUE)

df=cbind(
	topGO_res,
	get_counts_up_down_by_GO(
		sigGOs,
		sigGenes,
		ret_term_names = TRUE)[,3:4])
colnames(df)[10:11]=c('Dark Up','Light Up')
df$Term=sigGOs[match(df$GO.ID,sigGOs$ID),'term'] %>% as.character()
write.xlsx(
	df,
	file='Suppl_Tables.xlsx',
	append=TRUE,
	sheetName=sprintf("%s_GSEA_GOs",contrast),row.names = FALSE)

df=sigGenes[,c(2,5)] %>% unique
colnames(df)=c('GO','transcript')
df=df[,2:1]
write.xlsx(
	df,
	file='Suppl_Tables.xlsx',
	append=TRUE,
	sheetName=sprintf("%s_transcript2GO",contrast),row.names = FALSE)

###################################################################################
# This is to draw plots of transcript expression by contrast
# for desired GO-term of gene set
###################################################################################

df=cbind(
	topGO_res,
	get_counts_up_down_by_GO(
		sigGOs,
		sigGenes,
		ret_term_names = TRUE)[,3:4])
colnames(df)[10:11]=c('Dark Up','Light Up')
ldf=lapply(onts,function(o) subset(df, Ontology==o)[1:10,c(2,3,9:11)])

# BP.5  GO:1903861   positive regulation of dendrite extension 4.6e-04       2        1
# BP.6  GO:0009867    jasmonic acid mediated signaling pathway 5.7e-04       2        2
# BP.7  GO:0009657                        plastid organization 7.1e-04       4        6
# BP.8  GO:0006429                  leucyl-tRNA aminoacylation 7.9e-04       1        2
# MF.8  GO:0004823                leucine-tRNA ligase activity 8.3e-04       1        2
# CC.1  GO:0009524                                phragmoplast 0.00014       3        2
# CC.5  GO:0005886                             plasma membrane 0.00123      48       33
# CC.6  GO:0016021              integral component of membrane 0.00170      58       62
# CC.7  GO:0009706                  chloroplast inner membrane 0.00342       1        4


p=plot_counts_by_GO_ggplot(
	svalue=0.05,
	tgo=topGO_results$Dark_Light,
	dds=deseq_results$Dark_Light$dds,
	res=deseq_results$Dark_Light$res,
	gos=strsplit('GO:1903861,GO:0009867,GO:0009657,GO:0006429,GO:0004823',',') %>% unlist)
ggsave(p,filename = sprintf('%s_bilateral_GOs1.svg',contrast),height = 7,width=12)
p=plot_counts_by_GO_ggplot(
	svalue=0.05,
	tgo=topGO_results$Dark_Light,
	dds=deseq_results$Dark_Light$dds,
	res=deseq_results$Dark_Light$res,
	gos=strsplit('GO:0009524,GO:0005886,GO:0016021,GO:0009706',',') %>% unlist)
ggsave(p,filename = sprintf('%s_bilateral_GOs2.svg',contrast),height = 7,width=12)



#dsCYC2 gene is encoded by CL4842Contig1_1. 
#> which(assay(deseq_results$NExp_DSLT$dds) %>% rownames %in% 'CL4842Contig1_1')
# [1] 4526
# > assay(deseq_results$NExp_DSLT$dds)[4526,]
# NExp1    NExp2  DSLT1-0  DSLT2-0 DSLT1-20 DSLT2-20 DSLT1-40 
# 901      407      990      862     1184      774      536 
# > assay(deseq_results$Dark_Light$dds)[4526,]
# DSLT1-0  DSLT2-0 DSLT1-20 DSLT2-20 DSLT1-40 
# 990      862     1184      774      536 
# > deseq_results$Dark_Light$res[4526,]
# log2 fold change (MAP): Light Light vs Dark 
# 
# DataFrame with 1 row and 4 columns
# baseMean     log2FoldChange             lfcSE            svalue
# <numeric>          <numeric>         <numeric>         <numeric>
# 	CL4842Contig1_1 876.998168516948 -0.400626495919575 0.297853413258469 0.397567280852965
# > deseq_results$NExp_DSLT$res[4526,]
# log2 fold change (MAP): ctype NExp vs DSLT 
# 
# DataFrame with 1 row and 4 columns
# baseMean   log2FoldChange             lfcSE             svalue
# <numeric>        <numeric>         <numeric>          <numeric>
# 	CL4842Contig1_1 861.757790010396 1.31799036577506 0.398352157300559 0.0410845984503907