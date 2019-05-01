#If no universe annotations  (uni) and geneNames, look into enrichment.R script
#
#If no resShrink, go to analysis.R

source('functions.R')
onts=c('BP','MF','CC')

#these are stored in DE.Rdata file
load('DE.Rdata')

contrast='NExp_DSLT'

#Load data
dds=deseq_results[[contrast]]$dds
res=deseq_results[[contrast]]$res
topGO_data=topGO_results$NExp_DSLT$data
topGO_res=topGO_results$NExp_DSLT$results


###################################################################################
# Run GSEA with topGO
###################################################################################
geneList=geneNames %in% rownames(subset(res,svalue<=0.005)) %>% as.integer() %>% factor()
names(geneList)=geneNames
#(subset(resShrink,svalue<=0.005) %>% dim(.))[1]
#which(geneList==1) %>% length
#
# NOTICE: there are 5959 DE-genes with LFC>1 and s-value<0.005. 3215 have annotation and subjected to GSEA.

# NOTICE: 125 GO-terms are selected by q-value<0.01
topGO_data=get_topGO_data(gene2GO = uni, genelist = geneList)
topGO_res=get_topGO_res(
	topGO_data,
	alpha=0.01,
	#printDAG = sprintf('./topGO/topGO_%s_',contrast)
	)
#Save to topGO_results object
topGO_results[[contrast]]=list(data=topGO_data,res=topGO_res)

dds@colData$Labels=dds@colData$ctype

###################################################################################
# Draw GO-to-transcript infographics
###################################################################################
#sigGOs -- all enriched terms

sigGOs=plot_counts_by_GO_ggplot(
	tgo=topGO_results$NExp_DSLT,
	dds=dds,
	res=subset(res,svalue<=0.005),
	ret_ggplot = FALSE)

sigGenes=plot_counts_by_GO_ggplot(
	tgo=topGO_results$NExp_DSLT,
	dds=dds,
	subset(res,svalue<=0.005),
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

# For clustering transcripts by function, select GO-terms that have more than 10 transcripts
idx=which(sigGOs$genes
	%>% as.character
	%>% strsplit(.,', ')
	%>% lapply(.,function(l)length(l))
	%>% unlist >= 10)
gos=sigGOs[idx,c(1:3)]
gos$category %<>% as.character()
gos$ID %<>% as.character()
gos$term %<>% as.character()
rownames(gos)=gos$ID

###################################################################################
# Cluster GO-terms
###################################################################################
#Wang distance used here to cluster GO-terms.
#IC-metrics produce NAs in distance matrices as there are GO-terms with Inf Information content

simData=get_simData()
dist_measure = 'Wang'
go_cl=cluster_GOterms(gos = gos, simData = simData, dist_measure = dist_measure)
idx=match(go_cl$dist %>% rownames(),sigGOs$ID)

go_counts=get_counts_up_down_by_GO(
	sigGOs[idx,],
	sigGenes,
	ret_term_names = TRUE)
go_counts$label=
	substr(go_counts$term,1,60) %>% 
	sapply(.,function(i) 
		ifelse(nchar(i)<60,i,paste0(i,'...')))
colnames(go_counts)[(go_counts %>% colnames() %>% grep('up',.))]='NExp Up'
colnames(go_counts)[(go_counts %>% colnames() %>% grep('down',.))]='DSLT Up'

write.xlsx(
	sigGOs[idx,c(1:3)] %>% as.data.frame() %>% cbind(.,go_count_up_down[,3:4]),
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
	sz=8,
	contrast = contrast)

phm_data=list()
phm_data[[contrast]]=draw_GO_term_clusters_heatmap(
	go_cl = go_cl,
	go_counts = go_counts,
	sz=8,
	contrast = contrast,
	silent=TRUE)


###################################################################################
# Draw results of clustering transcripts by GO-terms and LFC (GOCluster)
###################################################################################
# ns=res[geneNames[geneList==1],] %>% as.data.frame()
# genes=data.frame(ID=rownames(ns), logFC=ns$log2FoldChange, adj.P.Val=ns$svalue)
# circ=circle_dat(sigGOs,genes)
# #Fix logFC -- some are NAs !!! there is likely a problem with circle_dat function
# circ$logFC=res$log2FoldChange[match(circ$genes,res %>% rownames %>% toupper)]
# 
# geneids=sigGOs[idx,'genes'] %>%
# 	as.character() %>%
# 	strsplit(.,', ') %>%
# 	unlist %>%
# 	unique %>%
# 	toupper()
# 
# #genes=circ[circ$genes %in% geneids,c(5,6)] %>% unique()
# colnames(genes)=c('ID','logFC')
# rownames(genes)=genes$ID
# 
# #reorder GO-terms as in ann_row for colors from ann_colors to be consistent with intup GOs
# gos=gos[rownames(ann_row),]
# chord=chord_dat(circ,genes,gos %>% rownames)
# #Fix chord logFC values
# #chord[,63]=circ$logFC[match(rownames(chord),circ$genes)]
# #check if chord is consistent
# print(
# 	(
# 		subset(
# 			sigGenes,
# 			(sigGenes$tgs %>% toupper) %in% geneids & 
# 				sigGenes$GO.ID %in% gos$ID
# 			)[,c(2,5,7)] %>% unique() %>% dim)[1] == 
# 	which(chord[,-63]==1) %>% length)
# circ1=circ[circ$ID %in% gos$ID,]
# 
# p1=GOCluster(
# 	data=circ1,
# 	process=gos$ID,
# 	term.col=ann_colors$Group[ann_row$Group %>% as.character()] %>% unname,
# 	lfc.col=viridis::inferno(3,begin = 0.3,end=0.9,direction = -1),
# 	#term.width = 10,
# 	lfc.min=-3,lfc.max=3
# )
# ggsave(p1,filename =sprintf("%s_62GOs_GOCluster2.svg",contrast),width=7,height = 10)

go_colors=phm_data$NExp_DSLT$ann_colors$Group[phm_data$NExp_DSLT$ann_row$Group %>% as.character()] %>% unname
goplot_data=list()
goplot_data[[contrast]]=draw_GO_term_clusters_goplot(
	sigGOs=sigGOs,
	tgos=gos[rownames(phm_data$NExp_DSLT$ann_row),],
	res=res,
	term.col = go_colors,
	ret='data')
p=draw_GO_term_clusters_goplot(
	sigGOs=sigGOs,
	tgos=gos[rownames(phm_data$NExp_DSLT$ann_row),],
	res=res,
	term.col=go_colors,
	lfc.col=viridis::inferno(3,begin = 0.3,end=0.9),
	ret='graph')
ggsave(p,filename =sprintf("%s_62GOs_GOCluster3.svg",contrast),width=15,height=10)
#Save GOplot objects
#goplot_GSEA=list()
goplot_GSEA[[contrast]]=list(
	sigGOs=sigGOs,
	sigGenes=sigGenes,
	gos=gos,
	pheatmap=list(
		go_cl=go_cl,
		ann_row=ann_row,
		ann_col=ann_col,
		ann_colors=ann_colors,
		go_count_up_down=go_count_up_down,
		sz=sz
	),
	circ=circ,
	chord=chord
)


###################################################################################
# This is to draw plots of transcript expression by contrast
# for desired GO-term of gene set
###################################################################################
#p=plot_counts_by_GO_ggplot(tgo=topGO_results$NExp_DSLT,dds=dds,res=resShrink,gos=c('GO:0046961','GO:0009765'))

#GO:0055114	oxidation-reduction process
#GO:0006732	coenzyme metabolic process
#GO:0042182	ketone catabolic process
#GO:0009083	branched-chain amino acid catabolic process
#GO:0015780	nucleotide-sugar transmembrane transport

plot_counts_by_GO_ggplot(tgo=topGO_results$NExp_DSLT,dds=dds,res=res,gos=c('GO:0042182'))

###################################################################################
# Enrichment fof KEGG categories transferred by orthology
###################################################################################
library(edgeR)
library(KEGGREST)
library(gage)
library(DESeq2)

kgs=kegg.gsets(species='ko')
kgs_sigmet=kgs$kg.sets[kgs$sigmet.idx]

koe=list()
attach('DE.Rdata')
contrast="NExp_DSLT"
dds=deseq_results[[contrast]]$dds
res=deseq_results[[contrast]]$res
koe[[contrast]]=prepare_data_for_kegg_analyses(dds,res,ann_merged)

contrast="Dark_Light"
dds=deseq_results[[contrast]]$dds
res=deseq_results[[contrast]]$res
koe[[contrast]]=prepare_data_for_kegg_analyses(dds,res,ann_merged)
detach()

#Here is kegg enrichment commands. 
#KEGG-pathways are weakly enriched, though, only several kos reveal significant response over both contrasts
gage_res=gage(koe[[contrast]]$expr_mat, gsets=kgs_sigmet, ref=c(1:2), samp=c(3:5), compare="as.group")
gage_res$greater %>% head(20)
gage_res$less %>% head(20)



get_transcripts_by_KEGG_pathway(ann_merged,koe$NExp_DSLT,kgs_sigmet,'photosynthesis')

ko04721=list()
contrast="NExp_DSLT"
ko04721[[contrast]]=get_transcripts_by_KEGG_pathway(ann_merged,koe[[contrast]],kgs_sigmet,'Synaptic vesicle cycle',alpha=0.005)
contrast="Dark_Light"
ko04721[[contrast]]=get_transcripts_by_KEGG_pathway(ann_merged,koe[[contrast]],kgs_sigmet,'Synaptic vesicle cycle',alpha=0.05)


get_transcripts_by_KEGG_pathway(ann_merged,koe$NExp_DSLT,kgs_sigmet,'Synaptic vesicle cycle',alpha=0.005)




kids=which( grepl('KO:K',ta_clean$Kegg) & (ta_clean$gene_id %in% rownames(res)) )
#1901 gene with KEGG ortholog mapped
kfc=res$log2FoldChange[ta_clean$gene_id[kids]]
korth=sapply(ta_clean$Kegg[kids],function(i) (strsplit(i,'KO:') %>% unlist)[2]) %>% unname
names(kfc)=korth
#kfc=kfc[which(names(kfc) %in% (kgs$kg.sets %>% unlist %>% unique))]
#allkos=kgs$kg.sets %>% unlist %>% unique
#allkosfc=rep(0,which(!allkos %in% names(kfc)) %>% length)
#names(allkosfc)=allkos[which(!allkos %in% names(kfc))]
#kfc=c(kfc,allkosfc)

gage_res=gage(kfc, gsets=kgs_sigmet, ref=NULL, samp=NULL)


korth_tps=sapply(seq(1,korth %>% length(), by=100), function(i) keggLink('tps',korth[i:(i+99)])) %>% unlist
names(korth_tps) %>% length()
names(korth_tps) %>% unique( ) %>% length()
#Only 553 genes are mapped to fcy and these are from 200 orthogroups
# 402 and 176 mapped to tps
# 347 and 176 mapped to pti
# !!!
# First problem -- foo few orthogroups could be mapped 
# Second problem -- one orthorgoup can contain many genes



###################################################################################
# Export data to supplementary tables
###################################################################################

df=subset(res,svalue<=0.005) %>% as.data.frame()
options("openxlsx.numFmt" = "$* #0,00")
write.xlsx(
	df,
	file='Suppl_Tables.xlsx',
	append=TRUE,
	sheetName=sprintf("%s_transcripts",contrast),
	row.names = TRUE)

df=topGO_res
df$Term=sigGOs$term
df=cbind(
	df,
	get_counts_up_down_by_GO(
		sigGOs,
		sigGenes,
		ret_term_names = TRUE)[,3:4])
colnames(df)[10:11]=c('NExp Up','DSLT Up')

write.xlsx(
	df,
	file='Suppl_Tables.xlsx',
	append=TRUE,
	sheetName=sprintf("%s_GOs",contrast),row.names = FALSE)


options("openxlsx.numFmt" = "$* #,#0.00")
wb=loadWorkbook('Suppl_Tables.xlsx')


