require(xlsx)
require(topGO)
require(magrittr)
require(DESeq2)
require(GOplot)
require(AnnotationDbi)
require(pheatmap)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	library(grid)
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	
	numPlots = length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
										 ncol = cols, nrow = ceiling(numPlots/cols))
	}
	
	if (numPlots==1) {
		print(plots[[1]])
		
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		
		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
																			layout.pos.col = matchidx$col))
		}
	}
}

split_eggnog <- function(x, hit = "eggnog"){
	y <- x[!is.na( get(hit) ), .( get(hit), gene_id, transcript_id, prot_id) ]
	
	# split multiple annotations in backtick-delimited list
	z <- strsplit(y$V1, "`" )
	n <- sapply(z, length)
	# and then split GO annotation into 3 columns
	
	z <- strsplit( unlist(z), "\\^" )
	
	x1 <- data.frame(  gene = rep(y$gene_id, n), 
										 transcript = rep(y$transcript_id, n) ,  
										 protein = rep(gsub(".*\\|", "", y$prot_id), n),
										 egg = sapply(z, "[", 1) , 
										 description = sapply(z, "[", 2) , 
										 stringsAsFactors=FALSE )
	message(nrow(x1), " ", hit, " annotations")
	data.table(x1)
}

# Function is to plot counts of genes belongong to specific GO-term by contrast variable.
# This is very ineffficient for large gene datasets
# Use plot_counts_by_GO_ggplot for this purpose
plot_counts_by_GO <- function(go=NULL,uni,dds,res,ns=NULL,filename=NULL){
	require(magrittr)
	require(DESeq2)
	if(!is.null(go)){
		ids=grep(go,uni)
		ddsns=dds %>% rownames
		ns=names(uni[ids])
	}
	#Check if all transcript IDs are in target dataset. If no, drop ID.
	ns=ddsns[which(ddsns %in% ns)]
	#Sort IDs by s-value ascending
	ns=ns[order(res[ns,'svalue'])]
	if(length(ns)>10){
		warning(sprintf("The raw list of genes to plot for %s is too long.\n  The first 10 genes with smallest s-values will be plotted.",go))
		ns=ns[1:10]
	}
	if(is.null(ns)){
		message('List of genes to plot is empty.')
		return()
	}
	
	nrows=ifelse(length(ns)<5,1,ceiling(length(ns)/5))
	if(!is.null(filename)){
		svg(filename = filename,width=10,height=nrows*3)
	}
	par(mfrow=c(nrows,5))
	
	#0.5 is a pseudocount used in DESeq2::plotCounts function. 
	#To have y-axis ticks, y-limits, counts, and count means consistent, 
	#I had to introduce additional 0.5 pseudocount.
	#For generated graphs, the resulted pseudocount is 1.
	assay(dds[ns,])=dds[ns,] %>% assay %>% ifelse(.==0,0.5,.)
	
	lys=cbind(
		'y1'=(dds[ns,] %>% assay)[,1:2] %>% rowSums/2+0.5,
		'y2'=(dds[ns,] %>% assay)[,3:5] %>% rowSums/3+0.5
	)
	lims=c(dds[ns,] %>% assay %>% min*.8,dds[ns,] %>% assay %>% max*1.2)
	lims=lims %>% ifelse(.<=0.5,1,.)
	
	for(n in ns){
		sval=res[n,'svalue']
		plotCounts(
			dds, n, 'Labels', transform = TRUE, normalize = FALSE,
			xlab = '',
			ylim = lims,
			main = paste0(
				LETTERS[which(ns %in% n)],
				ifelse(
					sval<=0.005,'***',ifelse(
						sval<=0.01,'**',ifelse(
							sval<=0.05,'*','')))
			)
		)
		title(
			sub=sprintf('%s\ns-value: %1.1e',
				n %>% gsub('_','',.) %>% gsub('TRINITY','',.) %>% gsub('Contig','Cnt',.),
				sval
				)
		)
		segments(x0=1,y0=lys[n,1],x1=2,y1=lys[n,2],col='red')
	}
	if(!is.null(filename)){
		dev.off()
	}
	par(mfrow=c(1,1))
}

plot_counts_by_GO_ggplot <- function(tgo,dds,res,gns=NULL,gos=NULL,filename=NULL,topN=10,svalue=NULL,ret_ggplot=TRUE){
	require(magrittr)
	require(DESeq2)
	require(ggplot2)
	require(tidyr)
	
	# require(AnnotationDbi)
	# require(topGO)
	
	#tgo -- is a list with slots 
	#	'data' -- list of three topGO_data ojects and 
	#	'res'  -- dataframe with enrched  GO -terms
	#	
	#dds and res -- are DESeq2 DESeqDataSet and DESeqResult objects
	#res is subsetted to transcripts for which plots should be drawn
	# First, create dataframe to plot
	sigGOs=tgo$res[,c(1,2,9)]
	colnames(sigGOs)[1]='category'
	sigGOs$term=AnnotationDbi::Term(sigGOs$GO.ID)

	#for each GO-term get genes
	#order genes by s-value
	idx=order(res$svalue)
	sigGOs$genes=apply(sigGOs, 1, function(go,gns) {
		t=topGO::genesInTerm(tgo$data[[ go[[1]] ]], go[[2]]) %>% 
			unlist %>%
			unname
		t=gns[match(t,gns) %>% na.omit() %>% sort()]
		return(t)
		},gns=rownames(res)[idx])
	
	if(!ret_ggplot){
		return(sigGOs)
	}
	
	if(!is.null(gos)){
		sigGOs=subset(sigGOs,GO.ID %in% gos)
	}
	
	if(!is.null(gns)){
		sigGOs=apply(sigGOs,1,function(go,gns) {
			idx=which(go$genes %in% gns)
			if(length(idx)==0){
				return()
			}else{
				go$genes=go$genes[idx]
				return(go)
				}
			},gns=gns)
	}
	if(!is.null(topN)){
		sigGOs$tgs=lapply(sigGOs$genes,function(genes,topN) {
			#subset topN genes if list is long
			if(length(genes)>topN){
				genes=genes[1:topN]
				}
			return(genes)
			},topN=topN)
		} else {
			sigGOs$tgs=sigGOs$genes
			}
	
	df=unnest(sigGOs,tgs)
	df$svalue=res[df$tgs,]$svalue
	#subset by s-value
	if(!is.null(svalue)){
		df=df[df$svalue<=svalue,]
	}
	df$LFC=res[df$tgs,]$log2FoldChange
	df=cbind(df,assay(dds)[df$tgs,]+1)
	df=reshape2::melt(df,id=c(1:7))
	colnames(df)[9]='count'
	contrast=as.character(dds@design)[2]
	if(contrast=='ctype'){
		df$contrast=gsub('-?\\d+','',df$variable)
	} else if(contrast=='Light'){
		df$contrast=gsub('DSLT\\d-0','Dark',df$variable)
		df$contrast=gsub('DSLT\\d-\\d0','Light',df$contrast)
	} else if(contrast=='salinity'){
	  df$contrast=dds@colData$salinity[df$variable] %>% as.character
	} else{
		message("Requested contrast cannot be handled. Please, write the contrast processing porcedure before move on.")
		return()
	}
	df$contrast=factor(
		df$contrast,
		levels=levels(dds@colData[[contrast]]))
	
	# df$tr_title=paste(ifelse(
	# 		df$svalue<=0.005,'***',ifelse(
	# 			.<=0.01,'**',ifelse(
	# 				.<=0.05,'*',''))),'\n')
	df$tr_title=paste(df$tr_title,
										paste(df$tgs %>% 
												gsub('_','',.) %>% 
												gsub('TRINITY','',.) %>% 
												gsub('Contig','Cnt',.),
											sprintf('\ns=%1.1e',df$svalue)
										))
	
	#This is to order facets by s-value and input GO-terms (if any)
	df$tr_title=factor(df$tr_title,levels=df$tr_title %>% unique())
	if(!is.null(gos)){
		df$GO.ID=factor(df$GO.ID,levels=gos)
	}
	
	if(is.null(gos) & is.null(gns)){
		message("No gene names or GO-terms are provided for plotting. Return the prepared dataframe.")
		return(df)
	}
	
	# Now plot with facet_grid ggdplot2 function
	p=ggplot(df,aes(x=contrast,y=count,group=GO.ID)) +
		geom_jitter(size=1.5,position = position_jitter(width=0.1))+
		#ggtitle(label=df$GO_title,subtitle = df$tr_title) +
		facet_wrap(GO.ID ~ tr_title,ncol = 10) +
		scale_y_log10() +
		stat_summary(fun.y=mean, geom="line", colour="red", size=0.5) +
		xlab(contrast) + theme_minimal()
		
	return(p)
}


#Calculate up- and down-regulated genes by GO-terms
#Input is data frame with rows as in siGOs
get_counts_up_down_by_GO <- function(sigGOs,sigGenes,ret_term_names=FALSE){
	require(magrittr)
	
	df=apply(sigGOs,1,function(go,genes) {
		id=go[2] %>% as.character()
		term=go[3] %>% as.character()
		gns=go[5] %>% as.character() %>% strsplit(.,', ') %>% unlist
		#cat(id, sep='\n')
		df=subset(genes, tgs %in% gns & GO.ID == id, c(2,5:7)) %>% unique()
		lfc=df$LFC > 0
		#Return number of genes with LFC<0 and LFC>0
		#i.e. first are genes upregulated in levels(dds@colData[[ dds$design[2] ]])[1]
		#  second are genes upregulated in   levels(dds@colData[[ dds$design[2] ]])[2]
		lfc=c(length(lfc)-length(which(lfc)),length(which(lfc)))
		return(c(id,term,lfc))
	},genes=sigGenes)
	df %<>% t %>% as.data.frame()
	colnames(df)=strsplit('id,term,up,down',',')  %>% unlist
	df$up %<>% as.character %>% as.numeric()
	df$down %<>% as.character %>% as.numeric()
	
	if(!ret_term_names){
		df=df[,c(1,3,4)]
	}
	return(df)
}

get_gids_GO_enriched <- function(go=NULL,uni,dds,res,ns=NULL,filename=NULL){
	if(!is.null(go)){
		ids=grep(go,uni)
		ddsns=dds %>% rownames
		ns=names(uni[ids])
	}
	#Check if all transcript IDs are in target dataset. If no, drop ID.
	ns=ddsns[which(ddsns %in% ns)]
	#Sort IDs ascending by s-value
	ns=ns[order(res[ns,'svalue'])]
	if(length(ns)>10){
		warning(sprintf("The raw list of genes to plot for %s is too long.\n  The first 10 genes with smallest s-values will be plotted.",go))
		ns=ns[1:10]
	}
	if(is.null(ns)){
		message('List of genes to plot is empty.')
		return()
	}
}

get_topGO_data <- function(gene2GO,genelist,onts='BP,MF,CC') {
	require(topGO)
	require(magrittr)
	onts=strsplit(onts,',') %>% unlist
	d=lapply(onts, function(ont) {
		GOdata=new(
			"topGOdata",
			annot = annFUN.gene2GO, 
			ontology = ont, 
			allGenes = geneList,
			gene2GO = uni)
		return(GOdata)
		})
	names(d)=onts
	return(d)
}

get_topGO_res <- function(GOdata,alpha,printDAG=NULL){
	require(topGO)
	require(magrittr)
	
	allRes=lapply(GOdata, function(d) {
		res.classic = runTest(d, algorithm = "classic", statistic = "Fisher")
		res.elim = runTest(d, algorithm = "elim", statistic = "Fisher")
		
		if(!is.null(printDAG)){
			printGraph(
				d,
				res.elim,
				firstSigNodes = 5,
				fn.prefix = paste0(printDAG,d@ontology),
				useInfo = "all",
				pdfSW = TRUE)
			
			printGraph(
				d,
				res.elim,
				firstSigNodes = 10,
				fn.prefix = paste0(printDAG,d@ontology),
				useInfo = "all",
				pdfSW = TRUE)
		}
		t=GenTable(
			d, 
			classic = res.classic,
			elim = res.elim,
			orderBy = "elim",
			topNodes = 100)
		return(t)
	})
	
	allRes=do.call('rbind',allRes)
	allRes$classic %<>% as.numeric()
	allRes$elim %<>% as.numeric()
	allRes$Ontology=rownames(allRes) %>% gsub('\\.\\d+$','',.)
	allRes=allRes[,c(9,1:8)]
	allRes=subset(allRes,elim<=alpha)
	return(allRes)
}

add.alpha <- function(col, alpha=1){
	if(missing(col))
		stop("Please provide a vector of colours.")
	apply(sapply(col, col2rgb)/255, 2, 
				function(x) 
					rgb(x[1], x[2], x[3], alpha=alpha))  
}

get_simData <- function(OrgDB='org.At.tair.db',IC=TRUE){
	library(GOSemSim)
	onts=strsplit('BP,MF,CC',',') %>% unlist
	simData=lapply(onts, function(ont) godata(OrgDb = OrgDB, ont=ont,computeIC = IC ))
	names(simData)=onts
	return(simData)
}

cluster_GOterms <- function(gos,simData,dist_measure='Wang'){
	library(pheatmap)
	onts=strsplit('BP,MF,CC',',') %>% unlist
	#go_dist=lapply(onts, function(ont){
	go_dist=list()
	for(ont in onts){
		d=simData[[ont]]
		idx=which(gos$category == ont)
		idx1=which(gos$ID %in% (d@IC %>% names))
		
		t=mgoSim(
			gos$ID[idx1],
			gos$ID[idx1],
			semData=d,
			measure=dist_measure,
			combine=NULL
			)
		
		gos_na=which(is.na(t[,1]))
		if(length(gos_na)>0){
			warning(sprintf(
				"Similarity of GO-terms %s is NA!\n Most likely, it means IC of these GO-terms is -Inf.\n You may wich to use IC-free distance metric for you GO-terms set.\n Forhierarchial clustering, assign -0.1 for these terms in the resulted distance matrix.",
				paste0(rownames(t)[gos_na],sep=',')))
			t[which(is.na(t))]=-0.1
		}
		
		if(dim(t)[1]!=length(idx)){
			g1=gos$ID[ which(gos$category == ont) ]
			g2=rownames(t)
			miss_ids=which(match(g1,g2) %>% is.na)

			warning(sprintf(
				"GO-terms %s not present in database!\n Assign -0.1 for these terms in the resulted distance matrix.",
				paste(g1[miss_ids],sep=', ',collapse='')))
			t1=rep(-0.1,length(g1)^2)
			dim(t1)=rep(length(g1),2)
			dimnames(t1)=list(g1,g1)
			# Names are temporary. These will be updated after clustering GO-terms.
			#diag(t1)=1
			l=dim(t)[1]
			t1[1:l,1:l]=t
			rownames(t1)[1:l]=rownames(t)
			colnames(t1)[1:l]=colnames(t)
			for(i in 1:length(miss_ids)){
				rownames(t1)[l+i]=g1[ miss_ids[i] ]
				colnames(t1)[l+i]=g1[ miss_ids[i] ]
			}
			t=t1
		}
		go_dist[[ont]]=t
	#	return(t)
	#})
	}
	names(go_dist)=onts
	#Init go_dist1 -- Matrix of clustered GO distances
	go_dist1=rep(NA,length(gos$ID)^2)
	dim(go_dist1)=rep(length(gos$ID),2)
	dimnames(go_dist1)=list(gos$ID,gos$ID) 
	# Names are temporary. These will be updated after clustering GO-terms.
	# Now cluster and populate go_dist1 matrix
	n=0
	clust_all=list()
	for(i in 1:3){
		d=go_dist[[i]]
		phm=pheatmap(d,silent=TRUE)
		clust=phm$tree_row
		sl=n+1:(dim(d)[1])
		go_dist1[sl,sl] <- d[clust$order,clust$order]
		gos1=clust$labels[clust$order]
		rownames(go_dist1)[sl] <- gos1
		colnames(go_dist1)[sl] <- gos1
		n <- n+length(sl)
		clust_all[[i]]=clust
	}
	names(clust_all)=onts
	return(list(dist=go_dist1,clust=clust_all,dist_measure=dist_measure))
}

# Function to calculate similarity between two sets of GO-terms
# Input is two-column dataframe with lists of GO-terms
# Sets are subsetted by ontology and molten list with desired distances/distance matrices is returned
# Takes considerable time to process large input data (~1 minute per 1000 rows)
sim_by_ontology=function(annots,simData,measure='Wang',combine='BMA'){
	gons=lapply(simData,function(d)names(d@IC))
	t=apply(annots,1,function(i,simData,gons,measure,combine) {
		g1=unlist(i$triGO)
		g2=unlist(i$tmGO)
		lapply(simData, function(d,g1,g2,gons,measure,combine){
			#print(d@ont)
			g11=g1[which(g1 %in% gons[[d@ont]])]
			g21=g2[which(g2 %in% gons[[d@ont]])]
			#cat('g11:',g11,'g21:',g21,sep=' ',collapse='\n')
			# if(length(g11)==0 | length(g21) == 0){
			# 	return(NA)
			# 	}
			mgoSim(g11,g21,semData = d,measure=measure,combine=combine)
		},g1=g1,g2=g2,gons=gons,measure=measure,combine=combine)
	},simData=simData,gons=gons,measure=measure,combine=combine)
	require(reshape2)
	t=melt(t)
	colnames(t)=c('sim','ont','id')
	return(t[,rev(1:3)])
}

# Function will update the name with the statistic of your choice
# https://stackoverflow.com/a/47549522
#
# Create a label for the facet
# mtcars$label  <- AddNameStat(mtcars, "carb", "cyl", stat = "count")
# 
# mtcars %>% 
# 	ggplot(aes(x = cyl)) + geom_bar()+
# 	facet_wrap(~label)
AddNameStat <- function(df, category, count_col, stat = c("sd","mean","count"), dp= 0){
	require(plyr)
	# Create temporary data frame for analysis
	temp <- data.frame(ref = df[[category]], comp = df[[count_col]])
	
	# Aggregate the variables and calculate statistics
	agg_stats <- plyr::ddply(temp, "ref", summarize,
													 sd = sd(comp),
													 mean = mean(comp),
													 count = length(comp))
	
	# Dictionary used to replace stat name with correct symbol for plot
	labelName <- plyr::mapvalues(stat, 
															 from=c("sd","mean","count"), 
															 to=c("\u03C3", "x", "n"))
	
	# Updates the name based on the selected variable
	agg_stats$join <- paste0(agg_stats$ref, ": ", labelName," = ",
													 round(agg_stats[[stat]], dp))
	
	# Map the names
	name_map <- setNames(agg_stats$join, as.factor(agg_stats$ref))
	return(name_map[as.character(df[[category]])])
}

#Draw an annotated heatmap of clustered GO-terms
draw_GO_term_clusters_heatmap <- function(go_cl,go_counts,sz=10,contrast=NULL,group=FALSE,silent=FALSE,fn=FALSE){
	require(magrittr)
	require(pheatmap)
	require(viridis)
	require(pals)
	
	onts=strsplit('BP,MF,CC',',') %>% unlist
	
	ann_row=data.frame(
		Ontology=factor(rep(
			onts,
			lapply(
				go_cl$clust,
				function(l) l$labels %>% length))))
	
	if(group){
		if(is.null(contrast)){
			message('No contrast specified. Can not read group variable.')
			return()
			}
		require(xlsx)
		group=read.xlsx(
			file='Suppl_Tables.xlsx',
			sheetName=sprintf("%s_GOCluster",contrast))$group
		ann_row=cbind(ann_row,Group=group)
	}
	
	rownames(ann_row)=go_cl$dist %>% rownames()
	
	colvec.onts=viridis::plasma(3,direction = -1,begin=.15,end=0.8)
	names(colvec.onts)=onts
	ann_colors=list(Ontology=colvec.onts)
	
	if(!identical(group,FALSE)){
  	l=ann_row$Group %>% levels
  	lns=l[ann_row$Group %>% unique() %>% as.numeric()] %>% as.character()
  	#colvec=inferno(l %>% length,direction = -1,begin=.15,end=0.8)
  	colvec.groups=pals::kelly()[1:(l %>% length)] %>% rev()
  	names(colvec.groups)=lns
  	ann_colors$Group=colvec.groups
	}
	
	up=grep('Up',go_counts %>% colnames)
	up_names=(go_counts %>% colnames)[up]
	up_names=gsub('\\.',' ',up_names)
	
	ann_row[[ up_names[1] ]]=(go_counts[,up[1]] + 1) %>% log()
	ann_row[[ up_names[2] ]]=(go_counts[,up[2]] + 1) %>% log()
	ann_colors[[ up_names[1] ]]=c("white", "forestgreen")
	ann_colors[[ up_names[2] ]]=c("white", "forestgreen")
	
	if(silent){
		return(list(
			ann_row=ann_row,
			ann_colors=ann_colors))
		
	}else{
	  
		pheatmap(
			go_cl$dist,
			cluster_rows = FALSE,
			cluster_cols = FALSE,
			color=viridis::viridis(255,begin = 0.25),
			labels_row = go_counts$label,
			annotation_row = ann_row,
			annotation_colors=ann_colors,
			cellwidth=sz,
			cellheight=sz,
			fontsize=sz,
			filename = fn
		)
	  if(!identical(fn,FALSE)){
  	  pheatmap(
  	    go_cl$dist,
  	    cluster_rows = FALSE,
  	    cluster_cols = FALSE,
  	    color=viridis::viridis(255,begin = 0.25),
  	    labels_row = go_counts$label,
  	    annotation_row = ann_row,
  	    annotation_colors=ann_colors,
  	    cellwidth=sz,
  	    cellheight=sz,
  	    fontsize=sz,
  	    filename = fn
	    )
	  }else{
	    pheatmap(
	      go_cl$dist,
	      cluster_rows = FALSE,
	      cluster_cols = FALSE,
	      color=viridis::viridis(255,begin = 0.25),
	      labels_row = go_counts$label,
	      annotation_row = ann_row,
	      annotation_colors=ann_colors,
	      cellwidth=sz,
	      cellheight=sz,
	      fontsize=sz
	    )
    }
	}
}

#Draw clustering of genes by GOs and LFC, genes are colored by GO-groups assigned in heatmap
draw_GO_term_clusters_goplot <- function(
	sigGOs,
	tgos,
	res,
	term.col,
	lfc.col,
	contrast=NULL,
	ret=c('graph','data'),
	lfc=3){
	
	require(GOplot)
	require(magrittr)
	require(DESeq2)
	require(viridis)

	idx=which(sigGOs$ID %in% tgos$ID)
	tgs=sigGOs[idx,'genes'] %>%
		as.character() %>%
		strsplit(.,', ') %>%
		unlist %>%
		unique
	idx=match(tgs,res %>% rownames)
	tres=res[idx,] %>% as.data.frame()
	
	genes=data.frame(
		ID=rownames(tres) %>% toupper(),
		logFC=tres$log2FoldChange)
	
	rownames(genes)=genes$ID
	circ=circle_dat(sigGOs,genes)
	#genes=circ[circ$genes %in% tgs,c(5,6)] %>% unique()
	#Fix logFC -- some are NAs !!! there is likely a problem with circle_dat function
	circ$logFC=res$log2FoldChange[match(circ$genes,res %>% rownames %>% toupper)]
	chord=chord_dat(circ,genes,tgos %>% rownames)
	
	tst=(
		subset(
			sigGenes,
			sigGenes$tgs %in% tgs & 
				sigGenes$GO.ID %in% gos$ID
		)[,c(2,5,7)] %>% unique() %>% dim)[1] == (which(chord[,-63]==1) %>% length)
	#circ1=circ[circ$ID %in% gos$ID,]
	if(!tst){
		warning('Chord object is not consistent with input target genes for plotting! Exiting.')
		return()
	}

	if(ret=='graph'){
		lfc=c(-1,1)*lfc
		p=GOCluster2(
			data=circ,
			process=tgos$ID,
			term.col=term.col,
			lfc.col=lfc.col,
			lfc.min=lfc[1],lfc.max=lfc[2],
			chord=chord
		)
		return(p)
	}
	if(ret=='data'){
		return(list(circ=circ,chord=chord))
	}
}

#Modified GOCluster function from GOplot
#Accepts specific 'chord' object instead of looking it in global environment.
GOCluster2 <- function (data, chord, process, metric, clust, clust.by, nlfc, lfc.col, 
					lfc.min, lfc.max, lfc.space, lfc.width, term.col, term.space, 
					term.width) 
{
	x <- y <- xend <- yend <- width <- space <- logFC <- NULL
	if (missing(metric)) 
		metric <- "euclidean"
	if (missing(clust)) 
		clust <- "average"
	if (missing(clust.by)) 
		clust.by <- "term"
	if (missing(nlfc)) 
		nlfc <- 0
	if (missing(lfc.col)) 
		lfc.col <- c("firebrick1", "white", "dodgerblue")
	if (missing(lfc.min)) 
		lfc.min <- -3
	if (missing(lfc.max)) 
		lfc.max <- 3
	if (missing(lfc.space)) 
		lfc.space <- (-0.5)
	else lfc.space <- lfc.space * (-1)
	if (missing(lfc.width)) 
		lfc.width <- (-1.6)
	else lfc.width <- lfc.space - lfc.width - 0.1
	if (missing(term.col)) 
		term.col <- brewer.pal(length(process), "Set3")
	if (missing(term.space)) 
		term.space <- lfc.space + lfc.width
	else term.space <- term.space * (-1) + lfc.width
	if (missing(term.width)) 
		term.width <- 2 * lfc.width + term.space
	else term.width <- term.width * (-1) + term.space
	if (clust.by == "logFC") 
		distance <- stats::dist(chord[, dim(chord)[2]], method = metric)
	if (clust.by == "term") 
		distance <- stats::dist(chord, method = metric)
	cluster <- stats::hclust(distance, method = clust)
	dendr <- dendro_data(cluster)
	y_range <- range(dendr$segments$y)
	x_pos <- data.frame(x = dendr$label$x, label = as.character(dendr$label$label))
	chord <- as.data.frame(chord)
	chord$label <- as.character(rownames(chord))
	all <- merge(x_pos, chord, by = "label")
	all$label <- as.character(all$label)
	if (nlfc) {
		lfc_rect <- all[, c(2, dim(all)[2])]
		for (l in 4:dim(data)[2]) lfc_rect <- cbind(lfc_rect, 
																								sapply(all$label, function(x) data[match(x, data$genes), 
																																									 l]))
		num <- dim(data)[2] - 1
		tmp <- seq(lfc.space, lfc.width, length = num)
		lfc <- data.frame(x = numeric(), width = numeric(), 
											space = numeric(), logFC = numeric())
		for (l in 1:(length(tmp) - 1)) {
			tmp_df <- data.frame(x = lfc_rect[, 1], width = tmp[l + 
																														1], space = tmp[l], logFC = lfc_rect[, l + 1])
			lfc <- rbind(lfc, tmp_df)
		}
	}
	else {
		lfc <- all[, c(2, dim(all)[2])]
		lfc$space <- lfc.space
		lfc$width <- lfc.width
	}
	term <- all[, c(2:(length(process) + 2))]
	color <- NULL
	termx <- NULL
	tspace <- NULL
	twidth <- NULL
	for (row in 1:dim(term)[1]) {
		idx <- which(term[row, -1] != 0)
		if (length(idx) != 0) {
			termx <- c(termx, rep(term[row, 1], length(idx)))
			color <- c(color, term.col[idx])
			tmp <- seq(term.space, term.width, length = length(idx) + 
								 	1)
			tspace <- c(tspace, tmp[1:(length(tmp) - 1)])
			twidth <- c(twidth, tmp[2:length(tmp)])
		}
	}
	tmp <- sapply(lfc$logFC, function(x) ifelse(x > lfc.max, 
																							lfc.max, x))
	logFC <- sapply(tmp, function(x) ifelse(x < lfc.min, lfc.min, 
																					x))
	lfc$logFC <- logFC
	term_rect <- data.frame(x = termx, width = twidth, space = tspace, 
													col = color)
	legend <- data.frame(x = 1:length(process), label = process)
	
	#https://gist.github.com/dgrtwo/38e70de658b48e166f90
	theme_blank <- function(...) {
		ret <- theme_bw(...)
		ret$line <- element_blank()
		ret$rect <- element_blank()
		ret$strip.text <- element_blank()
		ret$axis.text <- element_blank()
		ret$plot.title <- element_blank()
		ret$axis.title <- element_blank()
		ret$plot.margin <- structure(c(0, 0, -1, -1), unit = "lines", valid.unit = 3L, class = "unit")
		ret
	}
	
	ggplot() + 
		geom_segment(data = segment(dendr), aes(x = x, y = y, xend = xend, yend = yend)) + 
		geom_rect(data = lfc, aes(xmin = x - 0.5, xmax = x + 0.5, ymin = width, ymax = space, fill = logFC)) + 
		scale_fill_gradient2(
			"logFC", space = "Lab", low = lfc.col[3], mid = lfc.col[2], high = lfc.col[1], 
			guide = guide_colorbar(title.position = "top", title.hjust = 0.5),
			breaks = c(min(lfc$logFC), max(lfc$logFC)),
			labels = c(round(min(lfc$logFC)), round(max(lfc$logFC)))) + 
		geom_rect(
			data = term_rect,
			aes(xmin = x - 0.5, xmax = x + 0.5, ymin = width, ymax = space),
			fill = term_rect$col) +
		geom_point(
			data = legend,
			aes(x = x, y = 0.1, size = factor(label, levels = label),	shape = NA)) + 
		guides(size = guide_legend(
			"GO Terms", ncol = 4, byrow = T, 
			override.aes = list(shape = 22, fill = term.col, size = 8))) + 
		coord_polar() + 
		scale_y_reverse() + 
		theme(
			legend.position = "bottom", 
			legend.background = element_rect(fill = "transparent"),
			legend.box = "horizontal",
			legend.direction = "horizontal") + 
		theme_blank()
}


draw_COG_histogram <- function(triR, ids){
	#TriR -- trinotateR object from which to get mapping of eggNOG genes
	#ids  -- named list of sets of transcript ids for which to plot histogram
	require(trinotateR)
	require(ggplot2)
	require(viridis)
	require(reshape2)
	require(magrittr)
	
	#get egg and cog data
	ann_dir='~/matrosov/axSA_assemblies/drap/stranded/trinity_norm_altogether/trinotate/'
	egg=read.table(paste0(ann_dir,"NOG.annotations.tsv"), sep="\t", stringsAsFactors=FALSE, quote="")
	names(egg) <- c("db", "nog", "proteins", "species", "class", "description")
	data(cogs)
	#plot_NOGs(ta_clean)
	
	#prepare dataframe to plot histogram
	df=lapply(ids,function(l) {
		y=split_eggnog(ta_clean[l,])
		nogs <- gsub("(.*)\\^.*", "\\1", y$egg)
		n <- match(nogs, egg$nog)
		y <- table( unlist( strsplit(egg$class[n], "")) )
		#rev( y[ order( match(names(y) , cogs$code) ) ] )
		t=rep(0,cogs$code %>% length)
		names(t)=cogs$code
		t[match(names(y),names(t))]=y
		return(as.list(t))
	}) %>% melt
	names(df)=strsplit('count,cog,category',',') %>% unlist
	
	#remove empty categories
	del_codes=cogs$code[
		sapply(
			cogs$code,
			function(i) sum(df[df$cog==i,1])) == 0]
	df=df[!df$cog %in% del_codes,]
	
	#sort categories as these appear in ids list
	df$category=factor(df$category, levels=names(ids))
	
	#make labels and calculate percentages
	df$label=paste0(
		paste0(df$cog,'. '),
		cogs$name[match(df$cog,cogs$code)])
	df$freq=apply(
		df, 1, 
		function(r,df) 
			sum(df[df$category==r[3],'count']),df=df )
	df$freq=sprintf('%.2f',df$count/df$freq*100) %>% as.numeric()
	
	#make plot
	p=ggplot(data=df,aes(x=label,y=freq,fill=category)) +
		geom_bar(stat='identity',position=position_dodge()) +
		scale_fill_viridis(
			name='Set of transcripts',
			discrete = TRUE,
			begin=0.02,end=0.98,alpha=0.8) +
		theme_minimal() +
		#coord_cartesian(ylim=c(0,15)) +
		theme(
			axis.text.x = element_text(angle = 45, hjust = 1),
			legend.position =  'left') + 
		ylab('percent') + 
		xlab('') +
		scale_y_sqrt(breaks=c(1,2,5,10,20,30))
	return(p)
}

get_slim_GOs <- function(gos){
	require(GSEABase)
	
	onts=strsplit('BP,MF,CC',',') %>% unlist
	fl <- "http://current.geneontology.org/ontology/subsets/goslim_generic.obo"
	oboSlim <- getOBOCollection(fl)
	slimres=lapply(onts, function(x) goSlim(GOCollection(gos), oboSlim, x))
	names(slimres)=onts
	return(slimres)
}

get_slim_topGOs <-function(slimres,min=5){
	require(GO.db)
	require(magrittr)
	
	#onts=strsplit('BP,MF,CC',',') %>% unlist
	onts=names(slimres)
	GO2del=strsplit('biological_process,molecular_function,cellular_component,cell',',') %>% unlist
	slimres=lapply(names(slimres), function(x) {
		y=slimres[[x]]
		y=y[y$Percent>=min,]
		y$Term=sapply(rownames(y), function(z) Term(GOTERM[z]))
		y=y[!y$Term %in% GO2del, ]
		y[order(y$Count,decreasing = TRUE),]
	})
	names(slimres)=onts
	return(slimres)
}

###SNP functions
###These functions are related with analysis SNPs obtained by KISS analysis
#Get contingency table for distribution of SNPs across coding/non-coding etc
snp_get_good_snps <- function(f,pval=0.05,delta=0.1){
	library(magrittr)
	
	snp=read.csv(file=f,sep="\t",header = TRUE)
	snp$X.Component_ID %<>% as.character
	snp$KissDE_p.value %<>% as.character %>% as.double()
	snp$KissDE_DeltaF %<>% as.character %>% as.double()
	snp$Is_in_CDS %<>% as.logical(.)
	snp$Is_not_synonymous %<>% as.logical(.)
	snp$Is_condition_specific %<>% as.logical(.)
	snp$SNP_in_mutliple_assembled_genes %<>% as.logical(.)
	snp$SNP_in_mutliple_assembled_isoforms %<>% as.logical(.)
	snp$Possible_sequencing_error %<>% as.logical(.)
	snp=snp[which(snp$KissDE_p.value<=pval & snp$KissDE_DeltaF>delta),]
	
	return(snp)
}

snp_get_distr <- function(annot,res,snps,sval=0.05,tot=FALSE,counts=TRUE){
	library(magrittr)
	
	TRI_nc=annot$gene_id[is.na(annot$prot_id)]
	TRI_c=annot$gene_id[!is.na(annot$prot_id)]
	DE=rownames(res)[which(res$svalue<=sval)]
	SNP_c=snps$X.Component_ID[snps$Is_in_CDS]
	SNP_nc=snps$X.Component_ID[!snps$Is_in_CDS]
	SNP_c_nsyn=snps$X.Component_ID[snps$Is_in_CDS & snps$Is_not_synonymous]
	
	#Table 5x2
	#Table format is:
	#								1								2
	#						NDE-transcripts/DE-transcripts
	#	1 C-S
	#	2 NC-S
	#	3 NC+S
	#	4 C+S+Syn
	#	5 C+S+Nsyn
	d=list(
	##### Row 1
	which(!TRI_c %in% SNP_c & !TRI_c %in% DE),
	#C-S/NDE 21556
	which(!TRI_c %in% SNP_c & TRI_c %in% DE),
	#C-S/DE 296
	
	##### Row 2
	which(!TRI_nc %in% SNP_nc & !TRI_nc %in% DE),
	#NC-S/NDE 5414
	which(!TRI_nc %in% SNP_nc & TRI_nc %in% DE),
	
	##### Row 3
	which(TRI_nc %in% SNP_nc & !TRI_nc %in% DE),
	#NC+S/NDE 12
	which(TRI_nc %in% SNP_nc & TRI_nc %in% DE),
	#NC-S/DE 0
	
	#NC-S/DE 24
	# ##### Row 4 all SNPs in CDS
	# which(TRI_c %in% SNP_c & !TRI_c %in% DE)
	# #C+S/NDE 138
	# which(TRI_c %in% SNP_c & TRI_c %in% DE)
	# #C+S/DE 6
	
	##### Row 4
	which(TRI_c %in% SNP_c & !TRI_c %in% SNP_c_nsyn & !TRI_c %in% DE),
	#C+S+Syn/NDE 84
	which(TRI_c %in% SNP_c & !TRI_c %in% SNP_c_nsyn & TRI_c %in% DE),
	#C+S+Syn/DE 1
	
	##### Row 5
	which(TRI_c %in% SNP_c & TRI_c %in% SNP_c_nsyn & !TRI_c %in% DE),
	#C+S+Nsyn/NDE 54
	which(TRI_c %in% SNP_c & TRI_c %in% SNP_c_nsyn & TRI_c %in% DE)
	
	#C+S+NSyn/DE 5
	)
	if(!counts){
		d=list(
			'C-S'     =list('NDE'=TRI_c[ d[[1]] ],  'DE'=TRI_c[ d[[2]] ]),
			'NC-S'    =list('NDE'=TRI_nc[ d[[3]] ], 'DE'=TRI_nc[ d[[4]] ]),
			'NC+S'    =list('NDE'=TRI_nc[ d[[5]] ], 'DE'=TRI_nc[ d[[6]] ]),
			'C+S+Syn' =list('NDE'=TRI_c[ d[[7]] ],  'DE'=TRI_c[ d[[8]] ]),
			'C+S+Nsyn'=list('NDE'=TRI_c[ d[[9]] ],  'DE'=TRI_c[ d[[10]] ])
		)

		return(d)
	}
	
	d=lapply(d,function(x) lapply(x,function(y) length(y))) %>% unlist %>% unname
	m=matrix(d,byrow = TRUE,ncol=2)
	m=cbind(m,rowSums(m))
	m=rbind(m,colSums(m))
	rownames(m)=strsplit('C-S,NC-S,NC+S,C+S+Syn,C+S+Nsyn,Tot',',') %>% unlist
	colnames(m)=c('NDE','DE','Tot')
	if(tot){return(m)}
	else{return(m[-6,-3])}
}