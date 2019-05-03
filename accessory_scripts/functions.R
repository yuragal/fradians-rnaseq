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