require(magrittr)

check_uniprot_annot <-function(ta_clean,idx,ta_kos){

	# This is to add statistics on blast hits with query coverage per subject
	# Initial Trinotate NR_BBp field do not have this stat, though it is valuable for downstream manual annotation
	#   !!!Blast.out file contains by one query hit per subject. Multiple hits were filtered.
	nrBBpHs=read.csv("~/matrosov/axSA_assemblies/drap/stranded/trinity_norm_altogether/trinotate/trinotate.nr.blastp.qcovstat.out",header=FALSE, sep = "\t",stringsAsFactors = FALSE)
	colnames(nrBBpHs)=strsplit("qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp"," ") %>% unlist
	nrBBpHs$saccver %<>% sub('\\.\\d','',.)
	nrBBpHs$nr_blastp=apply(
		nrBBpHs[,c(2,3,11,13,14)],1,
		function(i) paste(i %>% sub(' ','',.),collapse='|'))	
	nrids=ta_clean$nr_BLASTP[idx] %>%
		lapply(.,function(i) strsplit(i,'\\^')) %>%
		unlist(.,recursive = FALSE) %>%
		lapply(.,function(i) i[1]) %>%
		sapply(.,function(i) i[1]) %>%
		unname %>% sub('\\.\\d$','',.)
	names(nrids)=ta_clean$gene_id[idx]
	
	# # Query NCBI for GB-files
	# library(rentrez)
	# b=nrids %>% unique
	# chunk=100
	# sapply(seq(1,length(b),by=chunk), function(i){
	# 	cat(i,'\n')
	# 	ef=entrez_fetch(
	# 		db='protein',
	# 		rettype="GB",
	# 		id=b[i:(i+chunk-1)]
	# 	)
	# 	write(ef,'tmp.gb',append=TRUE)
	# })
	
	#rm(list=c('a','b','idx','dia_kos','ta_kos'))

	# el=entrez_link(dbfrom='protein',id = es$ids,db='taxonomy')
	# 
	# ef=entrez_fetch(
	# 	db='protein',
	# 	rettype="GB",
	# 	#use_history = TRUE,
	# 	#retmax=100,
	# 	id=(nrids %>% unique)[1]
	# )
	
	#!!!run the following bash command from teminal:
	#cat tmp.gb | awk '/^ACCESSION/{printf("%s\t%d\t%s\t%s\n",a[1],a[3],a[4],a[2]);delete a;a[1]=$2;;next} /^SOURCE/{a[2]=$2" "$3;next} /\/db_xref=\"taxon:/{split($0,b,":");a[3]=substr(b[2],1,length(b[2])-1);next} /\/product=/{s=$0;while(s~/[^"]$/){getline;s=s $0};gsub(/\s+/," ",s);gsub(/- -/,"-",s);split(s,b,"\"");a[4]=b[2];next}' | sed '1d' > tmp.tx
	
	#Let's explore NR BBpHs by taxonomy. Here are hits' top-20 species:
	# 1     723 Phaeodactylum tricornutum					diatom
	# 2     306 Fragilariopsis cylindrus					diatom
	# 3     243 Thalassiosira oceanica						diatom
	# 4     205 Thalassiosira pseudonana					diatom
	# 5      17 synthetic construct							!!mostly human-related sequences; highly likely contaminant
	# 6      12 Guillardia theta									cryptophyte
	# 7      11 Ectocarpus siliculosus						filamentous brown alga
	# 8      10 Nannochloropsis gaditana					ochrophyte
	# 9       8 Chrysochromulina sp.							haptophyte
	# 10       8 Phytophthora parasitica					oomycete
	# 11       7 Emiliania huxleyi								haptophyte
	# 12       6 Aphanomyces invadans							oomycete
	# 13       6 Escherichia coli								!!bacteria -- highly likely HGT or contaminant
	# 14       6 Leptonychotes weddellii				!!Weddel seal, antarctic species; highly likely HGT (from baikal seal) or contaminant
	# 15       5 Aureococcus anophagefferens			pelagophytes, relatives of diatoms
	# 16       5 Saprolegnia parasitica						oomycete
	# 17       4 Aphanomyces astaci								oomycete
	# 18       4 Cricetulus griseus							!!chinese hamster -- highly likely contaminant
	# 19       4 Saprolegnia diclina							oomycete
	# 20       4 Strongylocentrotus purpuratus  !!sea urchin, echinodermata -- highly likely HGT or contaminant
	
	#These are species with acceptable taxonomy:
	taxok=strsplit('Phaeodactylum tricornutum|Fragilariopsis cylindrus|Thalassiosira oceanica|Thalassiosira pseudonana|Guillardia theta|Ectocarpus siliculosus|Nannochloropsis gaditana|Chrysochromulina sp.|Phytophthora parasitica|Emiliania huxleyi|Aphanomyces invadans|Aureococcus anophagefferens|Saprolegnia parasitica|Aphanomyces astaci|Saprolegnia diclina','\\|') %>% unlist
	
	#transfer annotations from Trinotate for items having reasonable NR BBp hit
	#with taxonomy in taxok
	#I will transfer annotation from nrids2tx$annot unless these are "predicted proteins"
	# if annotation in nr2txids$annot have no sense, I will check the gene name in Uniprot/SwissProt Trinotate hit
	# accs with hits in taxok
	
	nrids2tx=read.csv('tmp.tx',header=FALSE,sep='\t',stringsAsFactors = FALSE)
	names(nrids2tx)=c('acc','txid','annot','species')
	rownames(nrids2tx)=nrids2tx$acc
	nrids2tx$annot %<>% sub('(hypothetical|predicted|expressed|uncharacteri[zs]ed) protein.*',NA,ignore.case=TRUE,.) %>%
		sub('UNKNOWN','',ignore.case=TRUE,.)
	#nrids2tx[idx,] %>% View
	a=which(nrids2tx$species %in% taxok) %>% nrids2tx$acc[.]
	# length(nrids2tx$acc)
	# length(a)
	a=nrids %in% a %>% which %>% names(nrids)[.]
	# length(nrids)
	# length(a)
	# 1557 out of 1968 nr hits are in taxok
	# transcripts mapped to nrids: 2426 out of 2962
	#idx=idx[match(nrids,idx) %>% na.omit()]
	df=ta_clean[ta_clean$gene_id %in% a, c(5,7,9,14)]
	df$ko=ta_kos[a]
	df=cbind(df,nrids2tx[match(nrids[a],nrids2tx$acc),])
	#df=cbind(df[,c(1:5)],nrids2tx[match(df$V3, nrids2tx$acc),])
	df=cbind(df,df$sprot_Top_BLASTP_hit %>%
					 	sapply(.,function(i) {
					 		l=strsplit(i,'\\^') %>%
					 			unlist(.,recursive = FALSE)
					 		d=l[6]
					 		d %<>% sub('RecName: Full=','',.) %>%
					 			sub(';$','',.) %>%
					 			strsplit(.,' \\{ECO:') %>%
					 			unlist()
					 		if(length(d)==1){
					 			return(c(d[1],NA))
					 		}
					 		else{
					 			return(c(d[1],paste0('ECO:',d[2] %>% sub('}','',.))))
					 		}
					 	}) %>% 
					 	unname %>%
					 	t
	)
	colnames(df)[c(5,10,11)]=c('ko','annot_uniprot','ecode')
	
	#Now merge ann and df
	# add column source:
	#   descriptions which are obtained from diatom blast-hits with KEGG-annotations are marked by D
	#   those from uniprot but with NR BBp hit taxonomy assigned to protists are marked by U:protists
	#     NOTE: all 'U:protists'-marked transcripts have kegg orthologs,
	#     	but there are D-marked transcripts with no kegg orthologs
	#     	(only definition from kegg entry is available)
	
	df$Kegg %<>% gsub('(KEGG|KO):','',.)
	df$Kegg_hit=sapply(df$Kegg, function(i)
		strsplit(i,'`') %>% unlist %>% .[1])
	df$nr_blastp=apply(df,1,function(i){
		paste(c(
			nrBBpHs$nr_blastp[which(nrBBpHs$qaccver==i[1] & nrBBpHs$saccver==i[6])],
			i[9]),
			collapse='|')
	})
	df=df[,c(1,8,10,5,5,12,2,13)]
	df$src='u:u:p'
	df[,5]=NA
	df=df[,c(1,9,2:8)]
	names(df)=c("prot_id", "src", "annot", "annot_uniprot", "ko", "ko_desc", "Kegg_hit", "sprot_blastp", "nr_blastp")
	
	# It takes ~30' to query KEGG
	# Seems, such queries are better to execute in parallel mode
	# with 5 threads it takes ~6 minutes!
	# library(parallel)
	# cl=makeCluster(5); date()
	# kos=df$ko %>% unique()
	# clusterExport(cl,'kos')
	# kos=parSapply(cl,seq(1,length(kos), by=10), function(i) {
	# 	#cat(sprintf("KeggGet, %25s: %02.0f%% done\n",date(), i/(kos %>% length)*100))
	# 	#cat(start, 'keggGet', date(),i,"\n")
	# 	c(KEGGREST::keggGet(kos[i:(i+9)]))
	# })
	# date(); stopCluster(cl)
	# 
	#To load kos, run load('trinotate_merge_KEGG.RData')
	
	attach('Trinotate_merge_KEGG.RData')
	kos_d=sapply(kos %>% unlist(recursive = FALSE),function(i) i$DEFINITION)
	names(kos_d)=sapply(kos %>% unlist(recursive = FALSE),function(i) i$ENTRY)
	detach()
	
	df$ec=''
	df=df[,c(1:6,10,7:9)]
	df$ko_desc=kos_d[df$ko]
	df$ec=sub("^.+\\[EC:","EC:",df$ko_desc) %>% gsub("\\]","",.)
	df$ec[!grepl("\\[EC:",df$ko_desc)]=NA
	df$ko_desc %<>% sub(" \\[EC:.+","",.)
	return(df)
}

check_diatom_annots <-function(ta_clean){
	#Prepare queries for KEGG
	BBpHs=read.csv("~/matrosov/axSA_assemblies/drap/stranded/trinity_norm_altogether/trinotate/TD2SAgenome/vs_diatoms.ncbi/three_diatoms.BBpHits.out",header=FALSE, sep = "\t",stringsAsFactors = FALSE)
	colnames(BBpHs)=strsplit("qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp"," ") %>% unlist
	saccs=BBpHs$saccver %>% unique %>% sort %>% sub('\\.\\d$','',.)
	#save(saccs,file="kegg_accs.RData")
	
	#15920 genes to be queried
	#Quering KEGG takes quite a long time...
	
	# Run Rscript in background... It will take ~3-4 hours to complete
	# Alternatively, one cat run cmd in rstudio terminal window
	#shell(cmd="/home/yuragal/tools/dev/R-3.5.0/bin/Rscript kegg_genes.R", wait=FALSE)
	
	#This is to run all quieries as 5 parallel processes -- it takes ~ 1 hour for 5 threads
	#!!!Run it in terminal session from within rstudio
	#By some unknown reason related with curl it is failing to work from shell session outside rstudio
	#export offset=3200 by=100; seq 1 $offset 15920 | xargs -n 1 bash -c 'end=`echo $0 + $offset - 1 | bc`; echo -e "$0 $end $by"' | awk '{print NR%5*2 " " $0}' | xargs -n4 -P 5 bash -c 'sleep $0; /home/yuragal/tools/dev/R-3.5.0/bin/Rscript kegg_genes.R $1 $2 $3'
	#export offset=1600 by=100 nofitems=15920; seq 1 $offset $nofitems | xargs -n 1 bash -c 'end=`echo $0 + $offset - 1 | bc`; echo -e "$0 $end $by"' | awk -vOS=$offset -vNOI=$nofitems '{print NR%sprintf("%d",NOI/OS+0.5) " " $0}' | xargs -n4 -P10 bash -c 'sleep $0; /home/yuragal/tools/dev/R-3.5.0/bin/Rscript kegg_genes.R $1 $2 $3'
	# Then load data into current environment:
	
	# tkgids=list()
	# tkgs=list()
	# for(i in sprintf("kegg_genes.%d.RData",seq(1,15920,by=3200))){
	# 	load(i)
	# 	tkgids=c(tkgids,kgids)
	# 	tkgs=c(tkgs,kgs)
	# }
	# kgids=tkgids %>% unlist
	# kgs=tkgs %>% unlist(recursive = FALSE)
	# names(kgs) %<>% sub("ncbi-proteinid:","",.)
	# names(kgids) %<>% sub("ncbi-proteinid:","",.)
	# save(kgids,kgs,file="kegg_genes.RData")
	# rm(list=c('tkgids','tkgs'))
	
	attach("kegg_genes.RData")
	defs=do.call(
		rbind.data.frame,
		lapply(kgs, function(i)
			if(is.null(i$ORTHOLOGY)){
				list(i$ENTRY,i$DEFINITION,NA,NA)
			}else{
				list(i$ENTRY,i$DEFINITION,names(i$ORTHOLOGY),i$ORTHOLOGY)
			}
		)
	)
	detach()
	
	colnames(defs)=c('entry','definition','ko','orth_desc')
	defs=defs[
		!((grepl('hypothetical protein',defs$definition) |
			 	grepl('predicted protein',defs$definition)) &
				is.na(defs$ko)),]
	defs$entry %<>% sub(" CDS","",.)
	defs$definition %<>% sub("\\(.*\\) ","",.)
	defs$ec=sub("^.+\\[EC:","EC:",defs$orth_desc) %>% gsub("\\]","",.)
	defs$ec[!grepl("\\[EC:",defs$orth_desc)]=NA
	defs$orth_desc %<>% sub(" \\[EC:.+","",.)
	
	idx=gsub("\\.\\d","",BBpHs$saccver) %>% match(.,rownames(defs)) %>% na.omit()
	kegg_annots=cbind(
		BBpHs[
			setdiff(
				1:dim(BBpHs)[1],
				attributes(idx)$na.action),
			c(1:3,11,13,14)],
		defs[idx,])
	kegg_annots=kegg_annots[kegg_annots$qaccver %in% ta_clean$prot_id,]
	
	#merge multiple diatom hits
	ann=do.call(
		rbind.data.frame,
		lapply(kegg_annots[,1] %>% unique, function(i) {
			hits=kegg_annots[which(kegg_annots$qaccver == i),]
			
			#order descriptions by length
			
			gene_desc=combine_gene_descriptions(descr = hits$definition)
			
			c(
				i,
				gene_desc,
				hits[1,9:11],
				hits[,2:7] %>% apply(.,1,function(i) paste(i,collapse='|')) %>% paste(.,collapse="^")
			)
		})
	)
	names(ann)=strsplit('prot_id description ko orth_desc ec blast_stat',' ') %>% unlist
	ann$prot_id %<>% as.character
	ann$blast_stat %<>% as.character
	ann$description %<>% as.character
	ann$ko %<>% as.character
	ann$orth_desc %<>% as.character
	ann$ec %<>% as.character
	return(ann)
}

diatom2uniprot_annots <- function(ta_clean,dann,uann){
	require(magrittr)
	
	dann$sprot_stat=ta_clean$sprot_Top_BLASTP_hit[match(dann$prot_id, ta_clean$prot_id)]
	dann$src='n:n:d'
	dann$Kegg_hit=sapply(dann$blast_stat, function(i){
		h=strsplit(i,'\\^') %>%	unlist
		h=strsplit(h[1],'\\|') %>% unlist
		h=h[length(h)]
		if(	grepl('FRACY',h)){
			o='fcy'
		}
		else if(grepl('THAPS',h)){
			o='tps'
		}else{
			o='pti'
		}
		paste(c(o,h),collapse=':')
	}
	) %>% unname
	dann$annot_uniprot=sapply(dann$sprot_stat,function(i) {
		l=strsplit(i,'\\^') %>%
			unlist(.,recursive = FALSE)
		d=l[6]
		d %<>% sub('RecName: Full=','',.) %>%
			sub(';$','',.) %>%
			strsplit(.,' \\{ECO:') %>%
			unlist()
		d[1]
	}) %>% unname
	
	colnames(dann)[2]='annot'
	colnames(dann)[4]='ko_desc'
	colnames(dann)[6]='nr_blastp'
	colnames(dann)[7]='sprot_blastp'
	
	dann=dann[,match(colnames(uann),colnames(dann))]

	# !!! Turned out 916 transcripts were in both Diatom and SwissProt sets,
	# 			but have no KO there in dann.
	# In UniProt/SwissProt set annotations seem more complete for duplicated transcripts
	
	idx1=uann$prot_id %in% dann$prot_id
	idx2=dann$prot_id %in% uann$prot_id
	#These are duplicates
	#rbind(uann[idx1,],dann[idx2,]) %>% .[order(.$prot_id),] %>% View
	#rbind(uann[idx1,],dann[!idx2,]) %>% View
	
	# I will merge such transcripts setting annot
	# 
	# 
	# 
	z=uann
	apply(dann[which(idx2),], 1, function(d){
		idx=which(z[,'prot_id'] %in% d['prot_id'])
		#cat('D:',d[1:5],collapse='\n')
		src=deparse_source(z[idx,'src'])
		if(is.na(z[idx,'annot']) |
			 z[idx,'annot'] %>% tolower() == d['annot'] %>% tolower ){
			
			src=edit_source(src,action = 'replace',field='annot',value='n',mode='vec')
			src=edit_source(src,action = 'replace',field='tax',value='d')
			z[idx,'annot']<<-d['annot']
		}else{
			src=edit_source(src,action = 'add',field='annot',value='n',mode='vec')
			src=edit_source(src,action = 'replace',field='tax','d')
			z[idx,'annot']<<-paste(
				paste0('D:',d['annot']),
				paste0('U:',z[idx,'annot']),
				sep='^')
		}
		z[idx,'src']<<-src
		z[idx,'nr_blastp']<<-d['nr_blastp']
		return()
	})
	z=rbind(z,dann[!idx2,])
	
	#  Another 80 diatom-annotated transcripts have no KO assigned
	#  but KO are present for them in ta_clean
	#  Add these kegg annotations to transcripts, edit the src field
	ta_kos=ta_clean$prot_id[which(!is.na(ta_clean$prot_id) & grepl('KO:',ta_clean$Kegg))]
	idx=(z$prot_id %in% ta_kos) & is.na(z$ko)
	idx1=match(z$prot_id[idx],ta_clean$prot_id)
	
	
	koids=sapply(
		ta_clean$Kegg[which(ta_clean$prot_id %in% z$prot_id[idx])],
		function(i)  strsplit(i,'KO:') %>% unlist %>% .[2]
		)
	kos=get_kegg_data(koids)
	
	kos_d=sapply(kos,function(i) i$DEFINITION)
	names(kos_d)=sapply(kos,function(i) i$ENTRY)
	
	z$ko[idx]=names(kos_d)
	z$ko_desc[idx]=sub(" \\[EC:.+","",kos_d)
	z$ec[idx]=sub("^.+\\[EC:","EC:",kos_d) %>% gsub("\\]","",.)
	z$ec[!grepl('EC:',z$ec)]=NA
	z$src[idx] %<>% sapply(.,function(i) edit_source(i,action = 'replace',field='kegg',value='u'))
	
	return(z)
}

get_kegg_data <- function(ids){
	library(parallel)
	a=ids
	cl=makeCluster(5); date()
	clusterExport(cl,'a')
	kos=parSapply(cl,seq(1,length(a), by=10), function(i) {
		#cat(sprintf("KeggGet, %25s: %02.0f%% done\n",date(), i/(kos %>% length)*100))
		#cat(start, 'keggGet', date(),i,"\n")
		c(KEGGREST::keggGet(a[i:(i+9)]))
	})
	date(); stopCluster(cl)
	return(kos)
}

# These functions to work with src filed of dataframes
# returned by diatom2uniprot_annots and check_uniprot_annot functions
# Currently, the source (src) filed format is 
# Annotation : Kegg : Taxonomy
# name      annot      kegg   tax
# value [nu]{1,2}:[nu]{1,2}:[dpo]
# 
# Where annot -- annotation source (UniProt/SwissProt or/and NCBI GenBank)
#       kegg  -- how the kegg info was mapped (sources same as above)
#       tax   -- taxonomy of the blastp-hit (Diatoms/Protists/Other)
deparse_source <- function(src){
	src=strsplit(src,':') %>% unlist
	names(src)=strsplit('annot,kegg,tax',',') %>% unlist
	return(src)
}

edit_source <- function(
	src,
	action=c('add','clear','replace'),
	field=c('annot','kegg','tax'),
	value,
	mode='str'){

	src=deparse_source(src)

	if(action=='add'){
		if(field=='tax'){
			warning('Only clear and replace actions allowed for Taxonomy field!!!')
		} else{
			src[field]=paste0(src[field],value)
		}
	}
	else if(action=='clear'){
		src[field]=''
	}
	else{
		src[field]=value
	}
	if(mode=='str'){
		src=paste(src,collapse=':')
	}
	return(src)
}

#work with annotations
combine_gene_descriptions <-function(descr){
	ds=descr[!grepl("hypothetical|predicted|expressed protein",descr)] %>% 
		tolower() %>%
		unique
	#ds -- temp descriptions
	#o  -- indices of descriptions ordered by string lengths
	#d  -- distance matrix between strings (various metrics used)
	#      probably, algorithm and metrics used should be revised
	#        for more general-purpose applications
	o=order(nchar(ds),ds) %>% rev
	#use Jaro-Winker distance first
	#jw seems more robust for general fuzzy matching purposes
	d=stringdist::stringdistmatrix(ds[o],method='jw')
	if(length(ds)==0){
		gd='hypothetical protein'
	}
	else if(length(ds)==1){
		gd=descr[match(ds[o[1]],descr %>% tolower)[1]]
	}
	else if(min(d)<=0.5){
		# the first longest description is considered to be the best
		# as it 0.5 similar with at least one another description
		
		gd=descr[match(ds[o[1]],descr %>% tolower)[1]]
		#print(gd)
	}
	else{
		#calculate Levenshtein distance matrix
		#It seems more stable when the shorter query string is a substring of much longer subject string
		d=adist(ds[o],partial = TRUE)
		mind=min(d[-1,1])/nchar(ds[o[1]])
		
		if(mind<=0.35){
			# the first longest description is considered to be the best
			# as it is a substring of least one another description
			gd=descr[match(ds[o[1]],descr %>% tolower)[1]]
			#print(gd)
		}
		else{
			cat(sprintf(
				"Similar to %0.2f; This is lower than 0.35: %s.",
				mind,
				mind<=0.35))
			
			#Remove by ambiduous abbreviations in description
			if(grepl("hsp|nad",ds[o[1]]) & (grep("hsp|nad",ds[o[2:length(o)]]) %>% length > 0)){
				cat(" Remove ambiguities.")
				ds=ds[c(1,!grepl("hsp|nad",ds[o[2:length(o)]]))]	
			}
			else{
				ds=ds[o]
			}
			
			warning(sprintf("Can not find good unique description for transcript %s:\n",i))
			#print(t[o])
			#print(d)
			gd=paste(descr[match(ds,descr %>% tolower)],collapse='|')
			cat(' Final gene description:',gd,'\n')
		}
		return(gd)
	}
}

#work with annotations
parse_sprot_blasthit <- function(h){
	l=strsplit(h,'\\^') %>%
		unlist(.,recursive = FALSE)
	d=l[6]
	d %<>% sub('RecName: Full=','',.) %>%
		sub(';$','',.) %>%
		strsplit(.,' \\{ECO:') %>%
	unlist()
	return(d[1])
}

#work with annotations
convert_trinotate_records <-function(ta_clean,idx,nrBBpHs){
	
	df=ta_clean[idx,c(5,7,14)]
	names(df)[2]='sprot_blastp'
	df$src='u:u:o'
	df$annot=NA
	df$annot_uniprot=sapply(df$sprot_blastp,parse_sprot_blasthit)
	
	df$ko=sapply(
		df$Kegg,
		function(i)  strsplit(i,'KO:') %>%
			unlist %>%
			.[2]
	)
	df$Kegg_hit=sapply(
		df$Kegg,
		function(i)  strsplit(i,'KO:') %>%
			unlist %>%
			.[1] %>%
			sub('KEGG:','',.) %>%
			sub('`.*$','',.)
	)
	kos=get_kegg_data(df$ko)
	kos_d=sapply(kos %>% unlist(recursive = FALSE), function(i) i$DEFINITION)
	names(kos_d)=sapply(kos %>% unlist(recursive = FALSE), function(i) i$ENTRY)
	
	df$ko_desc=sub(" \\[EC:.+","",kos_d)
	df$ec=sub("^.+\\[EC:","EC:",kos_d) %>% gsub("\\]","",.)
	df$ec[!grepl('EC:',df$ec)]=NA
	df$nr_blastp=nrBBpHs$nr_blastp[match(df$prot_id,nrBBpHs$qaccver)]
	df=df[,colnames(dann)]
	return(df)
}

# KEGG enrichment and analysis
prepare_data_for_kegg_analyses <- function(dds,res,ann){
	require(DESeq2)
	rld=rlog(dds, blind = FALSE)
	expr=assay(rld)
	idx=which(!is.na(ann$ko) & ((ann$prot_id %>% sub('\\.p\\d+','',.)) %in% rownames(expr)))
	gene_ids=(ann$prot_id %>% sub('\\.p\\d+','',.))[idx]
	mat=expr[gene_ids,]
	rownames(mat)=ann$ko[idx]
	return(list(rld=rld,gene_ids=gene_ids,idx=idx,expr_mat=mat,res=res[gene_ids,c(2,4)]))
}

get_transcripts_by_KEGG_pathway <- function(
	ann,
	kegg_expr,
	kgs_sigmet,
	kegg_term,
	counts=FALSE,
	alpha=0.05){
	
	kos=kgs_sigmet[grep(kegg_term,names(kgs_sigmet),ignore.case = TRUE)]
	kos %<>% unlist
	idx=which(
		(kegg_expr$expr_mat %>% rownames %>% sub('\\.\\d+','',.)) %in% kos
		& kegg_expr$res$svalue <= alpha)
	if(counts){
		df=cbind.data.frame(
			kegg_expr$res,
			ann[kegg_expr$gene_ids,1:8],
			kegg_expr$expr_mat)
		cols=c(3,1,2,4:10,11:(11-1+dim(kegg_expr$expr_mat)[2]))
	}
	else{
		df=cbind.data.frame(
			kegg_expr$res,
			ann[kegg_expr$gene_ids,1:8])
		cols=c(3,1,2,4:10)
	}
	rownames(df)=c()
	df[idx,cols]
}


#Generate fasta file from Trinotate data frame
get_fasta_from_trinotate <- function(tri,gids,seqtype=c('n','p'),file='',app=FALSE){
	require(magrittr)
	cols=c(1,ifelse(seqtype=='n',17,18))
	apply(
		tri[tri$gene_id %in% gids,..cols],1,
		function(i) sprintf(">%s [moltype=transcribed_RNA] [tech=TSA]\n%s\n",i[1],i[2])) %>% cat(sep='',file=file,append = app)
}

get_nc_gff_from_trinotate <- function(tri,gids,file='',app=FALSE){
	require(magrittr)
	apply(
		tri[tri$gene_id %in% gids,c(1,17)], 1,
		function(i){
			g=list()
			g[[1]]=sprintf(
				"%s\tTrinotate\tgene\t1\t%d\t.\t+\t.\tID=%s\n",
				i[1],nchar(i[2]),i[1])
			g[[2]]=sprintf(
				"%s\tTrinotate\tmRNA\t1\t%d\t.\t+\t.\tID=%s.mRNA;Parent=%s\n",
				i[1],nchar(i[2]),i[1],i[1])
			g[[3]]=sprintf(
				"%s\tTrinotate\ttranscript\t1\t%d\t.\t+\t.\tID=%s.t1;Parent=%s.mRNA;transcribed_sequence\n\n",
				i[1],nchar(i[2]),i[1],i[1])
			paste0(g %>% unlist)
		}) %>% cat(sep='',file=file,append = app)
}

get_annotations_trinotate <- function(tri,man_ann,gids,file='',app=FALSE){
	require(magrittr)
	apply(
		tri[tri$gene_id %in% gids,],1,
		function(i){
			a=list()
			
			if(is.na(i[['prot_id']])){
				return()
			}
			
			a$prot_id=i[['prot_id']]
			if(i[1] %in% rownames(man_ann)){
				ann=man_ann[i[1],]
				if(!is.na(ann[['annot']])){
					if(grepl('\\^',ann[['annot']])){
						a$product=strsplit(ann$annot,'\\^') %>% unlist %>% .[[1]] %>% sub('[DU]:','',.)
					}
					else{
						a$product=ann$annot
					}
				}
				else{
					a$product=ann$annot_uniprot
				}
				if(!is.na(ann$ko)){
					a$ko=paste(ann$ko,ann$ko_desc,collapse = ' ')
				}
				if(!is.na(ann$ec)){
					a$ec=ann$ec
				}
			}
			if(!is.na(i[['Pfam']])){
				a$pfam=strsplit(i[['Pfam']],'`') %>%
					unlist() %>%
					sapply(.,function(j) sub('\\^.*$','',j)) %>%
					unname %>%
					c(.)
			}
			if(!is.na(i[['omcl']])){
				a$omcl=i[['omcl']]
			}
			a$taxon=2109649
			#print(a)
			res=c(sprintf('%s\tDbxref\ttaxon:%s\n',a$prot_id,a$taxon))
			if(!is.null(a$product)){
				res=c(res,sprintf('%s\tproduct\t%s\n',a$prot_id,a$product))
			}
			if(!is.null(a$pfam)){
				res=c(res,
							sapply(
								a$pfam,
								function(j) sprintf('%s\tDbxref\tPFAM:%s\n',a$prot_id,j)) %>%
								unname
				)
			}
			if(!is.null(a$omcl)){
				res=c(res,sprintf('%s\tDbxref\tOrthoMCL:%s\n',a$prot_id,a$omcl))
			}
			if(!is.null(a$ec)){
				res=c(res,sprintf('%s\tNote\t%s\n',a$prot_id,a$ec))
			}
			if(!is.null(a$ko)){
				res=c(res,sprintf('%s\tNote\tKEGG:%s\n',a$prot_id,a$ko))
			}
			#print(res)
			return(res)
			
		}) %>%
		unlist %>% 
		cat(sep='',file=file,append = app)
}


