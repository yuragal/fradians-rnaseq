library(KEGGREST)
library(magrittr)

args=commandArgs(TRUE) %>% as.numeric()

start=args[1]
end=args[2]
intstep=args[3]
# print(str(args))
# 
# print(seq(start,end, by=intstep))
# quit()

load("kegg_accs.RData")

kgids=lapply(seq(start,end, by=intstep), function(i) {
	#cat(start, 'keggConv', date(),i,"\n")
	cat(sprintf("%10s KeggConv, %25s: %03.0f%% done\n", start, date(), (i-start)/(end-start)*100))
	#print(paste0("ncbi-proteinid:",saccs[i:(i+intstep-1)]))
	keggConv("genes", paste0("ncbi-proteinid:",saccs[i:(i+intstep-1)]))
}) %>% unlist
kgs=sapply(seq(1,kgids %>% length(), by=10), function(i) {
	cat(sprintf("%10s KeggGet, %25s: %03.0f%% done\n", start, date(), i/(end-start)*100))
	#cat(start, 'keggGet', date(),i,"\n")
	keggGet(kgids[i:(i+9)]) %>% c()
})
#names(kgs)=kgids

save(kgids,kgs,file=sprintf("kegg_genes.%d.RData",start))