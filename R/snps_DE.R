library(kissDE)
#data_dir='~/matrosov/store/nodes/sm122/ygalachyants/'
data_dir='~/matrosov//home/ygalachyants/axSA_assemblies/drap/stranded/trinity_norm_altogether/snps/'
snps=kissplice2counts(paste0(data_dir,'results/results_k41_coherents_type_0a.fa'), pairedEnd=TRUE )
cond=c(rep('A6',2),rep('BK280',5))
snp_res=diffExpressedVariants(snps, cond, pvalue=1)

t=snp_res$finalTable[which(snp_res$finalTable$Adjusted_pvalue <=0.05),]
colnames(t)[3:16]%<>%gsub("Variant","V",.) %>% gsub("repl","r",.) %>% gsub("BK280","BK",.) %>% gsub("_Norm","",.)
View(t)
write.table(snp_res$finalTable, file=paste0(data_dir,'kissDE_output'), sep="\t", quote=FALSE)

snps2ref=snp_get_good_snps(f = paste0(data_dir,'mainOutput.tsv'))
snp_get_distr(annot=ta_clean,res=deseq_results$NExp_DSLT$res,snp=snps2ref, tot=TRUE)


#=====================================================
#280 only
#=====================================================
data_dir='~/matrosov//home/ygalachyants/axSA_assemblies/drap/stranded/trinity_norm_altogether/snps/'
snps=kissplice2counts(paste0(data_dir,'results_280/results_k41_coherents_type_0a.fa'), pairedEnd=TRUE )
cond=c(rep('dark',2),rep('light',3))
snp_res_280=diffExpressedVariants(snps, cond, pvalue=1)
t=snp_res_280$finalTable[which(snp_res_280$finalTable$Adjusted_pvalue <=0.05),]

colnames(t)[3:12] %<>% gsub("Variant","V",.) %>% gsub("repl","r",.) %>% gsub("light","L",.) %>% gsub("dark","D",.) %>% gsub("_Norm","",.)
#View(t)
write.table(snp_res_280$finalTable, file=paste0(data_dir,'kissDE_280_output'), sep="\t", quote=FALSE)
#Now run KissDE script
#And load results

snps2ref_280=snp_get_good_snps(f = paste0(data_dir,'mainOutput_280.tsv'))
snp_get_distr(annot=ta_clean,res=deseq_results$Dark_Light$res,snp=snps2ref_280,tot=TRUE)



ES_tr=snp_get_distr(annot=ta_clean,res=deseq_results$NExp_DSLT$res,snp=snps2ref,count=FALSE,tot=TRUE)
DL_tr=snp_get_distr(annot=ta_clean,res=deseq_results$Dark_Light$res,snp=snps2ref_280,count=FALSE,tot=TRUE)
ES_tr$`C+S+Syn`$NDE %in% DL_tr$`C+S+Syn`$NDE %>% which %>% length()
ES_tr$`C+S+Nsyn`$NDE %in% DL_tr$`C+S+Nsyn`$NDE %>% which %>% length()
ES_tr$`C+S+Syn`$DE %in% DL_tr$`C+S+Syn`$DE %>% which %>% length()
ES_tr$`C+S+Nsyn`$DE %in% DL_tr$`C+S+Nsyn`$DE %>% which %>% length()


which(!snps2ref$Is_in_CDS) %>% length
which(snps2ref$Is_in_CDS & snps2ref$Is_not_synonymous) %>% length
which(snps2ref$Is_in_CDS & !snps2ref$Is_not_synonymous) %>% length

which(!snps2ref_280$Is_in_CDS) %>% length
which(snps2ref_280$Is_in_CDS & snps2ref_280$Is_not_synonymous) %>% length
which(snps2ref_280$Is_in_CDS & !snps2ref_280$Is_not_synonymous) %>% length

a=which(snps2ref$Is_in_CDS & snps2ref$Is_not_synonymous)
b=which(snps2ref_280$Is_in_CDS & snps2ref_280$Is_not_synonymous)

aa=a[snps2ref$X.Component_ID[a] %in% snps2ref_280$X.Component_ID[b]]
bb=b[snps2ref_280$X.Component_ID[b] %in% snps2ref$X.Component_ID[a]]

#snps2ref$X.Component_ID[ aa ] %>% sort
#snps2ref_280$X.Component_ID[ bb ] %>% sort

View(snps2ref[ aa , c(1,5:9)])
View(snps2ref_280[ bb , c(1,5:9)])

View(snps2ref[ aa , ])
View(snps2ref_280[ bb , ])

