##############################################################################################
# Explore the annotated genes
##############################################################################################

outd=$HOME/axSA_assemblies/drap/stranded/trinity_norm_altogether

cd $outd/trinotate

ta=$outd/trinotate/trinotate.annotation_report.csv
cont=$outd/trinotate/filter_primates/primate.contig.ids
tr=$outd/trinotate/TD2SAgenome/transcripts_vs_sac_scaffolds.blastn.out
pr=$outd/trinotate/TD2SAgenome/proteins_vs_sac_scaffolds.blastn.out

#Annotated contigs! Do not account for SignalP and TmHMM hits!!!
cat ./trinotate.annot_feature_map.txt | sed 's/\^Tm[0-9]\+//;s/\^sigP//' | awk 'length($1)!=length($2)' | wc -l
# 15714

# TD peptides with hits in one of blastx #3, blastp#7 or pfam #8 analyses
#There are 1516 TD peptides that were not annotated; The rest 14198 were annotated
cat $ta | sed 1d | awk -F'\t' '$5!="." && ($3!="." || $7!="." || $8!="."){print $5}' | sort | uniq | wc -l
# 14198

#Out of 14198 ORFs annotated by Trinotate as proteins, there are 11056 ORFs that have highly-significant (e-value below  1e-10, query coverage per subject above 70%) matches with sequences of filtered gene models of Thaps3 and Phatr2 genomes.
cat splits/x???.vs_diatoms_proteins.blastp.out | awk '$11<=1e-10 && $13>70 {print $1}' | sort | uniq | wc -l 
# 11056
# Their annotation:
join -t $'\t' -1 1 -2 5 <(cat splits/x???.vs_diatoms_proteins.blastp.out | awk '$11<=1e-10 && $13>70 {print $1}' | sort | uniq) <(cat $ta | sed 1d | awk -F'\t' '$5!="."') | less

comm -1 -2 <(cat $ta | sed 1d | awk -F'\t' '$5!="." && ($3!="." || $7!="." || $8!="."){print $1}' | sort | uniq) <(cat $tr | awk '$11<=1e-50 && $13>=90 {print $1}' | sort | uniq) | wc -l
# 10554
comm -1 -2 <(cat $ta | sed 1d | awk -F'\t' '$5!="." && ($3!="." || $7!="." || $8!="."){print $1}' | sort | uniq) <(cat $pr | awk '$11<=1e-50 && $13>=90 {print substr($1,1,length($1)-3)}' | sort | uniq) | wc -l
# 11132

comm -1 -2 <(cat $ta | sed 1d| grep -vf $cont | awk -F'\t' '{print $1}' | sort | uniq) <(cat $tr | awk '$11<=1e-50 && $13>=90 {print $1}' | sort | uniq) | wc -l
#21140
comm -1 -2 <(cat $ta | sed 1d | grep -vf $cont | awk -F'\t' '$5!="." {print $1}' | sort | uniq) <(cat $tr | awk '$11<=1e-50 && $13>=90 {print $1}' | sort | uniq) | wc -l
#17634
comm -1 -2 <(cat $ta | sed 1d | grep -vf $cont | awk -F'\t' '$5!="." {print $1}' | sort | uniq) <(cat $pr  | awk '$11<=1e-50 && $13>=90 {print substr($1,1,length($1)-3)}' | sort | uniq) | wc -l
#17916
cat $ta | sed 1d | grep -vf $cont | awk -F'\t' '$5!="." && ($3!="." || $7!="." || $8!=".")' | wc -l
#13526
comm -1 -2 <(cat $ta | sed 1d | grep -vf $cont | awk -F'\t' '$5!="." && ($3!="." || $7!="." || $8!=".") {print $1}' | sort | uniq) <(cat $pr | awk '$11<=1e-50 && $13>=90 {print substr($1,1,length($1)-3)}' | sort | uniq) | wc -l
#11065

comm -1 -2 <(cat $ta | sed 1d | grep -vf $cont | awk -F'\t' '$5!="." {print $1}') <(cat splits/x???.vs_diatoms_proteins.blastp.out | awk '$11<=1e-10 && $13>=70 {print substr($1,1,length($1)-3)}' | sort | uniq) | wc -l
#10966
comm -1 -2 <(cat $ta | sed 1d | grep -vf $cont | awk -F'\t' '$5!="." && ($3!="." || $7!="." || $8!="."){print $1}') <(cat splits/x???.vs_diatoms_proteins.blastp.out | awk '$11<=1e-10 && $13>=70 {print substr($1,1,length($1)-3)}' | sort | uniq) | wc -l
#8563
comm -1 -2 <(comm -1 -2 <(cat $ta | sed 1d | grep -vf $cont | awk -F'\t' '$5!="." && ($3!="." || $7!="." || $8!="."){print $1}' | sort | uniq) <(cat $pr | awk '$11<=1e-50 && $13>=90 {print substr($1,1,length($1)-3)}' | sort | uniq)) <(cat splits/x???.vs_diatoms_proteins.blastp.out | awk '$11<=1e-10 && $13>=70 {print substr($1,1,length($1)-3)}' | sort | uniq) | wc -l
#7172

#Number of ORFs 
#with no annotatin either
cat $ta | sed 1d | grep -vf $cont | awk -F'\t' '$5!="." && $15=="." && $16=="."' | wc -l
#10661
#with no annotation and no hits in UniProt/Swissprot/Pfam but having hits against NR
cat $outd/trinotate/trinotate.annotation_report.nr.csv | sed 1d | grep -vf $cont | awk -F'\t' '$5!="." && $15=="." && $16=="." && ($3=="." || $7=="." || $10==".") && ($8!="." || $9!=".")' | wc -l
#7272

#Select Ex95 contigs with proteins for OrthoMCL clustering
#!!! See ../abundance_estimation folder for details
#There are 11102 such ORFs in total
cat ../abundance_estimation/salmon.isoform.TMM.EXPR.matrix.E-inputs | sed 1d | awk '$1<=95{print $2}' | wc -l
# 12446
join \
    <(cat ../abundance_estimation/salmon.isoform.TMM.EXPR.matrix.E-inputs \
        | sed 1d \
        | awk '$1<=95{print $2}' \
        | sort) \
    <(cat trinotate.annotation_report.csv \
        | sed 1d \
        | awk -F'\t' '$5!="."{print $1 "\t" $16}' \
        | sort -k1,1) \
    | awk '{print ">" $1 "\n" $2}' > trinotate.Ex95.faa

cat trinotate.Ex95.faa | grep '>' | wc -l
# 11102

# Out of 12446 Ex95 contigs, there are 11102 protein-coding contigs and 1344 contigs without ORFs. Among the protein-coding S.acus genes, 8125 are homologous (have orthology/outparalogy relationship) to Phatr2 or Thaps3 genes.
