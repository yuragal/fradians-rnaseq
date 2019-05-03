#export SIGNALP=$HOME/tools/annotation/signalp-4.1
#export PATH=$HOME/tools/annotation/signalp-4.1:$PATH
#export PATH=$HOME/tools/annotation/rnammer:$PATH
export PATH=$HOME/tools/annotation/tmhmm-2.0c/bin:$PATH

signalp=$HOME/tools/annotation/signalp-4.1/signalp
rnammer=$HOME/tools/annotation/rnammer/rnammer
tmhmm=$HOME/tools/annotation/tmhmm-2.0c/bin/tmhmm

trinotate_dir=$HOME/tools/annotation/Trinotate-v3.1.1

outd=/home/ygalachyants/axSA_assemblies/drap/stranded/trinity_norm_altogether
tempd=/store/$USER

qfna=$outd/trinotate/transcripts.fa
qfaa=$outd/trinotate/pfam_blastp_hits.transcripts.fa.transdecoder.pep


##############################################################################################
# This is to treat each resulted DRAP contig as unitig having single isoform
# Seems to be more appplicable for SA-RNAseq case as:
#!!! By using pfam and blastp hits with TD, number of ORFs per gene become 1. 
#!!! No contigs are now with multiple ORFS.
#!!!   This fact is in accordance with DRAP-algorithm as DRAP is designed to decrease redundancy 
#!!!   and generate assemblies with 1:1 gene-to-transcripts relationship.
#!!!   See https://peerj.com/articles/2988/#p-11 for details.
##############################################################################################
cd $outd/trinotate
cat $qfna \
    | grep  '>' \
    | cut -d' ' -f1 \
    | cut -c2- \
    | awk '{print $1 "\t" $1}' > tmp.map 
mv tmp.map transcripts.gene2tr.map

#This is to test gene2tr.map was correctly generated.
#join <(cat transcripts.gene2tr.map | sort -k1,1b) <(cat transcripts.fa | awk '{if(/>/){printf("\n%s\t",$1)}else{printf("%s",$1)}}' | cut -c2- | sort -k1,1b) | head

$signalp -f short -n signalp.out $qfaa &> signalp.log &
$tmhmm --short < $qfaa > tmhmm.out 2> tmhmm.log &
$trinotate_dir/util/rnammer_support/RnammerTranscriptome.pl --transcriptome $qfna --path_to_rnammer $rnammer &>rnammer.log &

$trinotate_dir/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate no_cleanup_flag
cp -f $trinotate_dir/Trinotate.sqlite $tempd && chmod u+w $tempd/Trinotate.sqlite

$trinotate_dir/Trinotate $tempd/Trinotate.sqlite init --gene_trans_map transcripts.gene2tr.map --transcript_fasta $qfna --transdecoder_pep $qfaa
$trinotate_dir/Trinotate $tempd/Trinotate.sqlite LOAD_swissprot_blastp trinotate.blastp.out
$trinotate_dir/Trinotate $tempd/Trinotate.sqlite LOAD_swissprot_blastx trinotate.blastx.out
$trinotate_dir/Trinotate $tempd/Trinotate.sqlite LOAD_pfam trinotate.pfam.out
$trinotate_dir/Trinotate $tempd/Trinotate.sqlite LOAD_tmhmm tmhmm.out
$trinotate_dir/Trinotate $tempd/Trinotate.sqlite LOAD_signalp signalp.out
$trinotate_dir/Trinotate $tempd/Trinotate.sqlite LOAD_rnammer transcripts.fa.rnammer.gff
$trinotate_dir/Trinotate $tempd/Trinotate.sqlite LOAD_custom_blast --outfmt6 trinotate.nr.blastp.out --prog blastp --dbtype nr
$trinotate_dir/Trinotate $tempd/Trinotate.sqlite LOAD_custom_blast --outfmt6 trinotate.nr.blastx.out --prog blastx --dbtype nr
$trinotate_dir/Trinotate $tempd/Trinotate.sqlite report --incl_pep --incl_trans > trinotate.annotation_report.csv
chmod a-w $tempd/Trinotate.sqlite
cp $tempd/Trinotate.sqlite ./


#Format of trinotate.annotation_report.csv
#Each line has 16 fields separated by tab
# 1         #gene_id
# 2         transcript_id
# 3         sprot_Top_BLASTX_hit
# 4         RNAMMER
# 5         prot_id
# 6         prot_coords
# 7         sprot_Top_BLASTP_hit
# 8         Pfam
# 9         SignalP
# 10        TmHMM
# 11        eggnog
# 12        Kegg
# 13        gene_ontology_blast
# 14        gene_ontology_pfam
# 15        transcript
# 16        peptide

#View this by the following command to have fileds separated by newlines and enumerated
cat trinotate.annotation_report.csv | sed 's/\t/\n\t/g' | awk '{if(NR%16==0){n=16}else{n=NR%16} print n " " $0}' | less
cat trinotate.annotation_report.csv | sed 's/\t/\n\t/g' | awk -vnf=16 '{if(NR%nf==0){n=nf}else{n=NR%nf} print n " " $0}' | less

$trinotate_dir/util/Trinotate_get_feature_name_encoding_attributes.pl trinotate.annotation_report.csv > trinotate.annot_feature_map.txt
$trinotate_dir/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls trinotate.annotation_report.csv -G --include_ancestral_terms > trinotate.go_annotations.txt
$trinotate_dir/Analysis/DifferentialExpression/run_GOseq.pl --factor_labeling  factor_labeling.txt --GO_assignments go_annotations.txt --lengths gene.lengths.txt

