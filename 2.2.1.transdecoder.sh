##############################################################################################
# Run Transdecoder with DRAP assembly
##############################################################################################

#export SIGNALP=$HOME/tools/annotation/signalp-4.1
#export PATH=$HOME/tools/annotation/signalp-4.1:$PATH
#export PATH=$HOME/tools/annotation/rnammer:$PATH
export PATH=$HOME/tools/annotation/tmhmm-2.0c/bin:$PATH

signalp=$HOME/tools/annotation/signalp-4.1/signalp
rnammer=$HOME/tools/annotation/rnammer/rnammer
tmhmm=$HOME/tools/annotation/tmhmm-2.0c/bin/tmhmm

outd=$HOME/axSA_assemblies/drap/stranded/trinity_norm_altogether
tempd=/store/$USER

TD=$HOME/tools/annotation/TransDecoder-v5.0.2
trinotate_dir=$HOME/tools/annotation/Trinotate-v3.1.1

mkdir $outd/trinotate && cd $outd/trinotate
ln -s ../f-rmbt_filtering/transcripts_fpkm_1.fa transcripts.fa

#TD with no hits from Pfam and uniprot_sprot
$TD/TransDecoder.LongOrfs -t ./transcripts.fa -S &> TD.longorfs.log &
$TD/TransDecoder.Predict -t ./transcripts.fa &> TD.predict.log &
prename 's/^/no_hits./' transcripts.fa.transdecoder.*

#TD with hits from Pfam
#Can be run on single computing node for reasonable time

qfna=$outd/trinotate/transcripts.fa
qfaa=$outd/trinotate/pfam_blastp_hits.transcripts.fa.transdecoder.pep

prepare_db.sh $trinotate_dir/Pfam-A.tar.bz2
prepare_db.sh $trinotate_dir/uniprot_sprot.tar.bz2

if [[ ! -d splits ]]; then \
    mkdir splits && cd splits; \
    cat $qfaa | awk '{if($1~/^>/){printf("\n%s\n",$1)}else{printf("%s",$1)}}' | sed 1d | sed 's/*$//' | split -d -a 3 -l 500 -; \
    cd ..; \
fi

if [[ ! -d splits_nt ]]; then
    mkdir splits_nt && cd splits_nt
    cat $qfna | awk '{if($1~/^>/){printf("\n%s\n",$1)}else{printf("%s",$1)}}' | sed 1d | sed 's/*$//' | split -d -a 3 -l 500 -
    cd ..
fi

find splits -type f -name x??? | head -36 | xargs -P36 -n1 bash -c 'hmmscan --cpu 1 --domtblout ${0}.pfam.out '$tempd'/Pfam-A.hmm ${0} &> ${0}.pfam.log'
d=$outd/trinotate/transcripts.fa.transdecoder_dir
if [[ ! -d $d/splits ]]; then \
    mkdir -p $d/splits && cd $d/splits; \
    cat ../longest_orfs.pep | awk '{if($1~/^>/){printf("\n%s\n",$1)}else{printf("%s",$1)}}' | sed 1d | sed 's/*$//' | split -d -a 3 -l 500 -; \
    cd ..; \
fi

(cat splits/x000.pfam.out | grep '^#'; cat splits/x???.pfam.out | grep -v '^#') > pfam.hits.out

t=$(stat ./transcripts.fa.transdecoder_dir.__checkpoints/orf_select.ok | grep Modify | cut -d' ' -f3 | cut -d'.' -f1)
find ./transcripts.fa.transdecoder_dir.__checkpoints/ -type f -newermt $t -exec rm {} \;
$TD/TransDecoder.Predict --retain_pfam_hits $d/pfam.hits.out -t ./transcripts.fa &> TD.predict.pfam.log &
prename 's/^/pfam_hits./' transcripts.fa.transdecoder.*

#TD with hits from uniprot_sprot
db=uniprot_sprot
find splits -type f -name x??? | xargs -P36 -n1 bash -c \
    'blastp -query $0 -db '$trinotate_dir/$db' -num_threads 1 -max_target_seqs 1 -outfmt 6 > $0.blastp.out 2>/dev/null'
find splits_nt -type f -name x??? | xargs -P36 -n1 bash -c \
    'blastx -query $0 -db '$trinotate_dir/$db' -num_threads 1 -max_target_seqs 1 -outfmt 6 > $0.blastx.out 2> /dev/null'
cat splits/x???.blastp.out > blastp.hits.out
cat splits_nt/x???.blastp.out > blastx.hits.out

t=$(stat ./transcripts.fa.transdecoder_dir.__checkpoints/orf_select.ok | grep Modify | cut -d' ' -f3 | cut -d'.' -f1)
find ./transcripts.fa.transdecoder_dir.__checkpoints/ -type f -newermt $t -exec rm {} \;
$TD/TransDecoder.Predict --retain_pfam_hits $d/pfam.hits.out --retain_blastp_hits $d/blastp.hits.out -t ./transcripts.fa &> TD.predict.pfam_blast.log &
prename 's/^/pfam_blastp_hits./' transcripts.fa.transdecoder.*

rm -rf splits splits_nt
#Get stats on ORFs found
#Seems that pfam+blastp has the best performance
echo no pfam pfam_blastp \
    | xargs -n1 bash -c 'l=$(cat $0_hits.transcripts.fa.transdecoder.pep | grep -v ">" | tr -d "\n" | wc -c); echo "$l -- length of $0 hits TD"' \
&& echo no pfam pfam_blastp \
    | xargs -n1 bash -c 'c=$(cat $0_hits.transcripts.fa.transdecoder.pep | grep ">" | cut -d" " -f1 | cut -c2- | wc -l); echo "$c -- count of $0 transcripts TD"' \
&& echo no pfam pfam_blastp \
    | xargs -n1 bash -c 'c=$(cat $0_hits.transcripts.fa.transdecoder.pep | grep ">" | cut -d" " -f1 | cut -c2- | grep -v ".p1$" | wc -l); echo "$c -- count of .p2+ $0 transcripts TD"'

