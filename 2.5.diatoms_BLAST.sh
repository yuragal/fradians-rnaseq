##############################################################################################
# BLAST annotated transcripts and predicted proteins to Synedra acus genome 
# and other diatom genomes
##############################################################################################

outd=$HOME/axSA_assemblies/drap/stranded/trinity_norm_altogether
tempd=/store/$USER
mkdir -p $outd/trinotate/TD2SAgenome && cd $outd/trinotate/TD2SAgenome

#Blast against Synedra acus genome scaffolds
#These are available at https://www.ebi.ac.uk/ena/data/view/CAAAJI010000000

mkdir -p $tempd && cd $tempd
#ln -s ~/axSA_assemblies/CAAAJI01.fna ./
wget http://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/ca/CAAAJI01.fasta.gz
gzip -d CAAAJI01.fasta.gz
makeblastdb -in CAAAJI01.fasta -dbtype nucl -out CAAAJI01
cd -

qfna=$outd/trinotate/transcripts.fa
qfaa=$outd/trinotate/pfam_blastp_hits.transcripts.fa.transdecoder.pep

d=splits_nt
in=$qfna
if [[ ! -d $d ]]; then
    mkdir -p $d && cd $d
    l=$(echo '(' $(cat $in | grep '>' | wc -l) '/ 36 + 1)*2' | bc)
    cat $in | awk '{if($1~/^>/){printf("\n%s\n",$1)}else{printf("%s",$1)}}' | sed 1d | sed 's/*$//' | split -d -a 3 -l $l -
    cd ..
fi

d=splits
in=$qfaa
if [[ ! -d $d ]]; then
    mkdir -p $d && cd $d
    l=$(echo '(' $(cat $in | grep '>' | wc -l) '/ 36 + 1)*2' | bc)
    cat $in | awk '{if($1~/^>/){printf("\n%s\n",$1)}else{printf("%s",$1)}}' | sed 1d | sed 's/*$//' | split -d -a 3 -l $l -
    cd ..
fi

outfmt='6 std qcovs qcovhsp'
find splits_nt -type f -name x???| xargs -P36 -n1 bash -c  'blastn -query $0 -db '$tempd'/CAAAJI01 -num_threads 1 -max_target_seqs 1 -outfmt "'"$outfmt"'" > $0.vs_sac_scaffolds.blastn.out 2>/dev/null' &
cat splits_nt/*.out > transcripts_vs_sac_scaffolds.blastn.out
find splits -type f -name x??? | xargs -P36 -n1 bash -c  'tblastn -query $0 -db '$tempd'/CAAAJI01 -num_threads 1 -max_target_seqs 1 -outfmt "'"$outfmt"'" > $0.vs_sac_scaffolds.tblastn.out 2>/dev/null' &
cat splits/*.out > proteins_vs_sac_scaffolds.tblastn.out
cp *blastn.out $outd/trinotate/TD2SAgenome

#Blast against genomes of other diatoms
tempd=/store/$USER
mkdir -p $tempd/vs_diatoms.ncbi/splits && cd $tempd
echo fracy phtri thapse | xargs -n1 -I% ln -s $HOME/repo/JGI/NCBI/%.faa %.ncbi.faa
echo fracy phtri thapse | xargs -n1 -I% do makeblastdb -in %.ncbi.faa -dbtype prot -out %.ncbi
cp -r $outd/trinotate/TD2SAgenome/splits ./

echo fracy phtri thaps; do 
    find splits -type f -name x??? | xargs -P36 -n1 bash -c  'blastp -query $0 -db '$p'.ncbi -num_threads 1 -outfmt 11 > $0.vs_'$p'.ncbi.blastp.asn 2>/dev/null'
done

for p in fracy phtri thapse; do
    ls splits/x???.vs_${p}.ncbi.blastp.asn | xargs -n1 -P36 bash -c 'b=${0%%.asn}; blast_formatter -archive $0 -max_target_seqs 1 -outfmt "'"$outfmt"'" > $b.out' &
done
mv ./splits./*.asn ./splits./*.out vs_diatoms.ncbi/splits
cd vs_diatoms.ncbi 
cat splits/*.blastp.out | awk '$11<=1e-20 && $13 >= 50' | sort sort -k1,1 -k11,11g > $outd/trinotate/TD2SAgenome/vs_diatoms.ncbi/three_diatoms.BBpHits.out
ls splits/x0??.*asn | xargs -n1 -P36 bash -c 'b=${0%%.asn}; blast_formatter -archive $0 -max_target_seqs 1 -outfmt "'"$outfmt"'" > $b.out'
cat splits/*out | sort -k1,1 -k11,11g | awk '{if($1==q && $2==s){next}else{print $0;q=$1;s=$2}}' > $outd/trinotate/trinotate.nr.blastp.qcovstat.out
rm -rf $tempd