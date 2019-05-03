##############################################################################################
# These commands are to submit jobs to queue
# BLAST against NCBI NR is required to run on HPC
##############################################################################################

db=$HOME/repo/NCBI/2016_Nov/nr.tar.bz2
export LD_LIBRARY_PATH=/share/apps/gcc/7.3/lib64:$LD_LIBRARY_PATH

outd=/home/ygalachyants/axSA_assemblies/drap/stranded/trinity_norm_altogether
tempd=/store/$USER
qfna=$outd/trinotate/transcripts.fa
qfaa=$outd/trinotate/pfam_blastp_hits.transcripts.fa.transdecoder.pep

if [[ ! -d $d/splits ]]; then \
    mkdir -p $d/splits && cd $d/splits; \
    cat ../$qnfa| awk '{if($1~/^>/){printf("\n%s\n",$1)}else{printf("%s",$1)}}' | sed 1d | sed 's/*$//' | split -d -a 3 -l 2000 -; \
    cd ..; \
fi

if [[ ! -d $d/splits_nt ]]; then \
    mkdir -p $d/splits_nt && cd $d/splits_nt; \
    cat ../$qfaa| awk '{if($1~/^>/){printf("\n%s\n",$1)}else{printf("%s",$1)}}' | sed 1d | sed 's/*$//' | split -d -a 3 -l 2000 -; \
    cd ..; \
fi

#In my settings blastx 1000 sequences against nr takes ~12-16 walltime hours at Intel 36-core node
#split dir -- splits or splits_nt
spld=splits
#blast program -- blastp or blastx
bprog=blastp
seq -f 'x%03g' 0 22 \
    | xargs -n1 -I FILE echo 'set -vex; cd \\$PBS_O_WORKDIR; prepare_db.sh '$db'; '$bprog' -query '$spld'/FILE -db /store/ygalachyants/nr/nr -num_threads 36 -max_target_seqs 100 -evalue 1e-5 -outfmt 11 2> '$spld'/FILE.'$bprog'.log | gzip > '$spld'/FILE.nr.'$bprog'.asn.gz; rm -rf /store/ygalachyants' \
    | xargs -L1 -I CMD bash -c 'echo "CMD" | /share/system/torque/master/bin/qsub -q reserve -l nodes=1:ppn=36 -N NR.'$bprog' -V; sleep 300'

spld=splits_nt
#blast program -- blastp or blastx
bprog=blastx
seq -f 'x%03g' 0 28 \
    | xargs -n1 -I FILE echo 'set -vex; cd \\$PBS_O_WORKDIR; prepare_db.sh '$db'; '$bprog' -query '$spld'/FILE -db /store/ygalachyants/nr/nr -num_threads 36 -max_target_seqs 100 -evalue 1e-5 -outfmt 11 2> '$spld'/FILE.'$bprog'.log | gzip > '$spld'/FILE.nr.'$bprog'.asn.gz; rm -rf /store/ygalachyants' \
    | xargs -L1 -I CMD bash -c 'echo "CMD" | /share/system/torque/master/bin/qsub -q reserve -l nodes=1:ppn=36 -N NR.'$bprog' -V; sleep 300'

#After jobs complete, merge results
#On allocated node:
db=$HOME/repo/NCBI/2016_Nov/nr.tar.bz2
outd=$HOME/axSA_assemblies/drap/stranded/trinity_norm_altogether
tempd=/store/$USER
prepare_db.sh $db &> prepare_db.log
cd $tempd
mkdir nr_res && cd nr_res 
cp $outd/trinotate/splits/x???.nr.blastp.asn.gz ./
cp $outd/trinotate/splits_nt/x???.nr.blastx.asn.gz ./
ls x???*gz | xargs -P36 -n1 -I% gzip -d %
seq -f 'x%03g' 0 22 | xargs -P36 -n1 bash -c 'blast_formatter -archive ${0}.nr.blastp.asn -outfmt 6 -max_target_seqs 1 > ${0}.nr.blastp.out' &
seq -f 'x%03g' 0 28 | xargs -P36 -n1 bash -c 'blast_formatter -archive ${0}.nr.blastx.asn -outfmt 6 -max_target_seqs 1 > ${0}.nr.blastx.out' &
cat x???.nr.blastp.out > trinotate.nr.blastp.out
cat x???.nr.blastx.out > trinotate.nr.blastx.out
cp trinotate.nr.blast?.out $outd/trinotate

cat splits_nt/x???.blastx.out > trinotate.blastx.out
cat splits/x???.blastp.out > trinotate.blastp.out
cat splits/x???.pfam.out > trinotate.pfam.out
rm -rf splits
