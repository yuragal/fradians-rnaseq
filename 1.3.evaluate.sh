######################################################################
#It seems, Trinity meta-assembly is better than Oases
#From altogether assemblies, trinity is MUCH better than Oases
#See TRANSRATE ASSEMBLY SCORE overall statistics and other figures
######################################################################
outd=$HOME/axSA_assemblies/drap
grep -B40 'p good contigs' metaRF_*/00-META-ASSEMBLY_RATING  | less
grep -B40 'p good contigs' *norm_altogether/00-ASSEMBLY_RATING  |less
grep -A5 'TRANSRATE ASSEMBLY SCORE' metaRF_*/00-META-ASSEMBLY_RATING  | less
grep -A5 'TRANSRATE ASSEMBLY SCORE' *norm_altogether/00-ASSEMBLY_RATING

######################################################################
#Compare the merged trinity assembly with altogether
######################################################################

bdir=$HOME/tools/annotation/busco
pyt=/share/apps/python/3.6.4/bin/python3
for d in trinity_norm_altogether metaRF_trinity; do
    for bdb in eukaryota_odb9 protists_ensembl; do
        bdb_p=$(echo $bdb | cut -d'_' -f1)
        cmd="$pyt $bdir/scripts/run_BUSCO.py \
            -i $outd/$d/transcripts_fpkm_1.fa \
            -m transcriptome \
            -o busco_${bdb_p}_${d} \
            -l $bdir/$bdb \
            -t ./tmp_${bdb_p}_$d \
            -c 9 -f > busco.${bdb_p}_${d}.log 2>&1 &"
        eval $cmd
    done
done
rm -rf tmp_*
grep -A2 'Results:' busco.*.log
cp -r run_busco_* busco.*.log $outd/

as=oases
outd=$HOME/axSA_assemblies/drap
for d in ${as}_norm_altogether metaRF_${as}; do
    for bdb in eukaryota_odb9 protists_ensembl; do
        bdb_p=$(echo $bdb | cut -d'_' -f1)
        cmd="${bdir}/scripts/run_BUSCO.py \
            -i $outd/$d/transcripts_fpkm_1.fa \
            -m transcriptome \
            -o busco_${bdb_p}_${d} \
            -l $bdir/$bdb \
            -t ./tmp_${bdb_p}_$d \
            -c 9 -f > busco.${bdb_p}_${d}.log 2>&1 &"
        eval "$cmd"
        #echo "$cmd"
    done
done

######################################################################
#Run Assessment between trinity assemblies
######################################################################
ln -s transcripts_fpkm_1.fa $outd/metaRF_trinity/metaRF_trinity_fpkm_1.fa 
ln -s transcripts_fpkm_1.fa $outd/trinity_norm_altogether/trinity_norm_altogether_fpkm_1.fa

cmd="runAssessment \
    --assemblies $outd/metaRF_trinity/metaRF_trinity_fpkm_1.fa,$outd/trinity_norm_altogether/trinity_norm_altogether_fpkm_1.fa \
    --R1 $tempd/reads/altogether/norm.R1.fq.gz \
    --R2 $tempd/reads/altogether/norm.R2.fq.gz \
    --busco-lineage eukaryota_odb9 \
    --outdir $tempd/trinity_assessment"
eval "$cmd &> runAssessment.trinity.log" &
#To restart assessment after failure, run eval again once error is fixed

