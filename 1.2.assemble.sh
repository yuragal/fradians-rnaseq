
######################################################################
#Run Drap AM Workflow
######################################################################

######################################################################
#Run Assembly part of Drap AM-Trinity Workflow
######################################################################
outd=$HOME/axSA_assemblies/drap
#tempd=/store/$USER/drap
tempd=/store2/$USER/drap
mkdir -p $tempd && cd $tempd
source $HOME/tools/assembly/drap/drap.rc
s=(`ls $outd/reads/normalized/*gz | xargs -n1 bash -c 'f=${0%_R[12].fq.gz}; echo ${f##*/20160805\.A-}' | uniq | xargs`)
cmd=''
for (( i=0;i<7;i++ )); do
    p=${s[$i]}
    cmd="${cmd}runDrap  \
        --R1 $outd/reads/normalized/20160805.A-${p}_R1.fq.gz \
        --R2 $outd/reads/normalized/20160805.A-${p}_R2.fq.gz \
        --dbg trinity \
        -s RF \
        --outdir $tempd/trinity_${p} \
        --no-norm --no-trim \
        --dbg-mem 96"
    cmd="$cmd &> runDrap.trinity_${p}.log\n"
done
echo -e $cmd | bash &

######################################################################
#Run Assembly part of Drap AM-Oases Workflow
######################################################################
s=(`ls $outd/reads/normalized/*gz | xargs -n1 bash -c 'f=${0%_R[12].fq.gz}; echo ${f##*/20160805\.A-}' | uniq | xargs`)
cmd=''
for (( i=0;i<7;i++ )); do
    p=${s[$i]}
    cmd="${cmd}runDrap  \
        --R1 $outd/reads/normalized/20160805.A-${p}_R1.fq.gz \
        --R2 $outd/reads/normalized/20160805.A-${p}_R2.fq.gz \
        --dbg oases -kmer 25,31,37,43,49 \
        --s RF \
        --outdir $tempd/oases_${p} \
        --no-norm --no-trim \
        --dbg-mem 96"
    cmd="$cmd &> runDrap.oases_${p}.log\n"
done
echo -e $cmd | bash &

######################################################################
#Merge Drap Assemblies
######################################################################
dbg=trinity
all_dirs=`echo ${s[@]} | xargs -n1 printf "$tempd/${dbg}_%s," | awk '{print substr($1,1,length($1)-1)}'`
cmd="runMeta \
    -s RF \
    --drap-dirs $all_dirs \
    --outdir $tempd/metaRF_${dbg}"
eval "$cmd &> runMeta.RF.$dbg.log" &

dbg=oases
all_dirs=`echo ${s[@]} | xargs -n1 printf "$tempd/${dbg}_%s," | awk '{print substr($1,1,length($1)-1)}'`
cmd="runMeta \
    -s RF \
    --drap-dirs $all_dirs \
    --outdir $tempd/metaRF_${dbg}"
eval "$cmd &> runMeta.RF.$dbg.log" &


######################################################################
#Run Drap MA Workflow
######################################################################
outd=$HOME/axSA_assemblies/drap
tempd=/store/$USER/drap
mkdir -p $tempd/reads/altogether && rsync -q $outd/reads/altogether/n*.R[12].fq.gz $tempd/reads/altogether
source $HOME/tools/assembly/drap/drap.rc
p=norm
as=trinity && asl='--dbg trinity'
cmd="runDrap  \
    --R1 $tempd/reads/altogether/${p}.R1.fq.gz \
    --R2 $tempd/reads/altogether/${p}.R2.fq.gz \
    $asl \
    --strand RF \
    --outdir $tempd/${as}_${p}_altogether \
    --no-norm --no-trim \
    --dbg-mem 96"
eval "$cmd &> runDrap.${as}_${p}_altogether.log" &

as=oases && asl='--dbg oases --kmer 25,31,37,43,49'
cmd="runDrap  \
    --R1 $tempd/reads/altogether/${p}.R1.fq.gz \
    --R2 $tempd/reads/altogether/${p}.R2.fq.gz \
    $asl \
    --strand RF \
    --outdir $tempd/${as}_${p}_altogether \
    --no-norm --no-trim \
    --dbg-mem 96"
eval "$cmd &> runDrap.${as}_${p}_altogether.log" &
