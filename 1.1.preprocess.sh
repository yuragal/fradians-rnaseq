#Should be run on HPC computing node, not on the head node!
#Before running, copy input files to local file system
#Copy the output files back to data storage when finish

#outd -- is path to data storage directory
#tempd -- path to local directory

outd=/home/ygalachyants/axSA_assemblies/drap
cd reads && mkdir -p raw cleaned normalized && cd -
cd reads/raw && ls -l ../../../trinity/reads/ | cut -d' ' -f 11 | xargs -I% -n1 ln -s % ./
tempd=/store/$USER/drap && cd $tempd && mkdir -p reads/raw reads/cleaned reads/normalized
######################################################################
#Clean reads
######################################################################
export PATH=/home/ygalachyants/tools/assembly/qc/bbmap:$HOME/local/compile/pigz-2.4:$PATH
cd rsync -qL $outd/reads/raw/*bz2 $tempd/reads/raw

#The input file names are:
# 20160805.A-13-1_R1.fastq.bz2
# 20160805.A-13-1_R2.fastq.bz2
# 20160805.A-13-2_R1.fastq.bz2
# 20160805.A-13-2_R2.fastq.bz2
# 20160805.A-3-1_R1.fastq.bz2
# 20160805.A-3-1_R2.fastq.bz2
# 20160805.A-3-2_R1.fastq.bz2
# 20160805.A-3-2_R2.fastq.bz2
# 20160805.A-3-3_R1.fastq.bz2
# 20160805.A-3-3_R2.fastq.bz2
# 20160805.A-3-4_R1.fastq.bz2
# 20160805.A-3-4_R2.fastq.bz2
# 20160805.A-3-6_R1.fastq.bz2
# 20160805.A-3-6_R2.fastq.bz2

#Three-pass filtering was used: 
#   first two passes to trim by adapters and quality,
#   the last pass is to filter out artifacts (Artificial contaminants filtered by JGI.) and phix

cd $tempd/reads/cleaned
for p in `ls $tempd/reads/raw/ | xargs -n1 bash -c 'echo ${0%[12].fastq.bz2}' | uniq`; do
    o=${p##*/}
    cmd="bbduk.sh in=${p}#.fastq.bz2 out=stdout.fq \
        ref=adapters \
        k=27 mink=16 ktrim=r hdist=2 tbo tpe ordered \
        2>$tempd/reads/cleaned/${o}.bbduk.1st.pass \
        | bbduk.sh in=stdin.fq int=t out=stdout.fq \
            ref=adapters \
            k=23 ktrim=r mink=10 hdist=1 hdist2=0 qtrim=rl trimq=10 minlen=25 tbo tpe ordered \
            2>$tempd/reads/cleaned/${o}.bbduk.2nd.pass \
        | bbduk.sh in=stdin.fq int=t out=$tempd/reads/cleaned/${o}#.fq.gz \
            k=31 ref=artifacts,phix ordered cardinality ow \
            2>$tempd/reads/cleaned/${o}.bbduk.3rd.pass"
    echo $cmd
    eval $cmd
done
rsync -qL $tempd/reads/cleaned/* /$outd/reads/cleaned

######################################################################
#Normalize reads to assemble by samples
######################################################################
mkdir $tempd/reads/normalized && cd $tempd/reads/normalized
for p in `ls $tempd/reads/cleaned/*gz | xargs -n1 bash -c 'f=${0%[12].fq.gz}; echo ${f##*/}' | uniq`; do
    o=${p#clean_default.}
    cmd="bbnorm.sh in=$tempd/reads/cleaned/$p#.fq.gz out=stdout.fq target=100 2>$o.bbnorm.log \
        | reformat.sh in=stdin.fq int=t out=$o#.fq.gz 2>/dev/null"
    echo $cmd
    eval $cmd
done
rsync -qL $tempd/reads/normalized/* $outd/reads/normalized

rm -rf $tempd/reads


######################################################################
#Merge raw reads, clean and normalize altogether
######################################################################
mkdir -p $tempd/reads/raw $tempd/reads/altogether
rsync --progress -L $outd/reads/raw/20160805.A-* raw/

ls $tempd/reads/raw/*bz2 \
    | xargs -n1 bash -c 'f=${0%[12].fastq.bz2}; echo ${f##*/}' \
    | uniq \
    | xargs -n1 bash -c 'reformat.sh in='$tempd'/reads/raw/${0}#.fastq.bz2 out='$tempd'/reads/altogether/raw.R#.fq.gz app=t' \
        2>$tempd/reads/altogether/reformat.log &
cd $tempd/reads/altogether
cmd="bbduk.sh in=raw.R#.fq.gz out=stdout.fq \
    ref=adapters \
    k=27 mink=16 ktrim=r hdist=2 tbo tpe ordered \
    2>clean.bbduk.1st.pass \
    | bbduk.sh in=stdin.fq int=t out=stdout.fq \
        ref=adapters \
        k=23 ktrim=r mink=10 hdist=1 hdist2=0 qtrim=rl trimq=10 minlen=25 tbo tpe ordered \
        2>clean.bbduk.2nd.pass \
    | bbduk.sh in=stdin.fq int=t out=clean.R#.fq.gz \
        k=31 ref=artifacts,phix ordered cardinality ow \
        2>clean.bbduk.3rd.pass"
echo $cmd
eval $cmd &

cmd="bbnorm.sh in=clean.R#.fq.gz out=stdout.fq target=100 2>norm.bbnorm.log \
        | reformat.sh in=stdin.fq int=t out=norm.R#.fq.gz 2>/dev/null"
echo $cmd
eval $cmd &


rsync -qr $tempd/reads/altogether $outd/reads/ &
