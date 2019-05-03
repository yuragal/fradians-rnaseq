##############################################################################################
# Filter contaminants
##############################################################################################

outd=$HOME/axSA_assemblies/drap/stranded/trinity_norm_altogether
tempd=/store/$USER

ncbi_tax=$HOME/repo/NCBI/taxonomy/2018_May
a2t=acc2taxid

cd $tempd
#NCBI taxonomy, prot.accession2taxid, and dead_prot.accession2taxid 
#should be downloaded from ftp://ftp.ncbi.nih.gov/pub/taxonomy.
#Next one need to prepate sqlite database like so
#####################################################
#SQLite commands:
#####################################################
sqlite3 accession2taxid.sqlite

CREATE TABLE prot_accession2taxid (
    accession   VARCHAR(128) NOT NULL,
    acc_ver     VARCHAR(128) NOT NULL,
    taxon_id    INT(10) NOT NULL,
    gi          INT(10),
    PRIMARY KEY (accession)
);

CREATE TABLE dead_prot_accession2taxid (
    accession   VARCHAR(128) NOT NULL,
    acc_ver     VARCHAR(128) NOT NULL,
    taxon_id    INT(10) NOT NULL,
    gi          INT(10),
    PRIMARY KEY (accession)
);

.separator "\t"
.import prot.accession2taxid prot_accession2taxid
.import dead_prot.accession2taxid dead_prot_accession2taxid


#ALTER TABLE accession2taxid RENAME TO prot_accession2taxid;

Ctrl+D
#####################################################
#End of SQLite session
#####################################################
cp $ncbi_tax/accession2taxid/$a2t.sqlite.bz2 ./
pbzip2 -dc $a2t.sqlite.bz2 > $a2t.sqlite

cd $ncbi_tax
zcat lineages-2018-03-12.csv.gz | grep -i primates > primates.csv
cat primates.csv | cut -d',' -f1 > primates.taxid
zcat lineages-2018-03-12.csv.gz | awk -F',' '$6=="Muridae"{print $1}' > muridae.taxid
zcat lineages-2018-03-12.csv.gz | awk -F',' '$6=="Bovidae"{print $1}' > bovidae.taxid
cd -

db=$HOME/repo/NCBI/2016_Nov/nr.tar.bz2 && prepare_db.sh $db &> prepare_db.log
mkdir -p nr_res && cd nr_res

cp $outd/trinotate/splits/x???.nr.blastp.asn.gz ./
ls x???.nr.bastp.asn.gz | xargs -P36 -n1 -I% gzip -d %
seq -f 'x%03g' 0 22 | xargs -P36 -n1 bash -c 'blast_formatter -archive ${0}.nr.blastp.asn -outfmt "6 std qlen slen qcovs qcovhsp" -max_target_seqs 10 > ${0}.nr.blastp.out' &

cd ..
cat nr_res/x???.nr.blastp.out > nr.blastp.out

cat nr.blastp.out \
    | cut -f2 \
    | sort \
    | uniq \
    | perl -ne 'chomp; $s=$_; $s=~s/\.\d+$//; print "$_\t$s\n"' > ids

cat ids \
    | cut -f2 \
    | xargs -n10000 \
    | xargs -L1 \
    | sed "s/\s/','/g; s/^/select * from prot_accession2taxid where accession in ('/; s/$/');/" > sql
sqlite3 acc2taxid.sqlite < sql > nr.taxids

cat ids \
    | cut -f2 \
    | xargs -n10000 \
    | xargs -L1 \
    | sed "s/\s/','/g; s/^/select * from dead_prot_accession2taxid where accession in ('/; s/$/');/" > sql
sqlite3 acc2taxid.sqlite < sql > nr.dead.taxids

comm -3 <(cut -f2 ids | sort) <(cat nr.taxids | cut -d'|' -f1 | sort) | grep '|' > unmapped.ids
cp unmapped.ids ~/

#On access2!!!
cat ~/unmapped.ids \
    | xargs -n100 \
    | sed 's/\s/,/g' \
    | xargs -L1 -I IDS curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=IDS&rettype=fasta&retmode=xml" > ~/unmapped.seqs.xml

#On computing node!!!
mv ~/unmapped.seqs.xml ./ && rm ~/unmapped.ids

cat unmapped.seqs.xml \
    | grep -A1 'TSeq_sid' \
    | sed 's/\s\+//; s/<\/\{0,1\}TSeq_.\{1,3\}id>//g' \
    | awk 'NR%3==1{printf("\n%s\t", $1)}NR%3==2{printf("%s",$1)}' \
    | sed 1d > unmapped.taxids

join -1 2 -2 1 -t$'\t' \
    <(cat ids | sort) \
    <({ \
        cat nr.*taxids \
            | sed 's/|/\t/g' \
            | cut -f 1-3; \
        cat unmapped.taxids \
        | awk '{print $1 "\t" $1 "\t" $2}'; \
        } \
        | sort -k1,1b) > acc2taxid.mapped

awk 'FILENAME=="acc2taxid.mapped" {a2t[$2]=$4; next} \
    {if(a2t[$2]!=""){a=a2t[$2]}else{a="NA"} print $0 "\t" a}' acc2taxid.mapped nr.blastp.out > nr.blastp.taxids.out

#677 contigs are primate-specific.
#$outd/trinotate/filter_primates/filter_BBH.pl $ncbi_tax/primates.taxid nr.blastp.taxid.out | less 
$outd/trinotate/filter_primates/filter_BBH.pl \
    $ncbi_tax/primates.taxid \
    nr.blastp.taxid.out \
    | cut -f1 \
    | sort \
    | uniq > primate.contig.ids

join -t$'\t' -1 1 -2 2 -o '2.1,2.2,2.3,2.4,2.5' \
    <(cat primate.contig.ids | sed 's/\.p[0-9]\+//' | sort) \
    <(cat $outd/abundance_estimation/salmon.isoform.TMM.EXPR.matrix.E-inputs | sort -k2,2) \
    | sort -k1n | less

#There is one highly-expressed contig with specificity to primates
# 30      CL748Contig1_1                  3727    1394.4  4811.7
# 86      TRINITY_DN5894_c0_g1_i1_1       1309    57.7    180.9
# 93      TRINITY_DN7392_c0_g1_i1_1       1097    30.6    78.8
# 94      TRINITY_DN14113_c0_g1_i1_1      1307    21.1    68.6
# 94      TRINITY_DN14917_c0_g1_i1_1      409     25.7    69.3
# 94      TRINITY_DN3034_c0_g1_i1_1       733     28.0    70.5
# 94      TRINITY_DN3478_c0_g1_i2_1       424     27.9    75.0
# 94      TRINITY_DN4439_c0_g1_i2_1       521     26.5    70.2
# 94      TRINITY_DN4877_c0_g1_i1_1       539     19.7    63.4
# 94      TRINITY_DN5792_c0_g1_i1_1       746     21.7    68.1
# 95      TRINITY_DN4207_c0_g1_i1_1       867     22.1    56.5
# 95      TRINITY_DN4257_c0_g2_i1_1       635     25.0    61.9
# 95      TRINITY_DN4465_c0_g1_i1_1       494     19.0    56.8
# 95      TRINITY_DN4894_c0_g1_i1_1       530     27.6    59.7

#The rest 816 are low-expressed as estimated by RSEM analysis. Nine contigs fall into expression 86-94 percentiles, whereas 808 are expressed below Ex95 level (see Supplementary tables).
cat $outd/trinotate/transcripts.fa | grep -A1 '>CL748Contig1_1' > CL748Contig1_1.fa

st=$HOME/tools/compare/samtools-1.3.1/samtools
fasta=$outd/abundance_estimation/transcripts_fpkm_1.fa.RSEM.idx.fa
bam=$outd/abundance_estimation/RSEM_lag-0_rep1/bowtie2.bam
region=CL748Contig1_1

cp $bam ./
$st sort -@36 bowtie2.bam > bowtie2.sorted.bam
$st index bowtie2.sorted.bam
$st mpileup -f $fasta -r $region $bam | less
#This contig seems to be chimeric. Its low-covered part starting from 2550bp has significant hit with human/bos sequences, 
#while 5'-end with high coverage has no hits against NR.

#OK. Let's prepare list of all contigs with BBH to primates, rodents, and cows.
for l in primates bovidae muridae do; 
    $outd/trinotate/filter_primates/filter_BBH.pl \
        $ncbi_tax/$l.taxid \
        nr.blastp.taxid.out \
        | cut -f1 \
        | sort \
        | uniq > $l.contig.ids
done
cat *contig.ids | grep -v CL748Contig1_1 > contaminant.contig.ids

#For contifg CL748Contig1_1, the 3'-end should be trimmed starting from 2200bp.

#816 contigs from contaminant set were removed for downstream analyses.

cp nr.blastp.*out nr.*taxids unmapped.* acc2taxid.mapped *contig.ids $outd/trinotate/filter_primates/