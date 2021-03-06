##############################################################################################
# Add annotations from TRAPID and Mercator
##############################################################################################

trinotate_dir=$HOME/tools/annotation/Trinotate-v3.1.1
outd=$HOME/axSA_assemblies/drap/stranded/trinity_norm_altogether
tempd=/store/$USER

mkdir -p $outd/trinotate/more_annotations && cd $outd/trinotate/more_annotations

#There are TRAPID and Mercator annotations computed for CDS sequences generated by Transdecoder 
# ../pfam_blastp_hits.transcripts.fa.transdecoder.cds
#These files were downloaded from TRAPID and Mercator web-sites:
# mercator.go.txt
#   classification of GO terms in Mercator -- BIN2GO.txt
# TRAPID.go.txt

#What I am going to do here is to push the missing GO terms into Trinotate results

#First, prepare two-column files with trID->GO mapping

#Mercator 2779 GOs
join -t$'\t' -o'2.2,1.2' \
    <(cat BIN2GO.txt \
        | sed 1d \
        | cut -f1,3 \
        | sed "s/'//g;s/go/GO/" \
        | awk -F '\t' '$2!=""' \
        | sort -k1,1b) \
    <(cat mercator.results.txt \
        | cut -f1,3 \
        | sed "1d;s/'//g" \
        | awk -F '\t' '$2!=""' \
        | sort -k1,1b) \
    | sed 's/trinity_dn/TRINITY_DN/;s/^cl\([0-9]\+\)contig/CL\1Contig/' \
    | sort k1,1b > mercator.go.txt

#TRAPID 31248 GOs
cat TRAPID.go.results.txt \
    | sed 1d \
    | awk '$6==0{print $2 "\t" $3}' \
    | sort -k1,1 > TRAPID.go.txt

#Then, collapse GOs by geneID
perl -F'\t' -alne 'next if /^$/;$h{$F[0]} = join(",",$F[1],$h{$F[0]}); END{print $_."\t".$h{$_} for sort keys %h}' TRAPID.go.txt  > TRAPID.go.collapsed
perl -F'\t' -alne 'next if /^$/;$h{$F[0]} = join(",",$F[1],$h{$F[0]}); END{print $_."\t".$h{$_} for sort keys %h}' mercator.go.txt  > mercator.go.collapsed

#Number of genes with GOs:
#   2444 mercator.go.collapsed
#   7426 TRAPID.go.collapsed

join -a2 -e NOT_IN_TRINOTATE -t$'\t' -o '0,2.2,1.2' \
    <(cat ../trinotate.go_annotations.txt | sort -k1,1b) \
    <(cat mercator.go.collapsed | sed 's/\.p[0-9]\+//;s/,$//' | sort -k1,1b) \
    | grep TRINOTATE | wc -l
#out of 2444, 30 proteins annotated by Mercator have not been annotated by Trinotate

join -a2 -e NOT_IN_TRINOTATE -t$'\t' -o '0,2.2,1.2' \
    <(cat ../trinotate.go_annotations.txt | sort -k1,1b) \
    <(cat TRAPID.go.collapsed | sed 's/\.p[0-9]\+//;s/,$//' | sort -k1,1b) \
    | grep TRINOTATE | wc -l
#TRAPID annotated 336 more additional genes
#5 new genes shared between TRAPID and mercator results
#there are 336+25=361 newly annotated genes in total

#Merge TRAPID and Mercator results
join -a1 -a2 -t$'\t' \
    <(cat TRAPID.go.collapsed | sort -k1,1b) \
    <(cat mercator.go.collapsed | sort -k1,1b) \
    | perl -F'\t' -MList::MoreUtils -alne \
        '@a=split(",",$F[1].$F[2]);@a=List::MoreUtils::uniq(@a); print $F[0],"\t",join(",",@a)' > both.go.collapsed

#Map TRAPID and Mercator GOs to new OBO
#Turned out, there were obsolete/alternative GO ids in this set

perl match_GO.pl $trinotate_dir/go-basic.obo both.go.collapsed 1> both.go.collapsed.matched 2>errlog
#errlog is empty; it seems that all deprecated GOs were successully mapped

#Annotations partly overlap; seems ok
#Alex M. suggests to measure semantic similarity of annotations obtained with various tools.
#Good idea. However, I'm not sure it could be properly implemented without considerable R&D effort.
#Maybe sometimes later I will come back here...
join -a1 -a2 -t$'\t' \
    <(cat ../trinotate.annotation_report.csv \
        | awk -F'\t' '$5!="."' \
        | cut -f5,13,14 \
        | sort -k1,1b) \
    <(cat both.go.collapsed.matched \
        | sort -k1,1b) \
    | sed 's/\t/\n\t/g'

#Now merge Trinotate GOs with TRAPID/Mercator
join -a1 -a2 -t$'\t' -e . -o '0,1.2,1.3,2.2' \
    <(cat ../trinotate.annotation_report.csv \
        | awk -F'\t' '$5!="."' \
        | cut -f5,13,14 \
        | sort -k1,1b) \
    <(cat both.go.collapsed.matched \
        | sort -k1,1b) \
    | perl merge_GO.pl > merged.go.txt

# Number of genes with updated annotation: 6314
# Number of new GOs: 18622
# merged.go.txt -- is 4-column file with geneID, blastGO, pfamGO and TRAPID/Mercator GOs not found in blast/pfamGO annotations

#GenerateTrinotate report file with new GOs from TRAPID/Mercator added to pfamGO field.

o=$(seq -f'1.%g' 2 16 | xargs printf '%s,')
o="0,${o}2.2"
#Check if script joins GO-field correctly
join -a1 -a2 -t$'\t' -e . -o${o} <(cat ../trinotate.annotation_report.csv | sort -k1,1b ) <(cat merged.go.txt | cut -f1,4 | perl -ne 'chomp;split("\t");$g=$_[0]; $g=~s/\.p\d+//; print join("\t",$g,$_[1])."\n"' | sort -k1,1b) | perl put_GO_to_trinotate_report.pl | sed 's/\t/\n\t/g' | awk -vnf=17 '{if(NR%nf==0){n=nf}else{n=NR%nf} print n " " $0}' | less

#Write
join -a1 -a2 -t$'\t' -e . -o${o} \
    <(cat ../trinotate.annotation_report.csv \
        | sed 1d \
        | sort -k1,1b ) \
    <(cat merged.go.txt \
        | cut -f1,4 \
        | perl -ne 'chomp;split("\t");$g=$_[0]; $g=~s/\.p\d+//; print join("\t",$g,$_[1])."\n"' \
        | grep -v prot_id \
        | sort -k1,1b) \
    | (head -1 ../trinotate.annotation_report.csv && perl put_GO_to_trinotate_report.pl) > trinotate.annotation_report.csv

#Check if all is correct
cat ../trinotate.annotation_report.csv | sed 1d | cut -d$'\t' -f13 | grep -v '^\.$' | sed 's/`/\n/g' | wc -l
# 119145
cat trinotate.annotation_report.csv | sed 1d | cut -d$'\t' -f13 | grep -v '^\.$' | sed 's/`/\n/g' | wc -l
# 137767
# 137767−119145=18622

$trinotate_dir/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls trinotate.annotation_report.csv -G --include_ancestral_terms > trinotate.go_annotations.txt


perl match_GO.pl $trinotate_dir/go-basic.obo TRAPID.go.collapsed > TRAPID.go.collapsed.matched
perl match_GO.pl $trinotate_dir/go-basic.obo mercator.go.collapsed > mercator.go.collapsed.matched

join -a2 -a1 -t$'\t' -e ',' -o'0,1.2,1.3,2.2' \
    <(join -a2 -a1 -t$'\t' -e ',' -o'0,1.2,2.2' \
        <(cat merged.go.txt \
            | awk -F$'\t' '$4{print $1 "\t" $4}' \
            | sed 's/\^[^`]\+`\?/,/g' \
            | sort -k1,1b) \
        <(cat mercator.go.collapsed.matched \
            | sed 's/\^[^`]\+`\?/,/g' \
            | sort -k1,1b) \
    ) \
    <(cat TRAPID.go.collapsed.matched \
        | sed 's/\^[^`]\+`\?/,/g' \
        | sort -k1,1b) > compare_merged_with_TRAPID_Mercator.txt
#Columns:
# geneId
# merged annotations
# mercator annotations
# TRAPID annotations

for i in `seq 2 4`; do 
    ./compare_merged_with_TRAPID_Mercator.pl < compare_merged_with_TRAPID_Mercator.txt \
        | cut -f${i} \
        |  sed 's/,/\n/g' \
        | grep -v '^$' \
        | wc -l
done
# 1258  Mercator GOs
# 17162 TRAPID GOs
# 202   both GOs

