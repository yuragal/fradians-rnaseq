#!/usr/bin/perl

#Load list of taxids to filter out
%p = map { $_ => 1 }  grep { chomp; } `cat $ARGV[0]`;

open BH,'<',$ARGV[1];
while(<BH>){
    chomp;
    $l=$_;
    @a=split("\t",$l);
    if($q eq $a[0]){
        push @bh, $l;
    }else{
        #Simple check for Hit id>90%,qcovs>95% and subj_taxid belongs to Primates
        if( $id[0]>=90 && $ev[0]<=1e-20 && $qcovs[0]>=40 && defined($p{ $taxid[0] }) ){
            if($ARGV[2]=~/bbh/i){
                print $bh[0]."\n";
            }else{
                print join("\n", @bh)."\n";
            }
        }
        @bh = ($l);
        $q=$a[0];
        (@s,@id,@al,@ev,@qcovs,@qcovh,@taxid)=();
    }
    push @s, $a[1];
    push @id, $a[2];
    push @al, $a[3];
    push @ev, $a[10];
    push @qcovs, $a[14];
    push @qcovh, $a[15];
    push @taxid, $a[16];
}