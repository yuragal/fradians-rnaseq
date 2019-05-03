#!/usr/bin/perl

use Array::Utils qw(:all);
use Data::Dumper;

while(<>){
    chomp;
    @l=split("\t");
    (%blgo, %pfgo, %tmgo)=();
    map {s/^.$//} @l;

    %blgo = map {@a=split('\^'); {$a[0] => [ @a[1,2] ]} } split('`',$l[1]) unless $l[1] eq '';
    %pfgo = map {@a=split('\^'); {$a[0] => [ @a[1,2] ]} } split('`',$l[2]) unless $l[2] eq '';
    %tmgo = map {@a=split('\^'); {$a[0] => [ @a[1,2] ]} } split('`',$l[3]) unless $l[3] eq '';
    @tm=( keys %tmgo );
    @trino=( keys %blgo, keys %pfgo );
    @newgo=array_minus(@tm,@trino);
    $count_go += scalar @newgo;
    $count_genes ++ if scalar @newgo > 0;
    @out= map { join('^',$_, @{ $tmgo{$_} }) } @newgo;
    #print 'blGO:', Dumper(\%blgo);
    #print 'pfGO:', Dumper(\%pfgo);
    #print 'tmGO:', Dumper(\%tmgo);
    #print 'newGO:', Dumper(\@newgo);
    print join("\t",@l[0..2],join('`',@out))."\n";
}

print STDERR "Number of genes with updated annotation: $count_genes\n";
print STDERR "Number of new GOs: $count_go\n";
