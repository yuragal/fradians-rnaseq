#!/usr/bin/perl

use Data::Dumper;

while(<>){
    chomp;
    @l=split("\t",$_);
    #print $_."\n";
    #print Dumper(\@l);
    @blgo = grep {/GO/} @l[12,16];
    #print 'blgo:',Dumper(\@blgo);
    $bg='.';
    $bg=join("`",@blgo) unless $#blgo<0;
    print join("\t",@l[0..11],$bg,@l[13..15])."\n"
}
