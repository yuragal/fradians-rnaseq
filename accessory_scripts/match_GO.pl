#!/usr/bin/perl

use List::MoreUtils;
use Data::Dumper;
use Switch;

$obo=$ARGV[0];
$mygos=$ARGV[1];

#Load GO obo file
open OBO,'<',$obo;
%go;
while(<OBO>){
    #last if scalar keys %go > 100; #for test
    if(/\[Term\]/){
        while(<OBO>){
            last if /^$/;
            chomp;
            /^([^:]+):\s(.+)$/;
            switch($1){
                case 'id'           { $id   =$2         }
                case 'name'         { $name =$2         }
                case 'namespace'    { $ns   =$2         }
                case 'alt_id'       { push @alt_id,$2   }
            }
        }

        $go{$id}= { 
            name        => $name,
            namespace   => $ns,
            alt_id      => [ @alt_id ]
        };
        @alt_id=();
    }
}

#Map alt_ids to ids
%alt;
for $id (keys %go){
    if($#{ $go{$id}{'alt_id'} } > -1){
        @a=@{ $go{$id}{'alt_id'} };
        for $aid(@a){
            $alt{$aid}=$id;
        }
    }
}

#print 'go:',Dumper(\%go); #die;
#print 'alt:',Dumper(\%alt);

#Match GOs, check if all GOs present in $go or %alt hashes.
open MYGOS,'<',$mygos;
while(<MYGOS>){
    chomp;
    @l=split;
    @out=();
    for $i(split(',',$l[1])){
        if(defined($go{$i})){
            push @out, join('^',$i,$go{$i}{'namespace'},$go{$i}{'name'});
        }elsif(defined($alt{$i})){
            $ti=$alt{$i};
            push @out, join('^',$ti,$go{$ti}{'namespace'},$go{$ti}{'name'});
        }else{
            print STDERR "Error: No term $i found in OBO!!!\n";
        }
    }
    print $l[0]."\t".join('`',@out)."\n";
}