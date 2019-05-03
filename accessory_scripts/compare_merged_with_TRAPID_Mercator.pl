#!/usr/bin/perl

use Data::Dumper;
use Array::Utils qw(:all);

while(<>){
    chomp;
    $l=$_;
    split("\t",$l);
    @a=split(",",$_[1]);
    @b=split(",",$_[2]);
    @c=split(",",$_[3]);
    
    @a=unique(@a);
    @b=unique(@b);
    @c=unique(@c);
    
    @m=array_minus(@b,@c);
    @m=intersect(@a,@m);
    
    @t=array_minus(@c,@b);
    @t=intersect(@a,@t);
    
    @mt=intersect(@b,@c);
    @mt=intersect(@a,@mt);
    
    
    #print 'A:',Dumper(\@a);
    #print 'B:',Dumper(\@b);
    #print 'C:',Dumper(\@c);
    #print "==$l\n";
    #print "Terms from Mercator:",Dumper(\@m)."\n" if scalar @m;
    #print "Terms from TRAPID:",Dumper(\@t)."\n" if scalar @t;
    #print "Terms from Both:",Dumper(\@mt)."\n" if scalar @mt;
    if(@a != (@m + @t + @mt)){
        print "ERROR when compairing GO-terms!!!\n";
        print 'Input:'.(scalar @a)."\tinMercator:".(scalar @b)."\tinTRAPID:".(scalar @c)."\n";
        print 'Input:'.(scalar @a)."\tMercator:".(scalar @m)."\tTRAPID:".(scalar @t)."\tBoth:".(scalar @mt)."\n";
    }
    print join("\t",
        $_[0],
        join(',',@m).',',
        join(',',@t).',',
        join(',',@mt).','
        )."\n";
    
}