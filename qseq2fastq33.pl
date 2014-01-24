#!/usr/bin/perl

use warnings;
use strict;

my $base;
my $qual;

while (<>) {
    chomp;
    my @parts = split /\t/;
    #print read ID
    print "@","$parts[0]:$parts[2]:$parts[3]:$parts[4]:$parts[5]:$parts[10]#$parts[6]/$parts[7]\n";
    
    #replace .'s with N's
    my $newread = "";
    foreach $base (split //, $parts[8]) {
    	if( $base =~ /\./ ) { $newread .= "N"; }	
    	else { $newread .= $base; }
    }
    print "$newread\n";
    
    # print 3rd line
    print "+\n";
    
    # print new phred score (sanger)
    my $newquality = "";
    foreach $qual (split //, $parts[9]) {
 		$newquality .= chr(ord($qual)-31);
    }
    print "$newquality\n";
}

# export lane=8; e=1; cat s_${lane}_${e}_????_qseq.txt | perl /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/qseq2fastq33.pl > s_${lane}_${e}.fq
