#!/usr/bin/perl

use warnings;
use strict;

while (<>) {
    chomp;
    my @parts = split /\t/;
    print "@","$parts[0]:$parts[2]:$parts[3]:$parts[4]:$parts[5]:$parts[10]#$parts[6]/$parts[7]\n";
    $parts[8] =~ tr/./N/;
    print "$parts[8]\n";
    print "+\n";
    $parts[9] = join('', map{chr(ord($_)-33)} split(//,$parts[9]));
    print "$parts[9]\n";
}

