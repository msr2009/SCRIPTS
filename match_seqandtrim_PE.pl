#!/usr/bin/perl

use strict;
use warnings;
#use Config::Tiny;
use Cwd;
use Getopt::Std;
use File::Copy;
use vars qw/ $opt_f $opt_o $opt_r $opt_s/;
my ($file, $rev, $outfile, $seq);
getopts('f:r:o:s:');

var_check();
my %acc;
open (OUT, ">junk.tmp") or die "WTF\n";                
open (IN, $file) or die;

my $count = 0;
while(my $line = <IN>){
    chomp($line);
    print OUT $line;
    $count = ($count + 1) % 4;
    
    if($count == 0){
        print OUT "\n"; 
    }else{
        print OUT "\t";
    }
}

close IN;
close OUT;
open (IN, "junk.tmp") or die "\n\nUsage: trim_5prime.pl <infile.fq> <outfile.fq> \n\n";
while(my $change = <IN>){
    chomp($change);
    if($change =~ m/$seq/ ){
        my @vals = split(/\t/, $change);
        my @id = split(/\#/, $vals[0]);
        $acc{$id[0]} = $change;
    }else{
        next;
    }
}
close IN;
system("rm junk.tmp");
open (OUT, ">junk.tmp") or die; 
open (IN, $rev) or die "\n\nUsage: trim_5prime.pl <infile.fq> <outfile.fq> \n\n";
$count = 0;
while(my $rev = <IN>){
    chomp($rev);
    print OUT $rev;
    $count = ($count + 1) % 4;
    
    if($count == 0){
        print OUT "\n"; 
    }else{
        print OUT "\t";
    }
}
close IN;
close OUT;
open (OUT, ">$outfile") or die;
open (IN, "junk.tmp") or die "\n\nUsage: trim_5prime.pl <infile.fq> <outfile.fq> \n\n";
while(my $changes = <IN>){
    chomp($changes);
    my @second = split(/\t/, $changes);
    my @idd = split(/\#/, $second[0]);
    if (exists($acc{$idd[0]})){
    print OUT "$second[0]\n";
    my $short = substr $second[1], 2, 19;
    print OUT "$short\n";   
    print OUT "$second[2]\n";
    my $short2 = substr $second[3], 2, 19;
    print OUT "$short2\n";
		my (@ID) = split /\t/, $acc{$idd[0]};
		print "$ID[0]\n";
                my $short3 = substr $ID[1], 15, 19;
		print "$short3\n";
		print "$ID[2]\n";
                my $short4 = substr $ID[3], 15, 19;
		print "$short4\n";
		next;
    }else{
        next;
    }
}         
system("rm junk.tmp");

#########################################################
# End Main body of Program                              #
#########################################################

#########################################################
# Start of Variable Check Subroutine "var_check"        #
# if we came from command line do we have information   #
# and is it the correct info.. if not go to the menu... #
# if we came from the web it should all check out...    #
#########################################################

sub var_check {

  if ($opt_f) {
   $file = $opt_f;
  } else {
   var_error();
  }

  if ($opt_o) {
   $outfile = $opt_o;
  } else {
   var_error();
  }
  if ($opt_s) {
   $seq = $opt_s;
  } else {
   var_error();
  }
  if ($opt_s) {
   $rev = $opt_r;
  } else {
   var_error();
  }

}

#########################################################
# End of Variable Check Subroutine "var_check"          #
#########################################################

#########################################################
# Start of Variable error Subroutine "var_error"        #
#########################################################
sub var_error {

  print "\n\n";
  print "  You did not provide all the correct command line arguments\n\n";
  print "  This script will match a specific sequence in the forward read, then create a reverse fastq file with the Paired end\n\nUsage:\n";
  print "  maergelines.pl -f forward.fq -r reverse_input.fq -o outfileReverseRead.fq  -s  sequence to match >outfileForwardRead.fq\n\n\n";
  print "\n\n\n";
  print "\n\n\n";
  exit 0;

}

#########################################################
# End of Variable error Subroutine "var_error"          #
#########################################################