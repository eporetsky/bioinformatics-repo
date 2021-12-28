#!/usr/bin/perl -w
# superfamily.pl
# http://supfam.org
# Copyright (c) 2001 MRC and Julian Gough; see http://supfam.org/SUPERFAMILY/license.html
# David Morais 28.09.10

use strict;

my $usage="USAGE:\nsuperfamily.pl <genome.fa>\n\noptionally\nnohup superfamily.pl <genome.fa> & \n";

unless (@ARGV==1)
{
print $usage;	
	
}   
$ARGV[0]=~/^(\w+)\.fa$/;
#my $file=$1;

print "Running fasta checker\n";
system "./fasta_checker.pl $ARGV[0] > scratch/$ARGV[0].torun.fa";


print "Running hmmscan\n";
system "./hmmscan.pl -o scratch/$ARGV[0].res -E 10 -Z 15438 hmmlib scratch/$ARGV[0].torun.fa --hmmscan hmmscan --threads 32 --tempdir scratch ";

print "Running assingments\n";
system "./ass3.pl -t n -f 32 -e 0.01 -r scop-cla-latest.txt scratch/$ARGV[0].torun.fa scratch/$ARGV[0].res $ARGV[0].ass ";

print "Running ass_to_html\n";
system "./ass_to_html.pl scop-des-latest.txt model.tab $ARGV[0].ass > $ARGV[0].html";

