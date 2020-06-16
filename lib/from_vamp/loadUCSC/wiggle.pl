#!/usr/bin/perl -w
use strict;
use warnings;

use vars(qw(%param %command %types %fields %fieldNames %refSeqs %reverseRefSeqs));

#callable parameters and commands in this script
defined $command{loadWiggle} or $command{loadWiggle} = ['singleThread', '1:00:00', 2000, 0];

# fixedStep  chrom=chrN  start=position  step=stepInterval  [span=windowSize]
  # dataValue1
  # dataValue2

sub loadWiggle {
    my ($name, $url) =  @_;
    # #because UCSC does not use a consistent path structure, user must supply full http or ftp path
    status("loading wiggle file $url\n");
    
    #download the file
    $url =~ m/.*\/(.*)/;
    my $inFile = $1;
    system("wget $url -O $inFile");
    
    #unzip if needed
    $inFile =~ m/(.*).gz$/ and system("gunzip $inFile");
    $inFile = $1;
    
    #convert from wiggle to flat file
    open my $inFileH, "<", $inFile or die "could not open $inFile\n";
    my $outFile = "$name.csv";
    open my $outFileH, ">", $outFile or die "could not open $outFile\n";
    my ($chrom, $pos, $step);
    while (my $line = <$inFileH>){
        chomp $line;
        $line or next;
        if ($line =~ m/^fixedStep chrom=(.*) start=(.*) step=(.*)/){
            ($chrom, $pos, $step) = ($1, $2, $3);
            $chrom = $refSeqs{$param{refSeqBase}}{$chrom}; #convert chromosome to vamp numbered values
        } elsif ($chrom){
            print $outFileH "$chrom,$pos,$line\n";
            $pos += $step;
        }
    }
    close $inFileH;
    close $outFileH;
    unlink $inFile;
    
    #load flat file into database
    my $outTable = newTable('Wig', "$name\_$param{refSeq}");
    loadData($outFile, $outTable, ',', $fieldNames{Wig});
}



1;





