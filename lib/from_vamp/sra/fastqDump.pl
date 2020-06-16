#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command));

#callable parameters and commands in this script
defined $param{sratoolsPath} or $param{sratoolsPath} = 0;  
defined $param{maxSpotId} or $param{maxSpotId} = 0;  
defined $command{fastqDump} or $command{fastqDump} = ['multiThread', '12:00:00', 2000, 0];

sub fastqDump { 
    my ($sample) = @_;
    $param{sratoolsPath} or die "must specify parameter sratoolsPath\n";
    getDirectories($sample, \my%dirs);
    getReadFiles($sample, 'fasta', \my%readFiles);
    my $sraFile = "$dirs{sample}/$sample.sra";
    -e $sraFile or die "could not find $sraFile\n";
    -d $dirs{read1} or mkdir $dirs{read1};
    -d $dirs{read2} or mkdir $dirs{read2};
    my $maxSpotID = "";
    $param{maxSpotId} and $maxSpotID = "--maxSpotId $param{maxSpotId}";
    system("$param{sratoolsPath}/fastq-dump --origfmt --fasta $maxSpotID -A $sample $sraFile"); #puts files into working directory
    system("mv $sample\_1.fasta $readFiles{read1}");
    system("mv $sample\_2.fasta $readFiles{read2}"); 

    
    # open my $inH, "<", $readFiles{read1} or die "could not open $readFiles{read1}: $!\n";
    # my %lanes;
    # while (my $line = <$inH>){
        # $line =~ m/(.*)\:(.*)\:\d*\:\d*\:\d*/ or next;
        # my $machine = $1;
        # my $lane = $2;
        # $lanes{$machine}{$lane}++;   
    # }
    # close $inH;
    # foreach my $machine(sort {$a cmp $b} keys %lanes){
        # print "$machine\n";
        # foreach my $lane(sort {$a cmp $b} keys %{$lanes{$machine}}){
            # print "  $lane:\t$lanes{$machine}{$lane}\n";
        # }  
    # }
     
}

1;
