#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %qseqColumns));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $command{qseq2fastq} or $command{qseq2fastq} = ['multiThread', '5:00:00', 1000, 0];        
                    
sub qseq2fastq {
    my ($sample) = @_;
    my ($inputSample, $outputSample) = splitSampleName($sample);
    getDirectories($inputSample, \my%inDirs);
    getReadFiles($outputSample, 'solexa_fastq', \my%fastqFiles);
    foreach my $readN (1..3) {
        my $read = "read$readN"; 
        my $qseq = "qseq$readN";   
        my $qseqDir = $inDirs{$qseq};
        -d $qseqDir or next;
        my $fastqFile = $fastqFiles{$read};
        if (-e $fastqFile){
            status("$fastqFile already exists, will not overwrite\n");
        } else {
            status("writing $fastqFile from qseq inputs\n");
            open my $fastqFileH, ">", $fastqFile or die "could not open $fastqFile: $!";
            foreach my $qseqFile (<$qseqDir/*qseq.txt>){ #read files from qseq subdir
                open my $qseqFileH, "<", $qseqFile or die "could not open $qseqFile: $!";
                while (<$qseqFileH>){
                    chomp $_;
                    my @in = split(/\t/, $_);                
                    if ($in[$qseqColumns{filter}]){                
#                        my $readID = join(":", ($in[$qseqColumns{machine}], 
#                                                @in[$qseqColumns{lane}..$qseqColumns{y}]) );  
                        #$readID .= "\#$in[$qseqColumns{index}]/$in[$qseqColumns{read}]";
                        my $readID = join(":", (@in[$qseqColumns{machine}..$qseqColumns{index}]) );  
                        print $fastqFileH '@'."$readID\n";
                        print $fastqFileH "$in[$qseqColumns{sequence}]\n";
                        print $fastqFileH '+'."$readID\n";
                        print $fastqFileH "$in[$qseqColumns{quality}]\n";     
                    }
                }
                close $qseqFileH;
            }
            close $fastqFileH;
        }
    }
}

#- Machine name: unique identifier of the sequencer.
#- Run number: unique number to identify the run on the sequencer.
#- Lane number: positive integer (currently 1-8).
#- Tile number: positive integer.
#- X: x coordinate of the spot. Integer (can be negative).
#- Y: y coordinate of the spot. Integer (can be negative).
#- Index: positive integer. No indexing should have a value of 1.
#- Read Number: 1 for single reads; 1 or 2 for paired ends.
#- Sequence (BASES)
#- Quality: the calibrated quality string. (QUALITIES)
#- Filter: Did the read pass filtering? 0 - No, 1 - Yes


1;



