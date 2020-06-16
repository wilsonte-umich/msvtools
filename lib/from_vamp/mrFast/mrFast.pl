#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %gffFields %bowtieFields));

#callable parameters and commands in this script
defined $param{mrFastPath} or $param{mrFastPath} = "$param{vampPath}/mrFast"; 
defined $command{mapReadsMrFast} or $command{mapReadsMrFast} = ['multiThread', '24:00:00', 10000, 0];

sub mapReadsMrFast {
    my ($sample) = @_;
    my $fastaOutFile = getMrFastIndex();
    getReadFiles($sample, $param{readType}, \my%readFiles);
    getMapFiles($sample, \my%mapFiles);
    

    my $mrFastCommand = "$param{mrFastPath}/mrfast --pe "; #Search will be done in paired-end mode
    $mrFastCommand .= "--search $fastaOutFile "; #Search the specified genome. Index file should be in same directory as the fasta file.
    $mrFastCommand .= "-o $mapFiles{sample} "; #Output of the mapped sequences (SAM format).
    $mrFastCommand .= "-e $param{maxDisc} "; #Maximum allowed edit distance (default 2).
    $mrFastCommand .= "--crop $param{readLength} "; #Crop the input reads at position [int].
    $mrFastCommand .= "--seq1 $readFiles{read1} --seq2 $readFiles{read2} ";
    system($mrFastCommand);
#--min [int]	Min inferred distance allowed between two pairend sequences.
#--max [int]	Max inferred distance allowed between two pairend sequences.
#--discordant-vh	To return all discordant map locations ready for the Variation Hunter program, and OEA map locations ready for the NovelSeq.


}



sub getMrFastIndex {
    my $windowSize = int( $param{readLength} / ($param{maxDisc} + 1) );
    my $fastaInFile = "$param{refSeqPath}/$param{refSeqBase}.fa";
    my $fastaOutFile = "$fastaInFile.$windowSize";
    my $indexfile = "$fastaOutFile.index";
    unless (-e $indexfile) {
        status("creating mrFast index for refSeq $param{refSeqBase} at window size $windowSize");
        system("cat $fastaInFile > $fastaOutFile");
        system("$param{mrFastPath}/mrfast --index $fastaOutFile --ws $windowSize");
        unlink($fastaOutFile);
        -e $indexfile or die "failed to create mrFast index file $indexfile\n";
    } 
    return $fastaOutFile;
}





1;
