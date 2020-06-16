use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs $paramCat));

sub mapReadsUnique_Bowtie {
    my ($faFile) = @_;
    (-e $faFile) or die "$faFile does not exist";    
    status("mapping unique reads using Bowtie...");
    if ($param{maxDisc} > 3){
        status("Bowtie allows up to 3 mismatches, changing maxDisc to 3...\n");
        $param{maxDisc} = 3;
    } 
    my $outputFile = "$faFile.bowtie";  
    unlink ($outputFile);  #clear any previous mapping output
    my $bowtieCommand = "$param{bowtiePath}/bowtie "; 
    $bowtieCommand .= "-f -B 1 ";  #fasta, report refSeq 1-referenced
    $bowtieCommand .= " -a -m 1 --best --strata";    
    $bowtieCommand .= " -v $param{maxDisc} ";
    $bowtieCommand .= " $param{refSeq} ";
    $bowtieCommand .= " $faFile ";
    $bowtieCommand .= " $outputFile "; 
    $bowtieCommand .= " > $outputFile.stdout";
    runBowtie($bowtieCommand);        
    status("mapping done");
}

1;

