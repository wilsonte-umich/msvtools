#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %fieldNames));

#three mate-pair samples/runs for mouse strain 129P2
#ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR026/ERR026749/ERR026749_1.fastq.gz
#ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR026/ERR026750/ERR026750_1.fastq.gz
#ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR026/ERR026751/ERR026751_2.fastq.gz

sub get129 {
    #for unclear reasons, the ENA mouse long insert runs have precious few reads!
    #combine all 129s
    #my @runs = ('ERR026749', 'ERR026750', 'ERR026751'); #129P2
    my @runs = (#'ERR026740', 'ERR026744', 'ERR026745',  #129S1/SvImJ
                 'ERR026745', 'ERR026741', 'ERR026742', 'ERR026743'); #129S5
    foreach my $run(@runs){
        getDirectories($run, \my%dirs);
        -d $dirs{sample} or mkdir $dirs{sample};
        foreach my $readN(1..2){
            my $read = "read$readN";
            -d $dirs{$read} or mkdir $dirs{$read};
            my $outFile = "$dirs{$read}/$run\_$readN.fq.gz";
            my $command = "wget -q ";
            $command .= "-O $outFile ";
            $command .= "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR026/$run/$run\_$readN.fastq.gz";
            system($command);
            system("gunzip $outFile");  
        }   
    }
}

1;

