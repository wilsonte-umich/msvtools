#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $command{checkStrandedness} or $command{checkStrandedness} = ['multiThread', '1:00:00', 1000, 0];

my %geneNames = (hg18=>'TP53',hg19=>'TP53',mm9=>'Trp53');

sub checkStrandedness { 
    my ($sample) = @_;
    my $hitsTable = getTableName('Hits', $sample);
    my $genesTable = "refgene_all_$param{refSeq}";
    my $geneName = $geneNames{$param{refSeq}};
    runSQL("select bin, 
                   sum(case when strand = 0 then count_ when strand = 1 then count_ else 0 end) top,
                   sum(case when strand = 2 then -count_ else 0 end) bottom
            from (
            select round(h.position/1000)*1000 bin, h.strand, h.count_
            from $hitsTable h, 
            (select chromosome, corrstart_, end_ from $genesTable where name2 like '$geneName') g
            where h.chromosome = g.chromosome
              and h.position >= g.corrstart_ - 5000
              and h.position <= g.end_ + 5000
            )
            group by bin
            order by bin",
            \my($bin,$top,$bottom));
    my $outFile = "$param{logPath}/checkStrandedness_$sample.csv";
    open my $outFileH, ">", "$outFile" or die "could not open $outFile: $!\n";
    print $outFileH "bin,top,bottom\n";
    while (fetchRow()){ print $outFileH "$bin,$top,$bottom\n" }
    close $outFileH;
    system("Rscript $param{vampPath}/bin/ipSeq/checkStrandedness.R $outFile $sample");
}

1;

