#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs));

my %discTypes = (1=>'A',2=>'C',3=>'G',4=>'T');

sub countMuts{
    my (%discs, %bases, %muts);
    runSQL("select chromosome, position, discrepancytype from (
            select chromosome, position, discrepancytype, ((s1_consistentcount + s2_consistentcount)/ (s1_readcount + s2_readcount)) F
            from discs_vac6wt_vac6mut
            where s1_readcount >= 3
              and s2_readcount >= 3
              and discrepancytype >= 1 and discrepancytype <= 4 )
            where F >= 0.9 or (F >= 0.35 and F <= 0.65)", \my($chrom, $pos, $disc));
    while(fetchRow()){$discs{$chrom}{$pos} = $discTypes{$disc}}
    foreach my $chrom(1..nChrom()){
        my $file = "$param{refSeqPath}/$param{refSeq}/$reverseRefSeqs{$param{refSeq}}{$chrom}.seq";
        open my $inH, "<", $file;
        my $pos = 1;
        while(read ($inH, my $base, 1)){
            $discs{$chrom}{$pos} and $muts{$base}{$discs{$chrom}{$pos}}++;
            $bases{$base}++; 
            $bases{total}++; 
            $pos++;
        }
        close $inH;
    }
    
    
    foreach my $base (sort {$a cmp $b} keys %bases){print "$base\t$bases{$base}\n"}
    print "\nwtBase\tA\tC\tG\tT\n";
    foreach my $wtBase('A','C','G','T'){
        foreach my $mutBase('A','C','G','T'){$muts{$wtBase}{$mutBase} or $muts{$wtBase}{$mutBase} = -1}
        print "$wtBase\t$muts{$wtBase}{A}\t$muts{$wtBase}{C}\t$muts{$wtBase}{G}\t$muts{$wtBase}{T}\n";
    }
}

1;




