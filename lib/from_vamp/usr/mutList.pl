#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs));

sub mutList{
    my %genome;
    foreach my $chrom(1..nChrom()){
        my $file = "$param{refSeqPath}/$param{refSeq}/$reverseRefSeqs{$param{refSeq}}{$chrom}.seq";
        open my $fileH, "<", $file or die "$file\n";
        my $pos = 1;
        while(read($fileH, my $base, 1)){
            $genome{$chrom}{$pos} = $base;
            $pos++;
        }
        close $fileH;
    }
    open my $fileH, ">", "Mutations.csv";
    print $fileH join(",", qw(Chromosome Position 
                              Reference_Base Mutant_Base Extra_Base 
                              Wild_Type_Pool_Count Wild_Type_Reads Mutant_Pool_Count Mutant_Reads
                              LOD Functional))."\n";
    runSQL("select chromosome, position, discrepancytype, extra, s1_consistentcount, s1_readcount, s2_consistentcount, s2_readcount, consistentlod, consequence
from discs_vac6wt_vac6mut d
where s1_readcount >= 3
  and s2_readcount >= 3
  and ((discrepancytype >= 1 and discrepancytype <= 4
  and (select min(position)
       from discs_vac6wt_vac6mut
       where s1_readcount >= 3
         and s2_readcount >= 3
         and discrepancytype > 5
         and ((s1_consistentcount + s2_consistentcount)/ (s1_readcount + s2_readcount)) >= 0.2
         and chromosome = d.chromosome
         and abs(position - d.position) < 10) is null)
    or discrepancytype > 5)
  and not (chromosome = 7 and position = 924485)
  and not (chromosome = 7 and position = 924486)
  and not (chromosome = 7 and position = 924487)
  and ((s1_consistentcount + s2_consistentcount)/ (s1_readcount + s2_readcount)) >= 0.35
order by chromosome, position",
\my($chrom, $pos, $discType, $extra, $s1cc, $s1rc, $s2cc, $s2rc, $lod, $consequence));
    while(fetchRow()){
        $genome{$chrom}{$pos} eq $types{RevDiscsFull}{$discType} and die "mutant and wild-type bases identical!\n";
        my $functional = '-';
        $consequence and $functional = 'Yes';
        print $fileH join(",", ($reverseRefSeqs{$param{refSeq}}{$chrom}, $pos, 
                                  $genome{$chrom}{$pos}, $types{RevDiscsFull}{$discType}, $types{RevDiscsFull}{$extra}, 
                                  $s1cc, $s1rc, $s2cc, $s2rc,
                                  sprintf("%.1f", $lod), $functional))."\n";
    }
    close $fileH;
}

1;

