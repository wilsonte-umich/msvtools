#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %samFields %samFlags));

#callable parameters and commands in this script
#binSize is required but already defined by VAMP

defined $command{mapMouseNormal} or $command{mapMouseNormal} = ['singleThread', '48:00:00', 2000, 0];

my @vampSamples = ('iMur','C57BL6');
# my @sangerSamples = ('129P2','129S1','129S5','AKR','A_J',
                     # 'BALB','C3H','C57BL','CAST','CBA',
                     # 'DBA','FVB','LP_J','NOD','NZO',
                     # 'PWK','SPRET','WSB');
my @sangerSamples = ('129P2');
my ($minSanger, $maxSanger) = (100, 700);

sub mapMouseNormal {

#2006x2486	3	87210129	87219289	87212601	87215515	87215746	87216817	4216	231
    my ($chrom,$start,$end,$p1,$p2,$p3,$p4) = (3,87210129,87219289,87212601,87215515,87215746,87216817);
    my $center = int($p2 + (($p3 - $p2)/2));
    my $width = $end - $start;
    my $halfWidth = int($width/2);
    ($start, $end) = ($center - $halfWidth, $center + $halfWidth);

    #processMouseVamp($chrom, $start, $end, $center);
    processMouseSanger($chrom, $start, $end, $center);
}

sub getMouseInsertionPaths {
    my ($sample, $chrom, $start, $end, $dirs, $readFiles) = @_;
    getDirectories($sample, $dirs);
    -d $$dirs{sample} or mkdir $$dirs{sample};
    my $region = "$chrom\_$start\_$end";
    getDirectories("$sample/$region", $dirs);
    -d $$dirs{sample} or mkdir $$dirs{sample};
    -d $$dirs{read1} or mkdir $$dirs{read1};
    -d $$dirs{read2} or mkdir $$dirs{read2};
    $$readFiles{read1} = "$$dirs{read1}/$sample\_$region\_1.fa";
    $$readFiles{read2} = "$$dirs{read2}/$sample\_$region\_2.fa"; 
}
sub processMouseVamp {
    my ($chrom, $start, $end, $center) = @_;    
    foreach my $sample (@vampSamples){
        runSQL("SELECT position1 low, position2 + length2 - 1 high
                FROM frags_$sample
                WHERE fragmenttype = $types{Frags}{Normal}
                  AND chromosome1 = $chrom
                  AND position1 <= $end
                  AND $start <= position2",
                \my($low,$high));   
        my %counts;
        foreach my $i($start..$end){ $counts{$i} = 0 } 
        while(fetchRow()){ processMouseNormal($low, $high, $start, $end, \%counts) }
        processMouseSample($chrom, $start, $end, $center, \%counts, $sample);
    }   
}
sub processMouseSanger {
    my ($chrom, $start, $end, $center) = @_;    
    fixBamSexChrom(\$chrom); #Sanger bams use 1,2,3...X,Y
    foreach my $sample (@sangerSamples){
        getMouseInsertionPaths($sample, $chrom, $start, $end, \my %dirs, \my%readFiles);
        -e $readFiles{read1} and next;  
        my $bamPath = "ftp://ftp-mouse.sanger.ac.uk/current_bams/$sample.bam";
        samtoolsView($bamPath, [[$chrom, $start, $end]], \my@samResults);
        sam2FastaArray(\@samResults, ">", \%readFiles, \my%pairs, \my$pairID);
        getPartnerRegions(\%pairs, \my@regions, \my%names);
        my $partners = samtoolsFastaIndex(\@regions, 1);
        partnerRegions2Fasta($partners, \%names, ">>",  \%readFiles, \%pairs, \$pairID);

        # my %counts;
        # foreach my $i($start..$end){ $counts{$i} = 0 } 
        # foreach my $line(@samResults){
            # $line =~ m/^\@/ and next;
            # my @line = split("\t", $line);
            # my $flag = $line[$samFields{FLAG}];
            # getSamFlagBit($flag, 'isReverse') and next;
            # getSamFlagBit($flag, 'nextIsReverse') or next;
            # my $low = $line[$samFields{POS}];
            # my $high = $line[$samFields{PNEXT}] + length($line[$samFields{SEQ}]);
            # my $tLen = $high - $low;
            # ($tLen > $minSanger and $tLen < $maxSanger) or next;
            # processMouseNormal($low, $high, $start, $end, \%counts); 
        # }
        # processMouseSample($chrom, $start, $end, $center, \%counts, $sample);   
    }
}

sub processMouseNormal {
    my ($low, $high, $start, $end, $counts) = @_;
    my $lowI = $low;
    $lowI < $start and $lowI = $start;
    my $highI = $high;
    $highI > $end and $highI = $end;
    $highI > $lowI or die "processMouseNormal: bad position order!\n";
    foreach my $i($lowI..$highI){ $$counts{$i}++ }    
}
sub processMouseSample {
    my ($chrom, $start, $end, $center, $counts, $sample) = @_;
    my ($nPos, $sum) = ($end - $start + 1);
    foreach my $i($start..$end){ $sum += $$counts{$i} } 
    my $mean = $sum / $nPos;
    foreach my $i($start..$end){ $$counts{$i} /= $mean } 
    my $outFile = "$param{inputPath}/mapMouseNormal_$sample.csv";
    open my $outFileH, ">", $outFile or die "could not open $outFile: $!";  
    print $outFileH "pos,normCount\n";
    foreach my $i($start..$end){ 
        my $pos = $i - $center;
        print $outFileH "$pos,$$counts{$i}\n"; 
    } 
    close $outFileH;
    system("Rscript $param{vampPath}/bin/mouseInsertion/mouseInsertion.R $outFile"); 
}



1;
