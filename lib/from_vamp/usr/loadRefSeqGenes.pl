#!/usr/bin/perl
use strict;
use warnings;

our (%param, %types, %fields, %refSeqs);

my %inCols = (
    bin=>0 ,
    name=>1  ,
    chrom=>2  ,	
    strand=>3  ,	
    txStart=>4  ,	
    txEnd=>5  ,	#i.e. chrom
    cdsStart=>6  ,	
    cdsEnd=>7  ,
    exonCount=>8  ,	
    exonStarts=>9  ,
    exonEnds=>10  ,	
    id=>11 ,	
    name2=>12 , 	
    cdsStartStat=>13,  
    cdsEndStat=>14,
    exonFrames=>15,
);
my %strands = ('+' => 1, '-' => 2);
my ($id, $cdsFileH, $ncdsFileH);

sub parseXCDS{
    my ($inputFile) = @_;;

    open my $inFileH, "<", $inputFile or die "could not open $inputFile";
    my $cdsTable = newTable('CDS', $param{refSeq});
    my $cdsFile = "$cdsTable.csv";
    open $cdsFileH, ">" , $cdsFile;
    my $ncdsTable = newTable('NCDS', $param{refSeq});
    my $ncdsFile = "$ncdsTable.csv";
    open $ncdsFileH, ">" , $ncdsFile;

    my %ignored;
    while (<$inFileH>){
        $_ =~ m/^#/ and next;     
        chomp $_;
        my @inFields = split("\t", $_);
        defined $refSeqs{$param{refSeq}}{$inFields[$inCols{chrom}]} or (${$ignored{chrom}}{$inCols{chrom}}++ and next);
        defined $strands{$inFields[$inCols{strand}]} or (${$ignored{strand}}{$inCols{strand}}++ and next);
        my $name1 = $inFields[$inCols{name}];  
        $name1 or $name1 = '';
        my $name2 = $inFields[$inCols{name2}];  
        $name2 or $name2 = '';
        my @common = ($id,
                      $name1,
                      $name2,
                      $refSeqs{$param{refSeq}}{$inFields[$inCols{chrom}]},      
                      $strands{$inFields[$inCols{strand}]});    
        my @exonStarts = split(",", $inFields[$inCols{exonStarts}]);
        my @exonEnds = split(",", $inFields[$inCols{exonEnds}]);
        my @exonFrames = split(",", $inFields[$inCols{exonFrames}]);
        my $cdsStart = $inFields[$inCols{cdsStart}];
        my $cdsEnd = $inFields[$inCols{cdsEnd}];
        my $isCoding = 0;
        foreach my $frame(@exonFrames){$frame > -1 and $isCoding = 1}
        my $carryOver = 0;    
        foreach my $exon(0..($inFields[$inCols{exonCount}] - 1)){  
            if ($isCoding){
                if($exonFrames[$exon] == -1){ #wholly non-coding exon
                    commitUTR($exon, $exonStarts[$exon] + 1, $exonEnds[$exon], @common);
                    next; 
                }
                unless($exonStarts[$exon] >= $cdsStart){ #first encountered coding exon
                    commitUTR($exon, $exonStarts[$exon] + 1, $cdsStart - 1 + 1, @common);
                    $exonStarts[$exon] = $cdsStart;                         
                } 
                unless ($exonEnds[$exon] <= $cdsEnd){ #last encountered coding exon
                    commitUTR($exon, $cdsEnd + 1, $exonEnds[$exon], @common);
                    $exonEnds[$exon] = $cdsEnd;
                }
                $id++; 
                my $line = join(",", @common, $exon + 1, $exonStarts[$exon] + 1, $carryOver, $exonEnds[$exon]); 
                print $cdsFileH "$line\n";
                $carryOver = ($exonEnds[$exon] - $exonStarts[$exon]) % 3;       
            } else {
                $id++; 
                my $line = join(",", @common, $exon + 1, $exonStarts[$exon] + 1, $exonEnds[$exon]); 
                print $ncdsFileH "$line\n";
            }   
        }
    }
    close $inFileH;
    close $cdsFileH;
    loadData($cdsFile, $cdsTable, ",",    "CDSID, NAME1, NAME2, CHROMOSOME, STRAND, EXON, START_, CARRYOVER, END_");
    loadData($ncdsFile, $ncdsTable, ",", "NCDSID, NAME1, NAME2, CHROMOSOME, STRAND, EXON, START_, END_");
    
    print "ignored chroms:\n";
    foreach my $ignored(keys %{$ignored{chrom}}){print "  $ignored\n"}          
    print "ignored strands:\n";
    foreach my $ignored(keys %{$ignored{strand}}){print "  $ignored\n"} 
}


sub commitUTR{
    my ($exon, $exonStart, $exonEnd, @common) = @_;
    $id++; 
    $common[2] .= "_UTR";
    my $line = join(",", @common, $exon + 1, $exonStart, $exonEnd); 
    print $ncdsFileH "$line\n";
}

1;
