#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#user.pl is a place for subs you wish to execute that are out of the main
#program content of flow of VAMP, e.g. single specific analysis queries, etc.
#These are executed using the 'run' command in an instructions file.
#########################################################################

use vars(qw(%param %types %fields %refSeqs));

sub simulateAllCNVs{
    my ($CNVsTable) = @_;
    my $nIterations = 3;
    my (%hsCovs);
    for my $iteration(1..$nIterations){
        print "$iteration \n";
        my $cnvSimTable = generateCNVSimulation($CNVsTable);
        mapAllCNVS($cnvSimTable, \%hsCovs);        
    }
    print "\n ";
    
    foreach my $hsCov(sort {$a <=> $b} keys %hsCovs){
        my $meanCount = $hsCovs{$hsCov}/$nIterations;
        print "$hsCov,$meanCount\n";
        
    }
}

sub generateCNVSimulation {
    my ($CNVsTable) = @_;
    my $nChrom = nChrom();
    $param{sex} and $param{sex} eq 'XX' and $nChrom--;
    my $chromInfoTable = "CHROMINFO_$param{refSeqBase}";
    runSQL("SELECT Sum(End_) FROM $chromInfoTable WHERE CHROMOSOME <= $nChrom", \my($genomeSize));
    fetchRow();
    runSQL("SELECT CHROMOSOME, END_ FROM $chromInfoTable WHERE CHROMOSOME <= $nChrom ORDER BY CHROMOSOME");
    my $chromsRef = fetchAllHashRef();
    my $maxI = scalar(@$chromsRef) - 1;
    my $cnvSimTable = "$CNVsTable\_SIM";
    dropTable($cnvSimTable);
    runSQL("CREATE TABLE $cnvSimTable (CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER)"); 
    my $cnvSimFile = "$cnvSimTable.csv";  
    open my $cnvSimFileH, ">", $cnvSimFile;    
    runSQL("SELECT END_ - START_ FROM $CNVsTable", \my($cnvSize));
    CNVSIM: while(fetchRow()){
        my $cnvIndex = int(rand $genomeSize);
        my ($cumGenomeSize, $cnvChrom, $cnvStart, $cnvEnd) = (0);
        foreach my $i(0..$maxI){
            if($cnvIndex <= ($cumGenomeSize + $$chromsRef[$i]{END_})){        
                $cnvChrom = $$chromsRef[$i]{CHROMOSOME};          
                $cnvStart = $cnvIndex - $cumGenomeSize;        
                $cnvEnd =  $cnvStart + $cnvSize;
                if($cnvEnd <= $$chromsRef[$i]{END_}){              
                    commitCNVSim($cnvSimFileH, $cnvChrom, $cnvStart, $cnvEnd);
                } else {
                    commitCNVSim($cnvSimFileH, $cnvChrom, $cnvStart, $$chromsRef[$i]{END_});
                    my $j = $i + 1;
                    $j > $maxI and $j = 0;
                    $cnvChrom = $$chromsRef[$j]{CHROMOSOME}; 
                    $cnvEnd = $cnvSize - ($$chromsRef[$i]{END_} - $cnvStart);                    
                    commitCNVSim($cnvSimFileH, $cnvChrom, 1, $cnvEnd);                    
                }
                next CNVSIM;
            }
            $cumGenomeSize += $$chromsRef[$i]{END_};
        }
    }
    close $cnvSimFileH;
    loadData($cnvSimFile, $cnvSimTable, ",", "CHROMOSOME, START_, END_"); 
    return $cnvSimTable; 
}

sub commitCNVSim{
    my ($cnvSimFileH, $cnvChrom, $cnvStart, $cnvEnd) = @_;
    print $cnvSimFileH "$cnvChrom,$cnvStart,$cnvEnd\n";
}

sub mapAllCNVS{
    my ($CNVsTable, $resultsRef) = @_;
    my $binSize = 1E4;
    my $windowSize = 1E6;
    my $halfWindowSize = $windowSize / 2;
    my $mapTable = "$CNVsTable\_MAP";
    dropTable($mapTable);
    runSQL("CREATE TABLE $mapTable (CHROMOSOME NUMBER, POSITION NUMBER, COVERAGE NUMBER)");    
    foreach my $chrom(1..nChrom()){
        my %map;   
        runSQL("SELECT ((trunc(START_/$binSize) * $binSize) + $binSize) AS LOWBIN, (trunc(END_/$binSize) * $binSize) AS HIGHBIN
                FROM $CNVsTable
                WHERE CHROMOSOME = $chrom",
                \my($lowBin, $highBin));
        while (fetchRow()){
            for (my $bin = $lowBin - $halfWindowSize; $bin <= $highBin + $halfWindowSize; $bin += $binSize){
                $map{$bin}++;
            }
        }
        my $mapFile = "$mapTable.csv";  #load the results into map table
        open my $mapFileH, ">", $mapFile;
        foreach my $bin (keys %map){
            print $mapFileH join(",", ($chrom, $bin, $map{$bin}))."\n";
        }
        close $mapFileH;
        loadData($mapFile, $mapTable, ",", "CHROMOSOME, POSITION, COVERAGE");         
    }
    my $hsTable = "$CNVsTable\_HS";
    dropTable($hsTable);
    runSQL("CREATE TABLE $hsTable (CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER, MAXCOVERAGE NUMBER)"); 
    foreach my $chrom(1..$refSeqs{$param{refSeqBase}}{nChrom}){
        runSQL("SELECT POSITION, COVERAGE
                FROM $mapTable
                WHERE CHROMOSOME = $chrom AND COVERAGE > 1
                ORDER BY POSITION", \my($pos, $coverage));  
        fetchRow();
        $pos or next;
        my ($start, $maxCoverage, $prevPos) = ($pos, $coverage, $pos);
        my $hsFile = "$hsTable.csv";  
        open my $hsFileH, ">", $hsFile;
        while (fetchRow()){
            if(($pos - $prevPos) > $binSize){ 
                print $hsFileH join(",", ($chrom, $start, $prevPos, $maxCoverage))."\n";
                $start = $pos;
                $maxCoverage = $coverage;  
            }
            $maxCoverage >= $coverage or $maxCoverage = $coverage;
            $prevPos = $pos;
        } 
        print $hsFileH join(",", ($chrom, $start, $prevPos, $maxCoverage))."\n";
        close $hsFileH;
        loadData($hsFile, $hsTable, ",", "CHROMOSOME, START_, END_, MAXCOVERAGE"); 
    }
    runSQL("SELECT MAXCOVERAGE FROM $hsTable", \my($maxCov));
    while(fetchRow()){$$resultsRef{$maxCov}++};   
}


1;
