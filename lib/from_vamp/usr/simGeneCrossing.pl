#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs));

require "$param{vampPath}/bin/findIntersection.pl";

my ($majorDiv, $minorDiv) = ("=" x 50, "-" x 50);
my ($CNVsTable, $idField, $nIterations, $maxFragSize, $minGeneSize, $maxGeneSize, 
    $cnvSimTable, $genomeSize, $nChrom, $maxI, $chromsRef, $cnvsRef, $nCNVs, %results, $bootH);

sub findCNVGeneCrossings {
    ($CNVsTable, $idField, $nIterations, $maxFragSize, $minGeneSize, $maxGeneSize) = @_;
    ######################################################
    # status("loading table");
    # dropTable($CNVsTable);
    # runSQL("CREATE TABLE $CNVsTable (CNVID VARCHAR2(255), CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER)");    
    # loadData("$CNVsTable.csv", $CNVsTable, ',', 'CNVID, CHROMOSOME, START_, END_'); 
    # status("table loaded");
    # exit;
    ###################################################
    $minGeneSize or $minGeneSize = 0;
    $maxGeneSize or $maxGeneSize = 1E9;
    initializeCNVGeneCrossings();       
    status("analyzing the actual data...\n$majorDiv\n");   
    countGeneCrossings($CNVsTable, \%{$results{Actual}}, 1); 
    status("$majorDiv\nanalyzing $nIterations simulations...\n");
    
    open $bootH, ">", "$CNVsTable.boot.data.csv" or die "could not open boot data file\n";
    print $bootH "count,bin5\n";
    
    for my $iteration(1..$nIterations){
        status("$iteration ");
        generateCNVGeneCrossings(); 
        countGeneCrossings($cnvSimTable, \%{$results{Simulated}}, 0);    
    } 
    status("\nreporting results...\n");
    
    close $bootH;
    
    printCNVGeneCrossings();
}
sub initializeCNVGeneCrossings  {
    $nChrom = nChrom();
    $param{sex} and $param{sex} eq 'XX' and $nChrom--;
    my $chromInfoTable = "CHROMINFO_$param{refSeqBase}";
    runSQL("SELECT Sum(End_) FROM $chromInfoTable WHERE CHROMOSOME <= $nChrom", \($genomeSize));
    fetchRow();
    runSQL("SELECT CHROMOSOME, END_ FROM $chromInfoTable WHERE CHROMOSOME <= $nChrom ORDER BY CHROMOSOME");
    $chromsRef = fetchAllHashRef();
    $maxI = scalar(@$chromsRef) - 1;
    runSQL("SELECT $idField, CHROMOSOME, START_, END_, (END_ - START_) CNVSIZE 
            FROM $CNVsTable
            WHERE CHROMOSOME <= $nChrom
              AND END_ - START_ < $maxFragSize");   
    $cnvsRef = fetchAllHashRef();   
    $nCNVs = scalar @$cnvsRef;
    status("there are $nCNVs total CNV regions\n");
    $cnvSimTable = "$CNVsTable\_SIM";  
}
sub countGeneCrossings {
    my ($CNVsTable, $resultsRef, $isActual) = @_;
    my %cnvCounts;
    $isActual and status("CHROMOSOME\tSTART\tEND\tCNV_COUNT\n");
    my $count = 0;
    foreach my $chrom(1..$nChrom){
        my $geneSQL = "SELECT corrstart_, end_ 
                       FROM refgene_unq_$param{refSeq} 
                       WHERE CHROMOSOME = $chrom
                         AND geneSize >= minSize
                         AND geneSize < maxSize";
        my $cnvSQL = "SELECT $idField, start_, end_ FROM $CNVsTable WHERE CHROMOSOME = $chrom AND END_ - START_ < $maxFragSize";
        my $joinSQL = "SELECT c.$idField FROM ($cnvSQL) c, ($geneSQL) g WHERE c.start_ <= g.end_ AND g.corrstart_ <= c.end_";
        my $groupSQL = "SELECT $idField FROM ($joinSQL) GROUP BY $idField";
        my $countSQL = "SELECT count(*) FROM ($groupSQL)"; #number of unique CNVs crossing any gene
        runSQL($countSQL, \my$chromCount);
        fetchRow();
        $chromCount or $chromCount = 0;
        $count += $chromCount;
    }
    $$resultsRef{straight}{$count}++;  
    my $bin5 = int($count/5+0.5)*5;
    $$resultsRef{bin5}{$bin5}++;  
    
    $isActual or print $bootH "$count,$bin5\n";
}
sub generateCNVGeneCrossings {
    dropTable($cnvSimTable);
    runSQL("CREATE TABLE $cnvSimTable ($idField VARCHAR2(255), CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER)"); 
    my $cnvSimFile = "$cnvSimTable.csv";  
    open my $cnvSimFileH, ">", $cnvSimFile;   
    CNVSIM: foreach my $cnvRef(@$cnvsRef){
        my $id = $$cnvRef{$idField}; 
        my $cnvSize = $$cnvRef{CNVSIZE};
        TRY_AGAIN:
        my $cnvIndex = int(rand $genomeSize);
        my ($cumGenomeSize, $cnvChrom, $cnvStart, $cnvEnd) = (0);
        foreach my $i(0..$maxI){
            if($cnvIndex <= ($cumGenomeSize + $$chromsRef[$i]{END_})){ 
                $cnvChrom = $$chromsRef[$i]{CHROMOSOME};          
                $cnvStart = $cnvIndex - $cumGenomeSize;        
                $cnvEnd =  $cnvStart + $cnvSize;
                if(inGap($cnvChrom, $cnvStart, $cnvEnd)){goto TRY_AGAIN}  #simulated CNVs cannot reside in genome gaps
                if(inGap($cnvChrom, $cnvStart, $cnvEnd, 'omniQuad')){goto TRY_AGAIN}  #simulated CNVs cannot reside in array gaps
                if($cnvEnd <= $$chromsRef[$i]{END_}){              
                    commitCNVSim($cnvSimFileH, $id, $cnvChrom, $cnvStart, $cnvEnd); 
                } else { #simulated CNVs cannot wrap chromosomes
                    goto TRY_AGAIN;                   
                }
                next CNVSIM;
            }
            $cumGenomeSize += $$chromsRef[$i]{END_};
        }
    }
    close $cnvSimFileH;
    loadData($cnvSimFile, $cnvSimTable, ",", "$idField, CHROMOSOME, START_, END_"); 
}
sub printCNVGeneCrossings {  
    status("CNVs crossing genes\tIteration count\n");  
    printCNVGeneCrossings_('straight', 'Actual');
    printCNVGeneCrossings_('straight', 'Simulated');
    printCNVGeneCrossings_('bin5', 'Actual');
    printCNVGeneCrossings_('bin5', 'Simulated');
}
sub printCNVGeneCrossings_ {
    my ($binType, $simType) = @_;
    status("$binType\n"); 
    status("$simType\n"); 
    foreach my $count(sort {$a <=> $b} keys %{$results{$simType}{$binType}}){
        status("$count\t$results{$simType}{$binType}{$count}\n");     
    }
}

1;
