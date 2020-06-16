#!/usr/bin/perl -w
use strict;
use warnings;

use vars(qw(%param %types %fields %fieldNames %refSeqs));

my ($chrom, $eventsFileH, $eventID, $eventStart, $eventEnd, $eventType, @sets);
my (%sets, %querySets, %chains, %coverage, $bestRefStrand, $cncID, $cncsFileH, $nAlleles, $fmapTable);

sub assembleEvents{
    my ($sample) = @_;
    status("assembling events...\n");
    my $setsTable = getTableName('Sets', $sample);
    my $marksTable = getTableName('Marks', $sample); 
    $fmapTable = getTableName('FMap', $sample);
    (tableExists($setsTable) and tableExists($marksTable)) or die "could not find $setsTable and/or $marksTable";
    my $eventsTable = newTable('Events', $sample);
    my $cncsTable = newTable('CNCs', $sample);
    status("chrom ");
    foreach my $chrom_(1..nChrom()){  
        $chrom = $chrom_;
        status("$chrom ");
        my $eventsFile = "$eventsTable.csv";
        open $eventsFileH, ">", $eventsFile;
        my $cncsFile = "$cncsTable.csv";
        open $cncsFileH, ">", $cncsFile;
        runSQL("SELECT s.SETID, s.SETTYPE, 
                       s.SPANSTART, s.SPANEND, s.OVERLAPSTART, s.OVERLAPEND,
                       s.STRAND1, s.STRAND2
                FROM $setsTable s, $marksTable m
                WHERE s.SETID = m.SETID
                  AND (s.SETTYPE = $types{Sets}{Deletion} 
                       OR s.SETTYPE = $types{Sets}{Inversion} 
                       OR s.SETTYPE = $types{Sets}{Duplication})
                  AND s.CHROMOSOME1 = $chrom
                ORDER BY SPANSTART");
        my $setsRef = fetchAllHashRef();
        scalar @$setsRef or next;
        resetEvent();
        foreach my $setRef(@$setsRef){     
            if ($eventEnd and $$setRef{SPANSTART} > $eventEnd){ 
                commitEvent(); 
                resetEvent();
            }
            defined $eventStart or $eventStart = $$setRef{SPANSTART};
            $eventEnd >= $$setRef{SPANEND} or $eventEnd = $$setRef{SPANEND};
            $eventType = $eventType | $$setRef{SETTYPE};             
            push @sets, $setRef;
        }
        defined $eventStart and commitEvent(); 
        close $eventsFileH;
        close $cncsFileH;
        loadData($eventsFile, $eventsTable, ",", $fieldNames{Events}); 
        loadData($cncsFile, $cncsTable, ",", $fieldNames{CNCs});        
    }
    status("\n");
}

sub resetEvent{
    ($eventStart, $eventEnd, $eventType) = (undef, 0, 0); 
    @sets = ();
}
  
sub commitEvent{
    $eventID++;
    print $eventsFileH join(",", ($eventID, $eventType, $chrom, $eventStart, $eventEnd, scalar @sets))."\n";  
    %sets = ();
    %chains = ();
    %coverage = ();
    foreach my $pos ($eventStart..$eventEnd){ $coverage{$pos} = 0 }
    parseEventSets();          
    chainEventSets();
    parseChainSegments();
    parseCNVCalls();
}

sub parseEventSets{
    foreach my $setRef(@sets){
        my $leftLimit = $$setRef{OVERLAPSTART};
        $$setRef{STRAND1} == 2 and $leftLimit = $$setRef{SPANSTART};
        my $rightLimit = $$setRef{OVERLAPEND};         
        $$setRef{STRAND2} == 1 and $rightLimit = $$setRef{SPANEND};
        $sets{$$setRef{SETID}} = {setType => $$setRef{SETTYPE}, 
                                  leftStrand => $$setRef{STRAND1}, rightStrand => $$setRef{STRAND2}, 
                                  leftLimit => $leftLimit, rightLimit => $rightLimit};
    }
}   

sub chainEventSets{
    my %setsRemaining;
    foreach my $refStrand(1..2){
        %querySets = %sets;
        my ($success, $chain) = (1);
        while ($success and scalar keys %querySets){
            ($success, $chain) = tryChain($refStrand); 
            $success and push @{$chains{$refStrand}}, $chain; 
        }    
        $setsRemaining{$refStrand} = scalar keys %querySets;
    }
    $bestRefStrand = 1;
    $setsRemaining{2} < $setsRemaining{1} and $bestRefStrand = 2;
    $chains{$bestRefStrand} or return;
    $nAlleles = scalar @{$chains{$bestRefStrand}};
}

sub tryChain{
    my ($refStrand) = @_;  
    my $refPos = 0; $refStrand == 1 and $refPos = 1E9; 
    my ($firstRefStrand, @chain) = ($refStrand);
    while(1){     
        my ($searchStrand, %pos) = (($refStrand % 2)  + 1); #search for a strand of polarity opposite to the reference  
        foreach my $setID (keys %querySets){ 
            if($querySets{$setID}{leftStrand} == $searchStrand and checkLimitPos($refPos, $refStrand, $querySets{$setID}{leftLimit}) ){
                $pos{$querySets{$setID}{leftLimit}} = {setID => $setID, 
                                                       partnerStrand => $querySets{$setID}{rightStrand},
                                                       partnerLimit => $querySets{$setID}{rightLimit}}
            }
            if($querySets{$setID}{rightStrand} == $searchStrand and checkLimitPos($refPos, $refStrand, $querySets{$setID}{rightLimit}) ){
                $pos{$querySets{$setID}{rightLimit}} = {setID => $setID, 
                                                       partnerStrand => $querySets{$setID}{leftStrand},
                                                       partnerLimit => $querySets{$setID}{leftLimit}}
            }
        }   
        my @pos = sort {$a <=> $b} keys %pos; #the search direction is determined by the reference read polarity; 
        $refStrand == 2 or @pos = sort {$b <=> $a} keys %pos;   #reverse = forward search, forwarch = reverse search
        scalar @pos or last;   
        my $indexPos = $pos[0]; 
        my $setID = $pos{$indexPos}{setID};        
        delete $querySets{$setID};
        $refPos = $pos{$indexPos}{partnerLimit};    
        $refStrand = $pos{$indexPos}{partnerStrand};    
        if($firstRefStrand == 2){
            push @chain, "$indexPos:$refPos";
        } else {
            unshift @chain, "$refPos:$indexPos";
        }     
    }
    my $success = ($firstRefStrand == $refStrand and scalar @chain);
    my $chain = join("::", @chain);
    return ($success, $chain);  
}  

sub checkLimitPos{
    my ($refPos, $refStrand, $testPos) = @_;
    if ($refStrand == 2){ return $testPos > $refPos } else { return $testPos < $refPos }
}

sub parseChainSegments{
    foreach my $chain(@{$chains{$bestRefStrand}}){
        my @chain = split("::", $chain);    
        $bestRefStrand == 2 or @chain = reverse @chain;  
        my $segStart = $eventStart - 1; 
        $bestRefStrand == 2 or $segStart = $eventEnd + 1;         
        foreach my $i(0..(scalar(@chain) - 1)){
            my ($indexPos, $refPos) = split(":", $chain[$i]);
            $bestRefStrand == 2 or ($indexPos, $refPos) = ($refPos, $indexPos);
            my $segEnd = $indexPos;         
            fillSegmentPos($segStart, $segEnd);        
            $segStart = $refPos;      
        }    
        my $segEnd = $eventEnd + 1;
        $bestRefStrand == 2 or $segEnd = $eventStart - 1;
        fillSegmentPos($segStart, $segEnd);
    }    
}

sub fillSegmentPos{
    my ($segStart, $segEnd) = @_;
    $segStart > $segEnd and ($segStart, $segEnd) = ($segEnd, $segStart);
    foreach my $pos ($segStart..$segEnd){ $coverage{$pos}++ } 
}

sub parseCNVCalls{ 
    my $startPos;
    foreach my $pos(sort {$a <=> $b} keys %coverage){
        if($startPos and $coverage{$pos} != $coverage{$startPos}){
            commitCNVCall($startPos, $pos - 1);
            $startPos = undef;
        }
        $startPos or $startPos = $pos;
    }
    #shouldn't need to commit the last one...
}

sub commitCNVCall{
    my ($startPos, $endPos) = @_;
    my $copyNumber = $coverage{$startPos};
    $nAlleles == 1 and $copyNumber++; 
    $copyNumber == 2 and return; 
    my $nBins = int($endPos / $param{binSize}) - int($startPos / $param{binSize});
    my $normCov = 999;    
    if($nBins){
        runSQL("SELECT nvl(Sum(NORMALIZEDCOVERAGE),0)/$nBins FROM $fmapTable 
                WHERE CHROMOSOME = $chrom AND POSITION >= $startPos AND POSITION <= $endPos", \($normCov));   
        fetchRow();
    }
    $nAlleles == 1 and $copyNumber == 1 and $normCov < 0.25 and $copyNumber--;
    $cncID++;
    print $cncsFileH join(",", ($cncID, $eventID, $chrom, $startPos, $endPos, $copyNumber, $normCov))."\n"; 
}


1;


