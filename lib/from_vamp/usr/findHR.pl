#!/usr/bin/perl
use strict;
use warnings;
use DBI;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs));

my ($trackedSample, $lineSize, $chromLinesRef, $queryStart, $queryEnd, $query, $subjectStart, $subjectEnd, $subject);
my ($setID, $setType, $overlap, $chrom1, $pos1, $strand1, $chrom2, $pos2, $strand2);
my $flanking = 1500;
my (%segments, %segment1, %segment2);
my ($counter, $timeStamp, $outputFileH, %counts);

sub findHR{
    ($trackedSample) = @_;
    #$lineSize = getLineSize();
    #$chromLinesRef = loadChromLines($param{refSeq});
    my $outputTable = "AAA_$trackedSample\_HR_HITS";
    my $outputFile = "$outputTable.csv";
    open $outputFileH, ">", $outputFile;
    runSQL(fetchSetsSQL(), \($setID, $setType, $overlap, $chrom1, $pos1, $strand1, $chrom2, $pos2, $strand2));
    while(fetchRow()){
        getDNASegmentsHR();
        blastSegmentsHR();
    }
    close $outputFileH;
    runSQL("DELETE FROM $outputTable WHERE 1 = 1");
    loadData($outputFile, $outputTable, ',', "SETID, SETTYPE, BLASTCUTOFF, 
                                                              QUERYSTART, QUERYEND, SUBJECTSTART, SUBJECTEND, 
                                                              EVALUE, PERCENTIDENTITY, SCORE, LENGTH, OFFSET, OVERLAP");
}

sub fetchSetsSQL{
    my $setsTable = getTableName('Sets', $trackedSample);
    my $marksTable = getTableName('Marks', $trackedSample);
    my $fragsTable = getTableName('Frags', $trackedSample);
    my $marksSQL = "SELECT SETID FROM $marksTable WHERE MARK = 1 AND USERINITIALS = 'TEW' GROUP BY SETID";
    my $setsSQL = "SELECT SETID, SETTYPE, CHROMOSOME1, DECODE(CHROMOSOME2, 0, CHROMOSOME1, CHROMOSOME2) CHROMOSOME2,
                           SPANSTART, SPANEND, OVERLAPSTART, OVERLAPEND, STRAND1, STRAND2,
                           DECODE(STRAND1, 1, OVERLAPSTART, SPANSTART) POSITION1,                         
                           DECODE(STRAND2, 2, OVERLAPEND, SPANEND) POSITION2                           
                    FROM $setsTable
                    WHERE (SETTYPE = 1 OR SETTYPE = 4 OR SETTYPE = 8 OR SETTYPE = 16)
                    ";
    my $markedSetsSQL =  "SELECT s.* FROM ($setsSQL) s, ($marksSQL) m WHERE s.SETID = m.SETID";     
    my $modeNormal = "(f.FRAGMENTSIZE + f.EVENTSIZE)";                                  
    my $fragExcess = "(Abs(s.POSITION1 - f.POSITION1) + Abs(s.POSITION2 - f.POSITION2))";
    my $fragOverlap = "$modeNormal - $fragExcess";   
    my $overlap = "Trunc(Avg($fragOverlap))";               
    my $overlapSQL = "SELECT s.SETID, s.SETTYPE, $overlap OVERLAP,
                             s.CHROMOSOME1, s.POSITION1, s.STRAND1,
                             s.CHROMOSOME2, s.POSITION2, s.STRAND2
                        FROM $fragsTable f, ($markedSetsSQL) s    
                        WHERE f.CHROMOSOME1 = s.CHROMOSOME1
                          AND DECODE(f.CHROMOSOME2, 0, f.CHROMOSOME1, f.CHROMOSOME2) = s.CHROMOSOME2
                          AND f.POSITION1 >= s.SPANSTART AND f.POSITION1 <= s.OVERLAPSTART
                          AND f.POSITION2 >= s.OVERLAPEND AND f.POSITION2 <= s.SPANEND
                          AND f.STRAND1 = s.STRAND1
                          AND f.STRAND2 = s.STRAND2
                          AND f.NSETSFRAG > 0
                          AND f.FRAGMENTTYPE = s.SETTYPE
                        GROUP BY s.SETID, s.SETTYPE, 
                                 s.CHROMOSOME1, s.POSITION1, s.STRAND1, 
                                 s.CHROMOSOME2, s.POSITION2, s.STRAND2";                       
    return $overlapSQL;                           
}

sub getDNASegmentsHR{
    getDNASegmentHR($chrom1, $pos1, $strand1, 1, \%segment1);
    getDNASegmentHR($chrom2, $pos2, $strand2, 2, \%segment2);
    %segments = (query => \%segment1, subject => \%segment2); 
    $strand1 == 2 and %segments = (query => \%segment2, subject => \%segment1); 
}

sub getDNASegmentHR{
    my ($chrom, $pos, $strand, $side, $segmentRef) = @_;
    $$segmentRef{chrom} = $chrom;
    my $direction = 1;
    $strand == 2 and $direction = -1;
    ($$segmentRef{start}, $$segmentRef{end}) = ($pos, $pos + ($direction * $flanking)); 
    $$segmentRef{start} > $$segmentRef{end} and ($$segmentRef{start}, $$segmentRef{end}) = ($$segmentRef{end}, $$segmentRef{start}); 
    #$$segmentRef{sequence} = getDNASegment($chrom, $$segmentRef{start}, $$segmentRef{end}, $lineSize, $chromLinesRef);
    $$segmentRef{sequence} = getDNASegmentHR2($chrom, $$segmentRef{start}, $$segmentRef{end});
    $strand1 == $strand2 and $side == 2 and $$segmentRef{sequence} = reverseComplement($$segmentRef{sequence});
}

sub getDNASegmentHR2{
    my ($chrom, $start, $end) = @_;
    my $length = $end - $start + 1;
    my $file = "$param{refSeqPath}/$param{refSeqBase}/$reverseRefSeqs{$param{refSeqBase}}{$chrom}.seq";
    open my $fileH, "<", $file;   
    seek($fileH, $start, 0);
    my $bytes = read($fileH, my $sequence, $length);
    $bytes or print "SEQUENCE RETRIEVAL FAILURE: $chrom, $start, $end\n";
    return $sequence;
    close $fileH;
}

sub blastSegmentsHR{
    $timeStamp = getTimeStamp();    
    my $queryFile = getBlastSegmentFile('query');
    my $subjectFile = getBlastSegmentFile('subject');
    my $resultsFile = "$param{blastPath}/$timeStamp\_results.seq";    
    my $blast = "$param{blastPath}/blastn -query $queryFile -subject $subjectFile -out $resultsFile ";
    $blast .= "-strand plus -perc_identity 50 -outfmt \"10 qstart qend sstart send evalue pident score length\" "; 
    system($blast);  
    foreach (my $pi = 50; $pi <= 100; $pi += 5){
        open my $resultsFileH, "<", $resultsFile;  
        my %hits;    
        my @closestHit = ($setID, $setType, $pi);
        while($resultsFileH and my $line = <$resultsFileH>){
            chomp $line;
            my ($qStart, $qEnd, $sStart, $sEnd, $eValue, $pIdent, $score, $length) = split(",", $line);
            $pIdent >= $pi or next;
            $length >= 75 or next;
            my $absQStart = $segments{query}{start} + $qStart;    
            my $absQEnd = $segments{query}{start} + $qEnd;    
            my $absSStart = $segments{subject}{start} + $sStart;
            my $absSEnd = $segments{subject}{start} + $sEnd;    
            #$segments{query}{chrom} == $segments{subject}{chrom} and $absQStart == $absSStart and print "======"; #blast match of exact same bases
            $segments{query}{chrom} == $segments{subject}{chrom} and $absQStart == $absSStart and next; #blast match of exact same bases
            my $sizeDelta = $qStart - ($sStart - ($flanking - $overlap));   
            $hits{abs($sizeDelta)} = [$absQStart, $absQEnd, $absSStart, $absSEnd, $eValue, $pIdent, $score, $length, $sizeDelta, $overlap];
        }
        close $resultsFileH;      
        if (scalar(keys %hits)){
            my @hits = sort {$a <=> $b} keys %hits;
            my $closestHit = $hits{$hits[0]};        
            push @closestHit, @{$closestHit};
        } else {
            push @closestHit, (0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        }
        print $outputFileH join(",", @closestHit)."\n";       
    }
    unlink $resultsFile;
    unlink $queryFile;
    unlink $subjectFile;
}

sub getBlastSegmentFile{
    my ($type) = @_;
    my $file = "$param{blastPath}/$timeStamp\_$type.seq";     
    open my $fileH, ">", $file;    
    print $fileH $segments{$type}{sequence};
    close $fileH;
    return $file;
}

#    foreach my $chrom(1..nChrom()){
#        my $inFile = "$param{refSeqPath}/$param{refSeqBase}/$reverseRefSeqs{$param{refSeqBase}}{$chrom}.fa";
#        my $outFile = "$param{refSeqPath}/$param{refSeqBase}/$reverseRefSeqs{$param{refSeqBase}}{$chrom}.seq";
#        open my $inFileH, "<", $inFile;
#        open my $outFileH, ">", $outFile;    
#        while(my $line = <$inFileH>){
#            $line =~ m/^>/ and next;
#            $line =~ s/\s//g;
#            print $outFileH $line;
#        }
#        close $inFileH;
#        close $outFileH;
#    }
#exit;

1;



