#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %samFields %samFlags));

#callable parameters and commands in this script
defined $param{miVampSamples} or $param{miVampSamples} = 0; #$testSample should be in miVampSamples
defined $param{miVampNormSample} or $param{miVampNormSample} = 0;
defined $param{miSangerSamples} or $param{miSangerSamples} = 0;
defined $param{miSangerNormSample} or $param{miSangerNormSample} = 0;
defined $param{minSanger} or $param{minSanger} = 100;
defined $param{maxSanger} or $param{maxSanger} = 700;
defined $param{roots} or $param{roots} = 0;
defined $param{miMinFracUnq} or $param{miMinFracUnq} = 0.8;
defined $param{miMinNUnq} or $param{miMinNUnq} = 2;
defined $param{miMaxFracProm} or $param{miMaxFracProm} = 0.5;
defined $param{miMinFragsEachSide} or $param{miMinFragsEachSide} = 2;
defined $command{extractMouseInsertions} or $command{extractMouseInsertions} = ['multiThread', '48:00:00', 2000, 0];
defined $command{mapMouseInsertions} or $command{mapMouseInsertions} = ['singleThread', '48:00:00', 2000, 0];

my (@vampSamples, @sangerSamples, @roots, $fragsDir, $covDir, $plotsDir, $workingRoot);

sub extractMouseInsertions {
    my ($testSample) = @_;
    checkMouseInsertionParameters();
    checkMouseInsertionDirectories();
    my $outFile = "$param{inputPath}/$testSample.mouseInsertions.csv";
    open my $outFileH, ">", $outFile or die "could not open $outFile: $!\n";
    print $outFileH "Root,Chromosome,Center,Start,End,P2,P3,nUnq,fProm\n";
    foreach my $root(@roots){
        getMouseInsertionTables(\my%tables, $root, $testSample, []);
        getStatistics($tables{$testSample}{statsTable}, \my%testStats);
        fetchMouseInsertions($tables{$testSample}{rootFragsTable}, $tables{$testSample}{rootSetsTable}, $testStats{maxNormal}, $outFileH, $root);
    }
    close $outFileH;
}

sub mapMouseInsertions {
    my ($testSample, @refSamples) = @_;
    checkMouseInsertionParameters();
    checkMouseInsertionDirectories(); 
    getMouseInsertionDiscards(\my%discard);
    my $inFile = "$param{inputPath}/$testSample.mouseInsertions.csv";
    open my $inFileH, "<", $inFile or die "could not open $inFile: $!\n";
    my $outFile = "$inFile.targets.csv";
    open my $outFileH, ">", $outFile or die "could not open $outFile: $!\n";
    my $header = <$inFileH>;
    print $outFileH $header;
    my ($eventCount, $keptCount) = (0,0);
    while (my $line = <$inFileH>) {
        chomp $line;
        my ($root,$chrom,$center,$start,$end,$p2,$p3,$nUnq,$fProm) = split(",", $line);
        my $eventKey = "$chrom\_$start\_$end";
        $eventCount++;
        $discard{$eventKey} and next;
        $keptCount++;
        $workingRoot = $root;
        getMouseInsertionTables(\my%tables, $root, $testSample, @refSamples);
        my $fragsFile = fetchMouseInsertionFragments(\%tables, $chrom, $start, $end, $center, $testSample, @refSamples);    
        processMouseVamp(\my%vampCounts, $chrom, $start, $end, $center, $p2, $p3, $testSample, \my$testIsHomo, \my$refIsHomo);
        processMouseSanger(\my%sangerCounts, $chrom, $start, $end, $center, $p2, $p3, \$refIsHomo);
        !$testIsHomo and !$refIsHomo and print $outFileH "$line\n";
        plotMouseSamples($chrom, $start, $end, $center, \%vampCounts, $param{miVampNormSample}, 
                         $fragsFile, 'vamp', \@vampSamples, $testIsHomo, $refIsHomo);
        plotMouseSamples($chrom, $start, $end, $center, \%sangerCounts, $param{miSangerNormSample}, 
                         $fragsFile, 'sanger', \@sangerSamples, $testIsHomo, $refIsHomo);    
    }
    print "\nprocessed $keptCount of $eventCount events\n";
    close $inFileH;
    close $outFileH;
}
sub checkMouseInsertionParameters {
    ($param{miVampSamples} or $param{miSangerSamples}) or die "must specify either parameter miVampSamples or miSangerSamples\n";
    $param{miVampSamples} and @vampSamples = split(",", $param{miVampSamples});
    $param{miSangerSamples} and @sangerSamples = split(",", $param{miSangerSamples});
    $param{miVampSamples} and ($param{miVampNormSample} or die "must specific parameter miVampNormSample\n");
    $param{miSangerSamples} and ($param{miSangerNormSample} or die "must specific parameter miSangerNormSample\n");
    $param{roots} or die "must specify parameter roots\n";
    @roots = split(",", $param{roots});
}
sub checkMouseInsertionDirectories {
    $fragsDir = "$param{inputPath}/frags";
    -d $fragsDir or mkdir $fragsDir;
    $plotsDir = "$param{inputPath}/plots";
    -d $plotsDir or mkdir $plotsDir;
    $covDir = "$param{inputPath}/coverages";
    -d $covDir or mkdir $covDir;
}
sub getMouseInsertionDiscards {
    my ($discard) = @_;
    my $inFile = "$param{inputPath}/markList.csv";
    -e $inFile or return;
    open my $inFileH, "<", $inFile or die "could not open $inFile: $!\n";
    while(my $line = <$inFileH>){
        chomp $line;
        my ($eventKey, $mark) = split(",", $line);
        $mark == -1 and $$discard{$eventKey}++;   
    }
    close $inFileH;
}
sub getMouseInsertionTables {
    my ($tables, $root, @samples) = @_;
    foreach my $sample(@samples){ 
        $$tables{$sample} = {};
        getMouseSampleInsertionTables($$tables{$sample}, $root, $sample);
    }
}
sub getMouseSampleInsertionTables {
    my ($tables, $root, $sample) = @_;
    my $rootSample = "$sample\_$root";
    $$tables{statsTable} = getTableName('Stats', $sample);
    $$tables{fragsTable} = getTableName('Frags', $sample);
    $$tables{setsTable} = getTableName('Sets', $sample);
    $$tables{rootStatsTable} = getTableName('Stats', $rootSample);
    $$tables{rootFragsTable} = getTableName('Frags', $rootSample);
    $$tables{rootSetsTable} = getTableName('Sets', $rootSample);
}
sub fetchMouseInsertions {
    my ($rootFragsTable, $rootSetsTable, $maxNormal, $outFileH, $root) = @_;
    #separately extract forward and reverse one-ended sets
    #p1 and p2 are the forward set endpoints, p3 and p4 are the reverse set endpoints, each in ascending order   
    #thus the junction should be between p2 and p3, while p1 and p4 are at the outer flanks   
    foreach my $chrom(1..nChrom()){
        my $setsFSQL = "SELECT chromosome1 chromosome, spanstart p1, overlapstart p2, nfragssample nfs, nfragstotal nft
                         FROM $rootSetsTable
                         WHERE strand1 = 1
                           AND chromosome1 = $chrom
                           AND nfragssample >= $param{miMinFragsEachSide}";
        my $setsRSQL = "SELECT chromosome1 chromosome, spanstart p3, overlapstart p4, nfragssample nfs, nfragstotal nft 
                         FROM $rootSetsTable
                         WHERE strand1 = 2
                           AND chromosome1 = $chrom
                           AND nfragssample >= $param{miMinFragsEachSide}";
        #join forward and reverse sets, seeking closely spaced flanking pairs       
        #theoretically the outer flanks should be separated by no more than 2 * maxNormal, use 2.5 * maxNormal to give a bit of padding to each end 
        my $joinSQL = "SELECT f.chromosome, p2+round((p3-p2)/2) center, p1, p2, p3, p4
                       FROM ($setsFSQL) f, ($setsRSQL) r  
                       WHERE p4 - p1 <= 2.5 * $maxNormal
                         AND p4 - p1 >= 0
                         AND (f.nfs + r.nfs)/(f.nft + r.nft) > $param{miMinFracUnq}";
        #obtain aggregrate information about the fragments in the uber-sets   
        my $fragsSQL = "SELECT nfrags, nsetsfrag, position1, strand1
                         FROM $rootFragsTable
                         WHERE chromosome1 = $chrom
                           AND fragmenttype = $types{Frags}{Rooted}";  
        my $aggSQL = "SELECT chromosome, center, p2, p3,
                             sum(decode(f.nfrags,1,1,0)) nunq,
                             avg(decode(f.nsetsfrag,1,0,1)) fprom
                      FROM ($fragsSQL) f, ($joinSQL) s
                      WHERE ((f.position1 >= s.p1 and f.position1 <= s.p2 and f.strand1 = 1)
                         OR  (f.position1 >= s.p3 and f.position1 <= s.p4 and f.strand1 = 2))
                      GROUP BY chromosome, center, p2, p3";
        #apply uniqueness and quality filters
        #first experssion is sample percent sample v. reference uniqueness filter
        #second expression demands a certain number of unique fragments
        #third expression sets an upper limit on promiscuous reads
        my $qualSQL = "SELECT chromosome, center, center - $maxNormal start_, center + $maxNormal end_, p2, p3, nunq, fprom
                       FROM ($aggSQL)
                       WHERE nunq >= $param{miMinNUnq}
                         AND fprom < $param{miMaxFracProm}"; 
        runSQL($qualSQL, \my($chrom,$center,$start,$end,$p2,$p3,$nUnq,$fProm));
        while(fetchRow()){ print $outFileH "$root,$chrom,$center,$start,$end,$p2,$p3,$nUnq,$fProm\n"  }
    }
}
sub fetchMouseInsertionFragments {
    my ($tables, $chrom, $start, $end, $center, @samples) = @_;
    my $posID = "$chrom\_$start\_$end";
    my $outFile = "$fragsDir/$posID.csv";
    open my $outFileH, ">", $outFile or die "could not open $outFile: $!";  
    print $outFileH "y,forwardProm,forwardUnq,reverseProm,reverseUnq,normalRooted\n";
    my $y = -0.1;
    foreach my $sample(@samples){
        fetchMouseRootedFragments($$tables{$sample}{rootFragsTable}, $chrom, $start, $end, $center, $y, $outFileH);
        fetchMouseNormalFragments($$tables{$sample}{rootFragsTable}, $chrom, $start, $end, $center, $y, $outFileH);
        $y -= 0.1;
    }
    close $outFileH;  
    return $outFile;
}
sub fetchMouseRootedFragments {
    my ($rootFragsTable, $chrom, $start, $end, $center, $y, $outFileH) = @_;
    my $fragsSQL = "SELECT position1 - $center pos, strand1, decode(nsetsfrag,1,0,1) prom
                    FROM $rootFragsTable
                    WHERE chromosome1 = $chrom
                      AND position1 >= $start 
                      AND position1 <= $end
                      AND fragmenttype = $types{Frags}{Rooted}";
    runSQL($fragsSQL,\my($pos,$strand,$prom));
    while(fetchRow()){
        if ($strand == 1){
            if($prom){
                print $outFileH "$y,$pos,,,,\n";
            } else {
                print $outFileH "$y,,$pos,,,\n";
            } 
        } else{
            if($prom){
                print $outFileH "$y,,,$pos,,\n";
            } else {
                print $outFileH "$y,,,,$pos,\n";
            } 
        }
    }
}
sub fetchMouseNormalFragments {
    my ($rootFragsTable, $chrom, $start, $end, $center, $y, $outFileH) = @_;
    my $fragsSQL = "SELECT position2 - $center pos
                    FROM $rootFragsTable
                    WHERE chromosome1 = $chrom
                      AND position2 >= $start 
                      AND position2 <= $end
                      AND fragmenttype = $types{Frags}{Normal}";
    runSQL($fragsSQL,\my($pos));
    while(fetchRow()){ print $outFileH "$y,,,,,$pos\n" }
}
sub processMouseVamp { #vamp sample coverage data recovered from database as Normal Fragments
    my ($binCounts, $chrom, $start, $end, $center, $p2, $p3, $testSample, $testIsHomo, $refIsHomo) = @_; 
    foreach my $sample (@vampSamples){
        my $fragsTable = getTableName('Frags', $sample);
        runSQL("SELECT position1 low, position2 + length2 - 1 high
                FROM $fragsTable
                WHERE fragmenttype = $types{Frags}{Normal}
                  AND chromosome1 = $chrom
                  AND position1 <= $end
                  AND $start <= position2 + length2 - 1", #allows fragments overlapping start and end
                \my($low,$high)); 
        initializeMouseCounts($sample, $start, $end, \my%rawCounts, $binCounts);      
        while(fetchRow()){ processMouseNormal($low, $high, $start, $end, $rawCounts{$sample}, $$binCounts{$sample}) } #increment counter for all fragment positions
        my $isHomo = processMouseSample($start, $end, $p2, $p3, $rawCounts{$sample}, $$binCounts{$sample}); #normalize sample counts to mean position count   
        if($sample eq $testSample){
            $$testIsHomo = $isHomo;
        } else {
            $$refIsHomo = ($$refIsHomo or $isHomo);  
        }
    }  
}
sub processMouseSanger { #sanger sample coverage data recovered via ftp using samtools
    my ($binCounts, $chrom, $start, $end, $center, $p2, $p3, $refIsHomo) = @_;    
    fixBamSexChrom(\$chrom); #Sanger bams use 1,2,3...X,Y
    foreach my $sample (@sangerSamples){ 
        my $bamPath = "ftp://ftp-mouse.sanger.ac.uk/current_bams/$sample.bam";
        samtoolsView($bamPath, [[$chrom, $start, $end]], \my@samResults);
        initializeMouseCounts($sample, $start, $end, \my%rawCounts, $binCounts);
        foreach my $line(@samResults){
            $line =~ m/^\@/ and next;
            my @line = split("\t", $line);
            my $flag = $line[$samFields{FLAG}];
            getSamFlagBit($flag, 'isReverse') and next;
            getSamFlagBit($flag, 'nextIsReverse') or next;
            my $low = $line[$samFields{POS}];
            my $high = $line[$samFields{PNEXT}] + length($line[$samFields{SEQ}]);
            my $tLen = $high - $low; #recalculate tLen so that it reflects full fragment length
            ($tLen > $param{minSanger} and $tLen < $param{maxSanger}) or next; #rough enforcement of Normal Fragment length restriction
            processMouseNormal($low, $high, $start, $end, $rawCounts{$sample}, $$binCounts{$sample}) #increment counter for all fragment positions
        }
        my $isHomo = processMouseSample($start, $end, $p2, $p3, $rawCounts{$sample}, $$binCounts{$sample}); #normalize sample counts to mean position count
        $$refIsHomo = ($$refIsHomo or $isHomo);  
    }
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
sub initializeMouseCounts {
    my ($sample, $start, $end, $rawCounts, $binCounts) = @_;
    $$rawCounts{$sample} = {};
    $$binCounts{$sample} = {};
    initializeMouseCounts2($start, $end, $$rawCounts{$sample});
    initializeMouseCounts2($start, $end, $$binCounts{$sample}, $param{binSize});
}
sub initializeMouseCounts2 {
    my ($start, $end, $sampleCounts, $binSize) = @_;
    if($binSize){
        my ($lowBin, $highBin) = getMouseILimits($start, $end);
        for (my $i = $lowBin; $i <= $highBin; $i += $binSize){ $$sampleCounts{$i} = 0 }    
    } else {
        foreach my $i($start..$end){ $$sampleCounts{$i} = 0 }     
    }
}
sub getMouseILimits {
    my ($low, $high) = @_;   
    my $lowI = binMouseLimit($low);
    my $highI = binMouseLimit($high);
    return ($lowI, $highI);
}
sub binMouseLimit {
    my ($pos) = @_;
    return int(($pos/$param{binSize}) + 0.5) * $param{binSize};
}
sub processMouseNormal { #increment the sample count for every position covered by the incoming fragment
    my ($low, $high, $start, $end, $rawCounts, $binCounts) = @_;
    $low < $start and $low = $start;
    $high > $end and $high = $end;
    $high > $low or (print "processMouseNormal: bad position order: $low, $high, $start, $end\n" and return);
    foreach my $i($low..$high){ $$rawCounts{$i}++ }   
    my ($lowBin, $highBin) = getMouseILimits($low, $high);
    for (my $i = $lowBin; $i <= $highBin; $i += $param{binSize}){ $$binCounts{$i}++ }    
}
sub processMouseSample { #normalize each sample to its own mean coverage across region
    my ($start, $end, $p2, $p3, $rawCounts, $binCounts) = @_;
    my ($nPos, $sum);
    my ($lowBin, $highBin) = getMouseILimits($start, $end);
    for (my $i = $lowBin; $i <= $highBin; $i += $param{binSize}){ 
        $sum += $$binCounts{$i};
        $nPos++;
    }   
    my $mean = $sum / $nPos;
    $$binCounts{mean} = int($mean*10 + 0.5)/10;
    for (my $i = $lowBin; $i <= $highBin; $i += $param{binSize}){ 
        if ($mean){
            $$binCounts{norm}{$i} = $$binCounts{$i} / $mean; 
        } else {
            $$binCounts{norm}{$i} = -999;
        }
    } 
    $mean > 3 or return 0;
    $p2 > $p3 and ($p2,$p3) = ($p3,$p2);
    foreach my $i($p2..$p3){ $$rawCounts{$i} == 0 and return 1 }
    return 0;
}
sub plotMouseSamples { #normalize each sample to the indicated refSample and send to R for plotting
    my ($chrom, $start, $end, $center, $counts, $normSample, $fragFile, $source, $samples, $testIsHomo, $refIsHomo) = @_;
    my @plotSamples;
    foreach my $sample (@$samples){
       $$counts{$sample}{mean} > 3 or next; #demand minimal coverage
       $sample eq $normSample and next; #never plot the reference sample normalized to itself
       push @plotSamples, $sample;
    }
    scalar(@plotSamples) or return;
    my $homo;
    $testIsHomo and $refIsHomo and $homo = "Hom_In_Ref";
    $testIsHomo and !$refIsHomo and $homo = "Hom_No_Ref";
    !$testIsHomo and $refIsHomo and $homo = "Het_In_Ref";
    !$testIsHomo and !$refIsHomo and $homo = "Het_No_Ref";
    my $cDir = "$covDir/$homo";
    -d $cDir or mkdir $cDir;
    my $pDir = "$plotsDir/$homo";
    -d $pDir or mkdir $pDir;
    my $posID = "$chrom\_$start\_$end.$source";
    my $outFile = "$cDir/$posID.csv";
    open my $outFileH, ">", $outFile or die "could not open $outFile: $!";  
    print $outFileH "pos";
    foreach my $sample (@plotSamples){ print $outFileH ",$sample" }
    print $outFileH "\n";
    my ($lowI, $highI) = getMouseILimits($start, $end);
    for (my $i = $lowI; $i <= $highI; $i += $param{binSize}){ 
        my $pos = $i - $center; #express coordinates so that center = 0
        print $outFileH "$pos";
        foreach my $sample (@plotSamples){
            my $sampleNorm = $$counts{$sample}{norm}{$i};
            my $normSampleNorm = $$counts{$normSample}{norm}{$i};
            my $denom = $sampleNorm + $normSampleNorm;
            my $normValue = "";
            $denom and $normValue = $sampleNorm / $denom; #normalize to refSample as fraction sample
            print $outFileH ",$normValue"; 
        } 
        print $outFileH "\n";
    }
    close $outFileH;
    my $plotFile = "$pDir/$posID.jpg";
    my $chr = $reverseRefSeqs{$param{refSeq}}{$chrom};
    system("Rscript $param{vampPath}/bin/mouseInsertion/mouseInsertion.R $outFile $fragFile $plotFile $workingRoot $chr $start $end");  #plot with R
}

1;

