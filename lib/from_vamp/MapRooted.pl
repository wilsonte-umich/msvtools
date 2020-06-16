#!/usr/bin/perl -w
use strict;
use warnings;

#########################################################################
#MapRooted.pl finds read pairs for which one and only one read maps to 
#a group of repeats or other elements found within a database specified 
#by parameter rooted.  Any partner reads which do not also map to this
#database are then mapped to the full reference genome to establish their
#positions.  Partner reads are stored in a Frags table as "one-ended" Fragments.
#########################################################################

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %bowtieFields %gffFields));

my (%inDirs, %outDirs, %inReadFiles, %outReadFiles, %inMapFiles, %outMapFiles);
my (%alignedFiles, %notAlignedFiles, %alignedIDs);
my $aligned = "aligned.fa";
my $notAligned = "not_aligned.fa";  
my ($pairsTable, %stats, $inFragsTable, $outFragsTable);
my %strands = ('+' => 1, '-' => 2, 1 => 'plus', 2 => 'minus');

sub getmapRootedDirs {
    my ($inSample, $outSample) = @_;
    getDirectories($inSample, \%inDirs);
    getDirectories($outSample, \%outDirs);
    getReadFiles($inSample, 'purgedFasta', \%inReadFiles);
    getReadFiles($outSample, 'purgedFasta', \%outReadFiles);
    getMapFiles($inSample, \%inMapFiles);
    getMapFiles($outSample, \%outMapFiles);     
    $inMapFiles{read1} .= ".$param{rooted}";
    $inMapFiles{read2} .= ".$param{rooted}";
    $outReadFiles{pairedReads} = "$param{inputPath}/$outSample/$outSample\_$param{rooted}\_pairedReads.csv"; 
}

sub getMapRootedStats {
    my ($inSample, $outSample) = @_;
    status("copying $inSample statistics\n");
    my $inStatsTable = getTableName('Stats', $inSample);
    my $outStatsTable = newTable('Stats', $outSample);
    updateTable('Stats', $outStatsTable, "SELECT * FROM $inStatsTable");        
    getStatistics($outStatsTable, \%stats);          
}

sub mapRooted{
    my ($inSample, $outSample) = @_;
    my $paramFile = "$outSample\_map\_parameters.txt";     
    printParameters($paramFile, 'mapRooted', $inSample, $outSample);   

    getmapRootedDirs($inSample, $outSample);
    getMapRootedStats($inSample, $outSample);

    status("mapping $inSample to $param{rooted}\n");
        #always use pass for first round mapping, much faster than bowtie on small refSeqs   
        mapReads_PASS_Rooted($inSample, 'read1', \%inReadFiles, \%inMapFiles, $param{rooted}, 1);        
        moveAlignmentFiles('read1');        
        mapReads_PASS_Rooted($inSample, 'read2', \%inReadFiles, \%inMapFiles, $param{rooted}, 1);
        moveAlignmentFiles('read2');    
    
    status("writing $outSample partner read file as read1\n"); 
        unlink ($outReadFiles{read1}, $outReadFiles{read2}, $outReadFiles{pairedReads}); 
        my $oneEnded =  writeOneEnded($alignedFiles{read1}, $notAlignedFiles{read2})
                      + writeOneEnded($alignedFiles{read2}, $notAlignedFiles{read1});
        status("  $oneEnded pairs had one and only one read mapped to $param{rooted}\n\n");    
    
    status("mapping $outSample partner reads to $param{refSeq}\n");    
        #read1 = partner/unrooted read => run standard mapping
        #uses whichever mapping algorithm the user specifies, but should use bowtie, faster and smaller footprint for this purpose
        my $maxDiscHold = $param{maxDisc};
        mapReads($outSample, $outSample, 'read1');  
        $param{maxDisc} = $maxDiscHold;  #since bowtie will reduce maxDisc to 3, but later steps should use mapReads_PASS_Rooted value
            
 
    findNewRoots($inSample, $outSample);
 
    # status("extracting partner reads directly into fragments\n");
        # $inFragsTable = getTableName('Frags', $inSample);
    	# $outFragsTable = newTable('Frags', $outSample);
        # runExtraction_Rooted();
        # markRootedNormal();  #pick best Normal?
        # purgeDuplicateFragments($outFragsTable);
}

sub mapReads_PASS_Rooted{
    my ($sample, $read, $readFilesRef, $mapFilesRef, $refSeq, $maxHits) = @_;
    status("\nmapping $read\n");    
    my $inputFile = $$readFilesRef{$read};  #reads as purged fasta       
    -e $inputFile or die "$inputFile does not exist";     
    my $iQuery = "-i $inputFile";
    my $oGFF = "-o $$mapFilesRef{$read}"; 
    my $output  = "-gff -info_gff -aligned -not_aligned"; #output in gff format, including mismatch info
    my $refSeqFile = "$param{refSeqPath}/$refSeq/$refSeq.fa";
    my $refdb = "$refSeqFile.$param{longWord}.db";
    -e $refdb or createDatabase($refSeq);
    my $R = "-R $refdb";  #reference sequence derived from inputs
    my $g = "-g 1";   #allow alignment gaps
    my $maxBestHits = "-max_best_hits $maxHits";  
    my $pstPath = "$param{passPath}PST/W$param{shortWord}M1m0G0X0.pst";
    my $matchThreshold = int((($param{readLength}-$param{maxDisc})/$param{readLength})*100);        
    my $fid = "-fid $matchThreshold";  #allowed discrepancies communicated as % match        
    my $pstThreshhold = (2*$param{shortWord})-$param{maxDisc};  #same as table provided by pass
    my $pst = "-pst $pstPath $pstThreshhold";         
    unlink ($$mapFilesRef{$read}, $notAligned, $aligned);  #clear any previous mapping output     
    runPass(($iQuery, $R, $fid, $g, $maxBestHits, $pst, $output, $oGFF));
}

sub mapReads_Bowtie_Rooted{
    my ($sample, $read, $readFilesRef, $mapFilesRef, $refSeq, $maxHits) = @_;
    status("  mapping $read\n\n");     
    if ($param{maxDisc} > 3){
        status("Bowtie allows up to 3 mismatches, changing maxDisc to 3...\n");
        $param{maxDisc} = 3;
    }
    my $inputFile = $$readFilesRef{$read};  
    -e $inputFile or die "$inputFile does not exist"; 
    -e "$param{refSeqPath}/$refSeq/$refSeq.1.ebwt" or createBowtieIndex($refSeq); 
    my $bowtieCommand = "$param{bowtiePath}bowtie "; 
    $bowtieCommand .= "-f -B 1 ";  #fasta, report refSeq 1-referenced
    $bowtieCommand .= "--un $notAligned "; 
    $bowtieCommand .= "--al $aligned "; 
    $bowtieCommand .= " -v $param{maxDisc} ";
    $bowtieCommand .= " -k $maxHits ";        
    $bowtieCommand .= " $param{refSeqPath}/$param{refSeq}/$param{refSeq} ";
    $bowtieCommand .= " $inputFile ";
    $bowtieCommand .= " $$mapFilesRef{$read} "; 
    $bowtieCommand .= " > $$mapFilesRef{$read}.stdout";
    unlink ($$mapFilesRef{$read}, "$$mapFilesRef{$read}.stdout", $notAligned, $aligned);  #clear any previous mapping output
    runBowtie($bowtieCommand);        
}

sub moveAlignmentFiles{
    my ($read) = @_;
    unlink $inMapFiles{$read};
    $alignedFiles{$read} = "$outDirs{$read}/$aligned.$param{rooted}.$read";
    system("mv $aligned $alignedFiles{$read}"); 
    $notAlignedFiles{$read} = "$outDirs{$read}/$notAligned.$param{rooted}.$read";
    system("mv $notAligned $notAlignedFiles{$read}"); 
}

sub writeOneEnded{
    my ($alignedFile, $notAlignedFile) = @_;
    my %aligned;    
    open my $alignedH, "<", $alignedFile;
    while (<$alignedH>){ 
        if ($_ =~ m/^>(.+)\n/){ 
            #$aligned{$1}++;
            $aligned{$1} = <$alignedH>; 
            chomp $aligned{$1}; 
        }
    }
    close $alignedH;
    open my $notAlignedH, "<", $notAlignedFile;    
    open my $partnerH, ">>", $outReadFiles{read1}; 
    open my $pairedReadsH, ">>", $outReadFiles{pairedReads}; 
    my $counter = 0;   
    while (<$notAlignedH>){
        if ($_ =~ m/^>(.+)\n/ and $aligned{$1}){
            $counter++;
            print $partnerH $_; #id line
            my $partnerSeq = <$notAlignedH>;
            print $partnerH $partnerSeq; #seq line 
            print $pairedReadsH "$1,$aligned{$1},$partnerSeq"; #pairID,alignedSequence,partnerSequence                  
        } 
    }
    close $notAlignedH;
    close $partnerH;
    close $pairedReadsH;
    return $counter;      
}

sub runExtraction_Rooted{
    status("  extracting\n");
    my $timeStamp = getTimeStamp();
    my $getLineSub = "getLine\_$param{mapType}";
    my $getLineSubRef = \&{$getLineSub};
    my $outFragsFile = "$outFragsTable\.csv";   
    open my $outFragsFileH, ">", $outFragsFile or die "could not open $outFragsFile";    
    open my $mapFileH, "<", $outMapFiles{read1} or die "could not open $outMapFiles{read1}";
    my ($pairID, $map, $prevPairID, @prevMaps);
    my $line = <$mapFileH>;
    ($prevPairID, $map) = &$getLineSubRef(\$line);   
    push @prevMaps, $map;
    while (<$mapFileH>){     
        ($pairID, $map) = &$getLineSubRef(\$_);
        if ($pairID) {
            if ($pairID eq $prevPairID){
                push @prevMaps, $map,
            } else {
                printRootedFrags($outFragsFileH, $timeStamp, $prevPairID, \@prevMaps);
                @prevMaps = ($map);
            }
            $prevPairID = $pairID; 
        }     
    }
    printRootedFrags($outFragsFileH, $timeStamp, $prevPairID, \@prevMaps);
    close $outFragsFileH;
    close $mapFileH;
    loadData($outFragsFile, $outFragsTable, ",", "FRAGMENTID, FRAGMENTTYPE,
                                            CHROMOSOME1, POSITION1, LENGTH1, STRAND1, DISCREPANCIES1, NHITS1,
                                            CHROMOSOME2, POSITION2, LENGTH2, STRAND2, DISCREPANCIES2, NHITS2, 
                                            FRAGMENTSIZE, PAIRID, NFRAGS,
                                            EVENTSIZE, STDEVNORMAL, ENDTOLERANCE,
                                            NSETSFRAG, NSETSPAIR"); 
}

sub printRootedFrags{
    my ($outFragsFileH, $timeStamp, $pairID, $mapsRef) = @_;
    my $fragID = ($pairID * 1E8) + $timeStamp;   
    my $nFrags = scalar @$mapsRef;            
    foreach my $pairMap(@$mapsRef){
        my ($chrom, $pos1, $length, $strand, $discs) = split(/:/, $$pairMap);
        my $pos2 = $pos1 + 10000;
        print $outFragsFileH join(",", $fragID, $types{Frags}{Rooted},
                                $chrom, $pos1, $length, $strand, $discs, $nFrags,
                                0, $pos2, $length, $strand, $discs, $nFrags,
                                0, $fragID, $nFrags,
                                $stats{modeNormal}, $stats{stDevNormal}, $stats{maxNormal},  
                                1, 0 )."\n";                     
    }     
}

sub markRootedNormal{
    status("  marking Normal\n");
    my $normalPairIDsSQL = "SELECT Trunc(PAIRID/100000000)*100000000 
                            FROM $inFragsTable 
                            WHERE FRAGMENTTYPE = $types{Frags}{Normal}
                               OR FRAGMENTTYPE = $types{Frags}{ReverseNormal}";
    runSQL("UPDATE $outFragsTable SET FRAGMENTSIZE = -1
            WHERE   Trunc(PAIRID/100000000)*100000000 IN ($normalPairIDsSQL)");
}

sub findNewRoots { 
    my ($inSample, $outSample) = @_;
    #my $outSample = "$inSample\_$param{rooted}";
    #getmapRootedDirs($inSample, $outSample);
    #my $inStatsTable = getTableName('Stats', $inSample);  
    #getStatistics($inStatsTable, \%stats);   
    getRootedPairedReads($outReadFiles{pairedReads}, \my%pairedReads); #read from prior analysis
    status("loading $param{refSeq} chromosome sequences...\n"); #subs found in findHomozygous.pl
    my $lineSize = getLineSize(); 
    my $chromLines = loadChromLines(); 
    seekExistingRootPairs(\%pairedReads, $lineSize, $chromLines, \my%newRootPairs); 
    printNewRootPairs(\%newRootPairs, $outSample);
}

sub getRootedPairedReads { #load root/nonRoot pairs and store in hash, key = pairID
    my ($pairedReadsFile, $pairedReads) = @_;
    status("loading paired reads (one root, one non-root per pair)...\n");
    open my $pairedReadsH, "<", $pairedReadsFile or die "could not open $pairedReadsFile\n"; 
    while (my $line = <$pairedReadsH>){
        chomp $line;
        my ($pairID, $rootRead, $nonRootRead) = split(",", $line);
        $$pairedReads{$pairID} = [$rootRead, $nonRootRead];
    }
    close $pairedReadsH;
}

sub seekExistingRootPairs { #for every nonRoot map, establish whether root read maps within Normal Fragment distance
    my ($pairedReads, $lineSize, $chromLines, $newRootPairs) = @_;
    status("seeking existing root pairs...\n");
    my $getLineSub = "getLine\_$param{mapType}";
    my $getLineSubRef = \&{$getLineSub};
    open my $mapFileH, "<", $outMapFiles{read1} or die "could not open $outMapFiles{read1}";
    #################
    #my $counter = 0;
    #################
    while (<$mapFileH>){     
        my ($pairID, $map1) = &$getLineSubRef(\$_);
        if ($pairID) {
            my ($rootRead, $nonRootRead) = @{$$pairedReads{$pairID}};
            my ($chrom1, $pos1, $length1, $strand1, $discs1) = split(/:/, $$map1); #genome map of the nonRootRead, strand corrected for circles
            my ($minPos, $maxPos) = getCandidateRootCoordinates($pos1, $strand1);          
            my $rootSegment = getDNASegment($chrom1, $minPos, $maxPos, $lineSize, $chromLines); 
            $rootSegment or next;
            #my ($pos2, $length2, $strand2, $discs2) = mapCandidateRoot($rootRead, $rootSegment);  #bowtie. MUCH MUCH slower!
            my ($pos2, $length2, $strand2, $discs2) = runPerlRootMatch($rootRead, \$rootSegment, $strand1); #all mtaching done in Perl 
            my $fType;
            if ($strand2 and $strand1 != $strand2) {
                $fType = $types{Frags}{Normal};
                correctRootPosition(\$pos2, $pos1, $strand1); 
                #NOTE: do NOT correct so that pos1 < pos2
                #this is different than standard Frags table
                #here, keep the original order so that it is known that read1 = genome partner read, read2 = root read
            }
            defined $fType or ($fType, $pos2, $length2, $strand2, $discs2) = 
                              ($types{Frags}{Rooted}, $pos1 + 10000, $length1, $strand1, $discs1);
            #store all pair mapping results in hash to allow purging of sr with Normal later
            push @{$$newRootPairs{$pairID}{$fType}}, [$chrom1, $pos1, $length1, $strand1, $discs1,
                                                            0, $pos2, $length2, $strand2, $discs2];                                                 
            #################
            #$counter ++;
            #$counter >= 100 and return;
            #################
        }     
    }   
    close $mapFileH;
}

sub printNewRootPairs { #export paired root mappings to frags table with preference to Normal
    my ($newRootPairs, $outSample) = @_;
    status("parsing and loading data into frags table...\n");
    my $timeStamp = getTimeStamp();
    #my $outFragsTable = newTable('Frags', "$outSample\_newRoots");  #########change this if newRoots becomes standard mapRooted
    my $outFragsTable =  newTable('Frags', $outSample);
    my $outFragsFile = "$outFragsTable.csv";
    open my $outH, ">", $outFragsFile or die "could not open $outFragsFile\n";
    foreach my $pairID(keys %$newRootPairs){
        my $fragID = ($pairID * 1E8) + $timeStamp; 
        my $sr = $$newRootPairs{$pairID}{$types{Frags}{Rooted}}; #single nonRoot reads, rootRead did not match reference
        my $nl = $$newRootPairs{$pairID}{$types{Frags}{Normal}}; #rootRead matched reference = Normal fragments
        my ($frags, $fragType) = ($sr, $types{Frags}{Rooted});
        $nl and ($frags, $fragType) = ($nl, $types{Frags}{Normal}); #keep all Normal and discard all SR if any Normal found
        my $nFrags = scalar(@$frags);
        foreach my $frag (@$frags){ 
            #discs1 as per vamp standard, but discs2 here is ONLY the discrepancy count for Normal pairs
            my ($chrom1, $pos1, $length1, $strand1, $discs1, $chrom2, $pos2, $length2, $strand2, $discs2) = @$frag;
            print $outH join(",", $fragID, $fragType,
                                $chrom1, $pos1, $length1, $strand1, $discs1, $nFrags,
                                $chrom2, $pos2, $length2, $strand2, $discs2, $nFrags,
                                0, $fragID, $nFrags,
                                $stats{modeNormal}, $stats{stDevNormal}, $stats{maxNormal},  
                                1, 0 )."\n";        
        }
    } 
    close  $outH;
    loadData($outFragsFile, $outFragsTable, ",", "FRAGMENTID, FRAGMENTTYPE,
                                            CHROMOSOME1, POSITION1, LENGTH1, STRAND1, DISCREPANCIES1, NHITS1,
                                            CHROMOSOME2, POSITION2, LENGTH2, STRAND2, DISCREPANCIES2, NHITS2, 
                                            FRAGMENTSIZE, PAIRID, NFRAGS,
                                            EVENTSIZE, STDEVNORMAL, ENDTOLERANCE,
                                            NSETSFRAG, NSETSPAIR");   
    #thus, have purged rooted with normal, otherwise kept all rooted
    #have not picked best Normals
}

sub getCandidateRootCoordinates { #use Normal peak to grab just the "sweet spot" of genome where root read should map
    my ($pos1, $strand1) = @_; #strand1 already corrected for isCircles
    my ($minPos, $maxPos);
    if ($strand1 == 1) { #forward nonRootRead, grab segment to the right
        ($minPos, $maxPos) = ($pos1 + $stats{minNormal}, $pos1 + $stats{maxNormal});  
    } else { #reverse nonRootRead, grab segment to the left
        ($minPos, $maxPos) = ($pos1 - $stats{maxNormal}, $pos1 - $stats{minNormal});         
    }
    $minPos >= 0 or $minPos = 0;    
    $maxPos >= 0 or $maxPos = 0;   
    return ($minPos, $maxPos);
}

sub correctRootPosition { #reverse the above process to restore segment coordinates to genome coordinates
    my ($pos2, $pos1, $strand1) = @_;
    if ($strand1 == 1) {
        $$pos2 += $pos1 + $stats{minNormal};  
    } else {
        $$pos2 += $pos1 - $stats{maxNormal};        
    }   
}

sub runPerlRootMatch { #execute root strand matching; don't get full discrepancy, but much faster matching
    my ($rootRead, $rootSegment, $nonRootStrand) = @_;
    $rootRead = "\U$rootRead"; #ensure both sequences uppercase for accurate matching
    $$rootSegment = "\U$$rootSegment";
    my $rootStrand = ($nonRootStrand % 2) + 1; #only look for Normal strand orientations
    my $workingStrand = $rootStrand;
    $param{isCircles} and $workingStrand = ($workingStrand % 2) + 1;  #remembering to invert if isCircles
    $workingStrand == 1 or $rootRead = reverseComplement($rootRead);
    my $readLength = length($rootRead);
    my $targetLength = length($$rootSegment);
    #all return pos corrected 1-referenced coordinates
    if ($$rootSegment =~ m/^(.*)$rootRead/) { #check exact match first for speed
        return (length($1) + 1, $readLength, $rootStrand, 0);
    } else { #check each substring of target for partial match within maxDisc limit
        my $maxI = int( ($targetLength-(2*$readLength)) / 2 );
        my $center = int($targetLength/2);
        foreach my $i(0..$maxI){
            #check both directions working from center of peak, where match likelihood is highest
            my ($pos, $nDiscs) = getRootMatchDiscs($i, 1, $center, $rootSegment, $rootRead, $readLength);
            $nDiscs <= $param{maxDisc} and return ($pos + 1, $readLength, $rootStrand, $nDiscs);
            ($pos, $nDiscs) = getRootMatchDiscs($i, -1, $center, $rootSegment, $rootRead, $readLength);
            $nDiscs <= $param{maxDisc} and return ($pos + 1, $readLength, $rootStrand, $nDiscs);
        }   
    }
    return undef; #no matches
}

sub getRootMatchDiscs { #extract target substring and count its mismatches to rootRead (does NOT permit gaps)
    my ($i, $direction, $center, $rootSegment, $rootRead, $readLength) = @_;
    my $pos = $center + ($i * $direction); 
    my $target = substr($$rootSegment,$pos,$readLength);
    my $nDiscs = $readLength - (($rootRead ^ $target) =~ tr/\0/\0/);
    return ($pos, $nDiscs); 
    #http://www.perlmonks.org/?node_id=371799  VERY SLICK!
    #This takes advantage of the fact that the exclusive-or (^) of two characters is zero iff 
    #the two characters are the same; regex tr returns the count of zeros/nulls i.e. perfect matches
    #Thoroughly tested and works well and fast for counting number of discrepancies between equal length strings
}

# sub mapCandidateRoot { #uses bowtie
    # my ($rootRead, $rootSegment) = @_;   
    # if ($param{maxDisc} > 3){
        # status("Bowtie allows up to 3 mismatches, changing maxDisc to 3...\n");
        # $param{maxDisc} = 3;
    # }
    # my $ewbt = 'candidateRoot';    
    # my $indexCommand = "$param{bowtiePath}bowtie-build -q -c $rootSegment $param{bowtiePath}/indexes/$ewbt"; 
    # system($indexCommand);
    # my $mapCommand = "$param{bowtiePath}bowtie --quiet -c -B 1 "; 
    # $mapCommand .= " -v $param{maxDisc} -k 1 ";      
    # $mapCommand .= " $ewbt $rootRead";
    # my $map = qx/$mapCommand/;
    # $map or return;
    # chomp $map;
    # my @in = split(/\t/, $map);
    # my $strand = $strands{$in[$bowtieFields{strand}]};
    # $param{isCircles} and $strand = ($strand % 2) + 1;    
    # my $pos = $in[$bowtieFields{position}];     
    # if ($pos < 0){$pos = 0}
    # my ($discs) = getDiscrepancies_bowtie($in[$bowtieFields{mismatches}]);
    # my $length = length($in[$bowtieFields{sequence}]);
    # return ($pos, $length, $strand, $discs);      
# }

1;
