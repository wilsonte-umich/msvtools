#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs));

#sub mapFixed{
#    my ($inSample) = @_;
#    my $outSample = $inSample."_fx";    
#    my $inFragsTable = getTableName('Frags', $inSample); 
#        
#    getMapFiles($inSample, \my %inMapFiles);     
#    getReadFiles($inSample, 'purgedFasta', \my %inReadFiles);      
#    getDirectories($outSample, \my %outDirs);      
#    -d $outDirs{sample} or mkdir $outDirs{sample};
#    -d $outDirs{read1} or mkdir $outDirs{read1};
#    -d $outDirs{read2} or mkdir $outDirs{read2}; 
#    getMapFiles($outSample, \my %outMapFiles_hold);  
#    $outMapFiles_hold{read1} .= ".hold";
#    $outMapFiles_hold{read2} .= ".hold";
#    getReadFiles($outSample, 'purgedFasta', \my %outReadFiles);
##    unlink $outReadFiles{read1};
##    unlink $outReadFiles{read2};   
#         
#    getNormalPairIDs($inFragsTable, \my %normPairIDs);
#    getFixedMaps(\%inMapFiles, \my %maps, \%normPairIDs);  
#    printFixedMaps(\%maps, \%outMapFiles_hold, \my %toMap);      
#    printFixedToMap(\%toMap, \%inReadFiles, \%outReadFiles);         

#    $param{mapType} = 'bowtie';
#    getMapFiles($outSample, \my %outMapFiles_new); 
#    mapFixedPartners(\%outReadFiles, \%outMapFiles_new, 'read1');
#    mapFixedPartners(\%outReadFiles, \%outMapFiles_new, 'read2');

#}

#sub getNormalPairIDs{
#    my ($inFragsTable, $normPairIDs) = @_;
#    status("getting normal pair IDs...\n");
#    runSQL("SELECT Trunc(PAIRID/1E8) 
#            FROM $inFragsTable 
#            WHERE FRAGMENTTYPE = $types{Frags}{Normal}
#               OR FRAGMENTTYPE = $types{Frags}{ReverseNormal}", \my($normPairID));
#    while(fetchRow()){ $$normPairIDs{$normPairID}++ }
#}

#sub getFixedMaps{
#    my ($inMapFiles, $maps, $normPairIDs) = @_;
#    status("getting fixed reads...\n");
#    $param{mapType} eq 'pass' or die "fixed analysis currently only supports PASS mapping data";
#    foreach my $read('read1','read2'){
#        open my $inFileH, "<", $$inMapFiles{$read};
#        while (my $line = <$inFileH>){    #####only written for PASS at present, easily adapted to reading bowtie or other mapping data     
#            $line =~ m/Name=(.+?);/ or next;
#            my $pairID = $1;
#            $$normPairIDs{$pairID} and next;
#            $line =~ m/Hits=(.*?);/ or next;
#            $$maps{$pairID}{$read}{nHits} = $1;
#            $1 <= $param{maxHits} and push @{$$maps{$pairID}{$read}{maps}}, $line;
#        }
#        close $inFileH;
#    }
#}

#sub printFixedMaps{
#    my ($maps, $mapFiles, $toMap) = @_;
#    status("printing previous good maps...\n");
#    open my $mapFileH1, ">", $$mapFiles{read1};
#    open my $mapFileH2, ">", $$mapFiles{read2};
#    foreach my $pairID(keys %$maps){
#        ($$maps{$pairID}{read1}{nHits} and $$maps{$pairID}{read2}{nHits}) or next;  
#        ($$maps{$pairID}{read1}{nHits} == 1 or $$maps{$pairID}{read2}{nHits} == 1) or next;   
#        printFixedMaps_Pair($maps, $pairID, 'read1', $mapFileH1, $toMap);
#        printFixedMaps_Pair($maps, $pairID, 'read2', $mapFileH2, $toMap);
#    }
#    close $mapFileH1;
#    close $mapFileH2;
#}

#sub printFixedMaps_Pair{
#    my ($maps, $pairID, $read, $mapFileH, $toMap) = @_;
#    if($$maps{$pairID}{$read}{nHits} <= $param{maxHits}){
#        foreach my $map (@{$$maps{$pairID}{$read}{maps}}){ print $mapFileH $map }
#    } else {
#        $$toMap{$pairID}{$read}++;
#    }   
#}

#sub printFixedToMap{
#    my ($toMap, $inReadFiles, $outReadFiles) = @_;
#    status("printing reads to remap...\n");
#    foreach my $read('read1','read2'){
#        open my $inFileH, "<", $$inReadFiles{$read};    
#        open my $outFileH, ">", $$outReadFiles{$read};  
#        while (my $name = <$inFileH>){
#            $name =~ m/^>(.*)\n/ or next;            
#            my $sequence = <$inFileH>;        
#            $$toMap{$1}{$read} and print $outFileH "$name$sequence";
#        }
#        close $inFileH;
#        close $outFileH;
#    }
#}

#sub mapFixedPartners{
#    my ($readFiles, $mapFiles, $read) = @_;
#    status("remapping promiscuous partner reads...\n");
#    my $nDisc = 3;
#    unlink $$mapFiles{$read};  #clear any previous mapping output
#    my $bowtieCommand = "$param{bowtiePath}/bowtie "; 
#    $bowtieCommand .= "-f -B 1 ";  #fasta, report refSeq 1-referenced
#    $bowtieCommand .= " -v $nDisc "; 
#    $bowtieCommand .= " -a ";  #return all hits  
#    $bowtieCommand .= " $param{refSeq} ";
#    $bowtieCommand .= " $$readFiles{$read} ";
#    $bowtieCommand .= " $$mapFiles{$read} "; 
#    $bowtieCommand .= " > $$mapFiles{$read}.stdout";
#    runBowtie($bowtieCommand);        
#}



#my @repLen = (1,2,3,4);
#my %fAllowed = (1=>0.75,2=>0.65,3=>0.55,4=>0.45);
#my %repeats = (1=>['A','C','G','T']);
#foreach my $repLen(2..4){
#    foreach my $base(@{$repeats{1}}){
#        foreach my $ext(@{$repeats{$repLen-1}}){
#            push @{$repeats{$repLen}}, "$base$ext"
#        }    
#    }
#}

#my $file = "C:/VAMP/090a_fx_1.fa.p";
#open my $fileH, "<", $file;
#my $counter = 0;
#my %maxRepeats;
#while(my $line = <$fileH>){
#    $counter > 10000 and last;
#    $line =~ m/^>/ or next;
#    my $read = <$fileH>;
#    
#    
#    $read = '';
#    foreach my $base(1..36){$read .= $repeats{1}[int(rand(4))]}
#    
#    getMaxRepeats($read, \%maxRepeats);
#    $counter++;
#}
#close $fileH;


##print "result = ".checkForRepeats($read)."\n";
##getMaxRepeats($read, \my %maxRepeats);
#foreach my $repLen(sort {$a <=> $b} @repLen){
#    foreach my $f (sort {$a <=> $b} keys %{$maxRepeats{$repLen}}){
#        print "$repLen\t$f\t${$maxRepeats{$repLen}}{$f}\n"
#    }
#}

#sub checkForRepeats{
#    my ($read) = @_;
#    my $readLen = length($read);
#    foreach my $repLen(@repLen){    
#        my $scalar = $repLen / $readLen;
#        my $fAllowed = $fAllowed{$repLen};
#        foreach my $repeat(@{$repeats{$repLen}}){
#            scalar(my @count = $read =~ m/$repeat/g) * $scalar > $fAllowed and return $repeat
#        }
#    }
#    return undef;
#}

#sub getMaxRepeats{
#    my ($read, $maxRepeats) = @_;
#    my $readLen = length($read);
#    foreach my $repLen(@repLen){    
#        my $scalar = $repLen / $readLen;
#        my $maxF = 0;
#        foreach my $repeat(@{$repeats{$repLen}}){
#            my $f = scalar(my @count = $read =~ m/$repeat/g) * $scalar;
#            $f > $maxF and $maxF = $f;
#        }
#        $$maxRepeats{$repLen}{(int($maxF/0.05))*0.05}++;
#    }
#}



##sub getFixedPairIDs{
##    my ($inMapFiles, $normPairIDs, $fixedPairIDs) = @_;
##    status("getting fixed pair IDs...\n");
##    foreach my $read('read1','read2'){
##        open my $inFileH, "<", $$inMapFiles{$read};
##        while (my $line = <$inFileH>){         
##            $line =~ m/Name=(.+?);/ or next;
##            my $pairID = $1;
##            $$normPairIDs{$pairID} and next;
##            $line =~ m/Hits=(.*?);/ or next;
##            my $hits = $1;
##            $hits == 1 or next;  #todo allow hits > 1 if can be stratified           
##            $$fixedPairIDs{$read}{$pairID}++;
##        }
##        close $inFileH;
##    }
##}

##sub prepareFixedPairs{
##    my ($inReadFiles, $outReadFiles, $fixedPairIDs) = @_;
##    status("extracting fixed pairs...\n");
##    foreach my $fixedRead_(1,2){    
##        my $partnerRead = ($fixedRead_ % 2) + 1;
##        my $fixedRead  = "read$fixedRead_";
##        $partnerRead  = "read$partnerRead";
##        open my $inFileH, "<", $$inReadFiles{$fixedRead};        
##        open my $outFileHF, ">>", $$outReadFiles{$fixedRead};
##        open my $outFileHP, ">>", $$outReadFiles{$partnerRead};
##        while (my $name = <$inFileH>){
##            $name =~ m/^>(.*)\n/ or next;   
##            my $pairID = $1;         
##            my $sequence = <$inFileH>;
##            $$fixedPairIDs{$fixedRead}{$pairID} and print $outFileHF "$name$sequence";
##            $$fixedPairIDs{$partnerRead}{$pairID} and print $outFileHP "$name$sequence";
##        }
##        close $inFileH;
##        close $outFileHF;
##        close $outFileHP;  
##    }   
##}

##sub extractFixed{
##    my ($sample) = @_;
##    $sample =~ m/.*_fx$/ or $sample = "$sample\_fx";    
##    status("extracting pairs...\n");
##    	getMapFiles($sample, \my %mapFiles);
##    	my $pairsTable = newTable('Pairs', $sample);    
##    
##        
##    


###      	runExtraction($mapFiles{read1}, $mapFiles{read2}, $pairsTable); 
### should just be able to set param rooted to prevent inversion of reads
### it is map reads that enforces maxHits, address this to recover L1 reads, etc.


##}

##    my (%normPairIDs, %fixedPairIDs);  
##    
##       
##    my $inFragsTable = getTableName('Frags', $inSample);    
##    getNormalPairIDs($inFragsTable, \%normPairIDs);   
##    getMapFiles($inSample, \my %inMapFiles);
##    getFixedPairIDs(\%inMapFiles, \%normPairIDs, \%fixedPairIDs);        
##    getReadFiles($inSample, 'purgedFasta', \my %inReadFiles);
##    getDirectories($outSample, \my %outDirs);
##    -d $outDirs{sample} or mkdir $outDirs{sample};
##    -d $outDirs{read1} or mkdir $outDirs{read1};
##    -d $outDirs{read2} or mkdir $outDirs{read2};        
##    getReadFiles($outSample, 'purgedFasta', \my %outReadFiles);
##    unlink $outReadFiles{read1};
##    unlink $outReadFiles{read2};
##    prepareFixedPairs(\%inReadFiles, \%outReadFiles, \%fixedPairIDs); 



##sub printFixedPairReads{
##    my ($inMapFiles, $outMapFiles, $fixedPairIDs) = @_;
##    status("printing fixed pair reads...\n");
##    foreach my $read('read1','read2'){
##        open my $inFileH, "<", $$inMapFiles{$read};
##        open my $outFileH, ">", $$outMapFiles{$read};              
##        while (my $line = <$inFileH>){         
##            $line =~ m/Name=(.+?);/ or next;
##            my $pairID = $1;
##            $$fixedPairIDs{$pairID} or next;
##            print $outFileH $line;
##        }
##    }   
##}  


##sub findFixed{
##    my (@inSamples) = @_;

##    foreach my $inSample (@inSamples){
##        my $outSample = $inSample."_an";
##        
##        status("getting statistics...\n");    
##        my $statsTable = getTableName('Stats', $inSample);
##        getStatistics($statsTable, \my %stats);
##        $param{minFragsSet1} = $stats{minCoverageFMap} / 2; #/2 to account for mono-allelic events
##        $param{minFragsSet1} >= 2 or $param{minFragsSet1} = 2; 
##        $param{maxFragsSet1} = $stats{maxCoverageFMap};
##        ($param{minFragsSet2}, $param{maxFragsSet2}, $param{minFragsSet}, $param{maxFragsSet}) =
##            (999999, 0, $param{minFragsSet1}, $param{maxFragsSet1}); 
##            
##    status("finding sets...\n");
##        $param{fragsTable1} = getTableName('Frags', $outSample);
##        $param{strataFilter} = '';
##        runFind($outSample);

##    status("establishing mismatch strata of fragments in sets...\n");
##        establishFragmentStrata($param{fragsTable1}, 1);
##        
##    status("re-finding sets within lowest mismatch strata...\n");
##        $param{strataFilter} = ' AND NSETSFRAG = -1 ';
##        runFind($outSample);

##    status("loading sets...\n");
##        loadFSData(1); 
##        
##    status("calculating fraction of promiscuous pairs in sets...\n");
##        calculateFracBadFrags(1); 

##    status("filling set counts...\n");
##        fillSetCounts($param{fragsTable1}, 1); 
##        
##    }
##} 



##        
##        

###        
###        printFixedPairReads(\%inMapFiles, \%outMapFiles, \%fixedPairIDs);
###    	my $outPairsTable = newTable('Pairs', $outSample); 
###        status("extracting pairs...\n");       
###      	runExtraction($outMapFiles{read1}, $outMapFiles{read2}, $outPairsTable); 
###        my $statsTable = getTableName('Stats', $inSample);
###        getStatistics($statsTable, \my %stats);
###        my $outFragsTable = newTable('Frags', $outSample);
###        status("parsing fragments...\n");
###        callExecuteParse($outPairsTable, $outFragsTable, \%stats);
###        purgeDuplicateFragments($outFragsTable);



1;

