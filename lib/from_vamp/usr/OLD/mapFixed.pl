#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %bowtieFields));

##########################
##need vamp.pl to assemble the following job threads
#mapFixed, one thread for each read
#prepareUnfixed, one thread total
#mapUnfixed, one thread for each read
#mergeFixedMaps (can just cascade from mapUnfixed)
##########################

sub mapFixed{
    my ($sample, $read) = @_;
    checkFixedOptions();
    getFixedReadFiles($sample, \my %readFiles);
    getFixedMapFiles($sample, \my %mapFiles);
    unlink ($mapFiles{fx}{$read}, $readFiles{fx}{$read}, $readFiles{ufx}{$read});  #clear any previous mapping output      
    my $bowtieCommand = "$param{bowtiePath}/bowtie "; 
    $bowtieCommand .= " -f -B 1 ";  #fasta, report refSeq 1-referenced
    $bowtieCommand .= " --max $readFiles{ufx}{$read} --al $readFiles{fx}{$read} ";      
    $bowtieCommand .= " -v $param{maxDisc} "; 
    $bowtieCommand .= " -k 1 -m 1 --best ";  #the one best hit, suppress all if >1 best hit    
    $bowtieCommand .= " $param{refSeq} "; 
    $bowtieCommand .= " $readFiles{$read} ";
    $bowtieCommand .= " $mapFiles{fx}{$read} ";   
    $bowtieCommand .= " > $mapFiles{fx}{$read}.stdout";
    runBowtie($bowtieCommand);     
}

sub prepareUnfixed{
    my ($sample) = @_;
    getFixedReadFiles($sample, \my %readFiles);
    getFixedMapFiles($sample, \my %mapFiles);
    foreach my $fixedRead_(1,2){    
        my $unfixedRead = ($fixedRead_ % 2) + 1;
        my $fixedRead  = "read$fixedRead_";
        $unfixedRead  = "read$unfixedRead";
        my %fixedReads;
        open my $fxReadsFileH, "<", $readFiles{fx}{$fixedRead};
        while(my $nameLine = <$fxReadsFileH>){ $nameLine =~ m/^>(.*)\n/ and $fixedReads{$1}++ }
        close $fxReadsFileH;
        open my $ufxReadsFileH, "<", $readFiles{ufx}{$unfixedRead};        
        my $ufxTmpFile = $readFiles{ufx}{$unfixedRead}.".tmp";
        open my $ufxTmpFileH, ">", $ufxTmpFile;
        while(my $nameLine = <$ufxReadsFileH>){         
            $nameLine =~ m/^>(.*)\n/ or next;
            if($fixedReads{$1}){
                my $readLine = <$ufxReadsFileH>;
                #unless(checkSimpleRepeats($readLine)){ print $ufxTmpFileH "$nameLine$readLine" } 
                print $ufxTmpFileH "$nameLine$readLine"
            } 
        }
        close $ufxReadsFileH;
        close $ufxTmpFileH;
        system("mv $readFiles{ufx}{$unfixedRead} $readFiles{ufx}{$unfixedRead}.bu");
        system("mv $ufxTmpFile $readFiles{ufx}{$unfixedRead}");  
    }    
}

sub mapUnfixed{
    my ($sample, $read) = @_;
    checkFixedOptions();
    getFixedReadFiles($sample, \my %readFiles);
    getFixedMapFiles($sample, \my %mapFiles);    
    unlink $mapFiles{ufx}{$read};  #clear any previous mapping output 
    my $maxHits = $param{maxHits} + 1;
    my $bowtieCommand = "$param{bowtiePath}/bowtie "; 
    $bowtieCommand .= " -f -B 1 ";  #fasta, report refSeq 1-referenced
    $bowtieCommand .= " -v $param{maxDisc} ";     
    $bowtieCommand .= " -k $maxHits "; #return up to maxHits + 1; suppress only those maps OVER maxHits + 1
    $bowtieCommand .= " $param{refSeq} ";     
    $bowtieCommand .= " $readFiles{ufx}{$read} ";
    $bowtieCommand .= " $mapFiles{ufx}{$read} ";   
    $bowtieCommand .= " > $mapFiles{ufx}{$read}.stdout";
    runBowtie($bowtieCommand);     
}

#this is problematic and unnecessary
#extractPairs ought to be able to handle an array of input files
sub mergeFixedMaps{
    my ($inSample, $outSample) = @_;
    checkFixedOptions();
    getFixedMapFiles($inSample, \my %inMapFiles);   
    getFixedMapFiles($outSample, \my %outMapFiles); 
    getDirectories($outSample, \my %outDirs);      
    -d $outDirs{sample} or mkdir $outDirs{sample};  
    foreach my $read('read1', 'read2'){    
        -d $outDirs{$read} or mkdir $outDirs{$read};   
        system("cat $inMapFiles{fx}{$read} $inMapFiles{ufx}{$read} > $outMapFiles{$read}");
        #unlink ($inMapFiles{fx}{$read}, $inMapFiles{ufx}{$read});  
    }
}

sub checkFixedOptions{
    unless ($param{mapType} eq 'bowtie'){
        status("mapFixed requires Bowtie, changing mapType to bowtie...\n");
        $param{mapType} = 'bowtie'; 
    } 
    if ($param{maxDisc} > 3){
        status("Bowtie allows up to 3 mismatches, changing maxDisc to 3...\n");
        $param{maxDisc} = 3;
    }     
}

sub getFixedReadFiles{
    my ($sample, $readFiles) = @_;
    getReadFiles($sample, 'purgedFasta', $readFiles); 
    $$readFiles{fx}{read1} = $$readFiles{read1}.".fx";
    $$readFiles{fx}{read2} = $$readFiles{read2}.".fx";
    $$readFiles{ufx}{read1} = $$readFiles{read1}.".ufx";
    $$readFiles{ufx}{read2} = $$readFiles{read2}.".ufx";
}

sub getFixedMapFiles{
    my ($sample, $mapFiles) = @_;
    getMapFiles($sample, $mapFiles); 
    $$mapFiles{fx}{read1} = $$mapFiles{read1}.".fx";
    $$mapFiles{fx}{read2} = $$mapFiles{read2}.".fx";
    $$mapFiles{ufx}{read1} = $$mapFiles{read1}.".ufx";
    $$mapFiles{ufx}{read2} = $$mapFiles{read2}.".ufx";
}

#sub countUnfixedMaps{
#    my $file = '/home/wilsonte/vamp/data/090a/1/090a_1.bowtie.ufx';
#    my (%pairIDs, $counter);
#    open my $fileH, "<", $file;
#    while (my $line = <$fileH>){
#        chomp $line;
#        my @in = getFields($line);
#        my $pairID = $in[$bowtieFields{name}];  
#        $pairIDs{$pairID} or $counter++;        
#        $pairIDs{$pairID}++;
#        $counter >= 10000 and last;
#    }
#    close $fileH;    
#    my %counts;
#    foreach my $pairID(keys %pairIDs){
#        my $count = $pairIDs{$pairID};
#        $count > 100 and $count = 101;
#        $counts{$count}++;
#    }
#    open my $outH, ">", "$param{vampPath}/ufxHitCounts.xls";
#    foreach my $count(sort {$a <=> $b} keys %counts){print $outH "$count\t$counts{$count}\n" }
#    close $outH;
#}

    #my $outSample = $inSample."_fx";    
    #my $inFragsTable = getTableName('Frags', $inSample); 
#      
#     
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

1;



