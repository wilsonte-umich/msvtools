#!/usr/bin/perl -w
use strict;
use warnings;

#########################################################################
#mapReadsSplit.pl extends mapReads.pl by providing the ability to 
#split the input reads into small changes which are submitted as
#multiple parallel jobs, to maximize server use and analysis speeed.
#########################################################################

use vars(qw(%param %types %fields %refSeqs $paramCat));

my %dirs;
my %fastqFiles;
my %qseqFiles;
my %fastaFiles;
my %gffFiles;

sub prepareReadsSplit{
    my ($sample) = @_;
    my $readType = createFastaFile($sample);
    if ($readType eq 'fasta') {  #true if readType was fasta or fastq was converted to fasta
        status("purging duplicate reads...\n");
            purgeDuplicatesSplit($sample);
            $readType = 'purgedFasta';
    } else {
        die "'$readType' is not a recognized readType"
    }    
}

sub purgeDuplicatesSplit{  #reject identical read pairs, within limitation that must be 100% identical
    my ($sample) = @_;
    getReadFiles($sample, 'fasta', \%fastaFiles);
    getReadFiles($sample, 'purgedFasta', \my%purgedFiles);
    getDirectories($sample, \%dirs);
    my ($inCounter, $outCounter, $countsRef) = purgeDuplicatesPairedReadsSplit(\%purgedFiles);
    status("\n  $inCounter\tread pairs processed, yielding...\n");
    status("  $outCounter\tunique read pairs\n");
    my $countFile = "$dirs{sample}/pairCounts.xls";
    status("  writing pair count frequencies to\n    $countFile\n");
        open my $countFileH, ">", $countFile or die "could not open $countFile";
        foreach my $count (sort {$a <=> $b} keys %$countsRef){print $countFileH "$count\t$$countsRef{$count}\n"}
        close $countFileH;           
}

sub purgeDuplicatesPairedReadsSplit{  #reject identical read pairs, within limitation that must be 100% identical  
    my ($purgedFilesRef) = @_;
    status("  collecting unique read pairs...\n");
        open my $fasta1H, "<", $fastaFiles{read1} or die "could not open $fastaFiles{read1}";
        open my $fasta2H, "<", $fastaFiles{read2} or die "could not open $fastaFiles{read2}";
        my $read; 
        my %pairs = ();
        my $inCounter = 0;
        while (<$fasta1H>){
            my $line1 = $_;
            my $line2 = <$fasta2H>;
            if ($line1 =~ m/^>/){
                my ($read1) = split("/", $line1);
                my ($read2) = split("/", $line2);
                ($read1 eq $read2) or die "corresponding read order in fasta files violated!";
                $read = $read1;
                $inCounter++;
            } elsif ($read) {
                my ($seq1, $seq2) = sort {$a cmp $b} ($line1, $line2);  #order sequences to account for reversed duplicates
                $pairs{"$seq1:$seq2"}++;  #hold unique read pairs in hash and count the number of occurrences
                $read = 0
            }
        }
        close $fasta1H;
        close $fasta2H;
    
    -d $dirs{prepared_read1} or mkdir $dirs{prepared_read1};
    -d $dirs{prepared_read2} or mkdir $dirs{prepared_read2};
    
    status("  writing new purged fasta files\n");
        my $pairID = 0;
        my %counts = ();
        foreach my $pair (keys %pairs){
            $pairID++;
            my $pairIDGroup = int($pairID / $param{mapChunk});
            my ($seq1, $seq2) = split (":", $pair); 
            my $fileH = getPurgedFileH('read1', $pairIDGroup);
            print $fileH ">$pairID\n";  #lose old pair names, just use an incrementing numerical pairID
            print $fileH "$seq1";  #newline is propogated from original file
            $fileH = getPurgedFileH('read2', $pairIDGroup);
            print $fileH ">$pairID\n";
            print $fileH "$seq2";
            $counts{$pairs{$pair}}++;
        }  
        closePurgedFileHs();
        
    return ($inCounter, $pairID, \%counts);      
}

my %purgedFileHs;
sub getPurgedFileH{
    my ($read, $pairIDGroup) = @_;
    $purgedFileHs{$read}{$pairIDGroup} and return $purgedFileHs{$read}{$pairIDGroup};
    my $prepared = "prepared_$read";
    my $dir = "$dirs{$prepared}/$pairIDGroup";
    -d $dir or mkdir $dir;
    my $file = "$dir/$pairIDGroup.fa";
    open my $fileH, ">", $file or die "could not open $file\n";
    $purgedFileHs{$read}{$pairIDGroup} = $fileH;
    return $fileH; 
}

sub closePurgedFileHs{
    foreach my $read('read1','read2'){
        foreach my $pairIDGroup(keys %{$purgedFileHs{$read}}){
            my $fileH = $purgedFileHs{$read}{$pairIDGroup};
            close $fileH;
        }
    }
}

#############
my %wallTimes = (pass=>'24:00:00',bowtie=>'1:30:00');
my %mems = (pass=>'20000',bowtie=>'4000');
#############
sub mapReadsSplit{
    my ($inputSample) = @_; 
    getDirectories($inputSample, \%dirs);
    my $perlTarget = "$param{vampPath}/bin/ExecuteInstruction.pl";
    foreach my $read ('read1','read2') {
        my $prepared = "prepared_$read";
        my $prepDir = $dirs{$prepared};
        foreach my $groupDir (<$prepDir/*>){
            -d $groupDir or next;
            foreach my $faFile (<$groupDir/*.fa>){
                my $vampCommand = "perl $perlTarget $paramCat mapReadsFile $faFile";
                my $jobFile = getJobFile($faFile);
                open my $jobFileH , ">", $jobFile;
                print $jobFileH ":\n";  #required first line
                print $jobFileH "#PBS -d $groupDir\n";  #where the job runs from and log files are deposited
                print $jobFileH "#PBS -j oe\n";  #merge STDERR and STDOUT into single log file
                print $jobFileH "#PBS -l walltime=$wallTimes{$param{mapType}},mem=$mems{$param{mapType}}"."M\n";  #CPU usage limits to assist with job queueing 
                print $jobFileH "$vampCommand\n";  #send the propagated information to ExecuteInstruction.pl
                close $jobFileH;   

                #queue it and return the resulting job ID
                my $jobID = qx/qsub $jobFile/;
                chomp $jobID;      
            }
        }
    }
}

sub getJobFile{
    my ($fastaFile) = @_;
    my $i = 1;
    my $jobFile;
    while (!$jobFile or -e $jobFile){
        $jobFile = "$fastaFile\_mapReads_.$i";
        $i++;  
    }
    return $jobFile;
}

sub mapReadsFile{  #not called directly, called by mapReadsSplit...
    my ($faFile) = @_;
    my $paramFile = "$faFile\_parameters.txt";     
    printParameters($paramFile, 'mapReadsFile', $faFile);     
    if ($param{mapType} eq 'pass'){
        mapReadsFile_PASS($faFile)
    } elsif ($param{mapType} eq 'bowtie') {
        mapReadsFile_Bowtie($faFile)
    } else {
        die "$param{mapType} is not a recognized mapType";
    }
}

sub mapReadsFile_Bowtie {
    my ($faFile) = @_;
    (-e $faFile) or die "$faFile does not exist";    
    my $inputFile = $faFile; 
    
    open my $statusFileH, ">", "progress.txt";
    reportStatus($statusFileH, "mapReads_Bowtie begun...");    
    
    if ($param{maxDisc} > 3){
        status("Bowtie allows up to 3 mismatches, changing maxDisc to 3...\n");
        $param{maxDisc} = 3;
    }

    #parameters which are held constant   
    my $outputFile = "$inputFile.bowtie";  #final catenated data destinations
    my $notAligned = "$outputFile.not_aligned";
    my $tooManyHits = "$outputFile.too_many_hits";
    unlink ($outputFile, $notAligned, $tooManyHits);  #clear any previous mapping output

    #loop through possible discrepancy thresholds in DESCENDING order
    my (@outputFiles, @notAlignedFiles, @tooManyHitsFiles);
    my $tooManyHits_nDisc;
    for (my $nDisc = $param{maxDisc}; $nDisc >= 0; $nDisc--){
        
        reportStatus($statusFileH, "  nDisc = $nDisc");
        
        my $notAligned_nDisc = "$notAligned.$nDisc";
        $tooManyHits_nDisc = "$tooManyHits.$nDisc";        
        my $output_nDisc = "$outputFile.$nDisc";
        
        push @notAlignedFiles, $notAligned_nDisc;
        push @tooManyHitsFiles, $tooManyHits_nDisc;        
        push @outputFiles, $output_nDisc;

        #calculate pass parameters dependent on discrepancy threshold and run pass
        print "-" x 60; print "\n$nDisc mismatches\n"; print "-" x 60;
        my $bowtieCommand = "$param{bowtiePath}/bowtie "; 
        $bowtieCommand .= "-f -B 1 ";  #fasta, report refSeq 1-referenced
        $bowtieCommand .= "--un $notAligned_nDisc ";   
        $bowtieCommand .= "--max $tooManyHits_nDisc ";
        $bowtieCommand .= " -v $nDisc ";
        $bowtieCommand .= " -k $param{maxHits} ";        
        $bowtieCommand .= " -m $param{maxHits} ";
        $bowtieCommand .= " $param{refSeq} ";
        $bowtieCommand .= " $inputFile ";
        $bowtieCommand .= " $output_nDisc "; 
        $bowtieCommand .= " > $outputFile.stdout";
        runBowtie($bowtieCommand);        
        
        reportStatus($statusFileH, "    Bowtie done");
        
        $inputFile = $tooManyHits_nDisc;
    }
    
    system("cat @notAlignedFiles > $notAligned");
    unlink (@notAlignedFiles);  
    system("cat $tooManyHits_nDisc > $tooManyHits");
    unlink (@tooManyHitsFiles);    
    system("cat @outputFiles > $outputFile");
    unlink (@outputFiles);      
    
    reportStatus($statusFileH, "mapReads done\n");
    close $statusFileH;
}

1;

