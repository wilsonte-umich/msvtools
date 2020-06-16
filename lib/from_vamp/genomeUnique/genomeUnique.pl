#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs $paramCat));

#callable parameters and commands in this script
#binSize is required but already defined by VAMP
#defined $param{upstreamBins} or $param{upstreamBins} = 5; #number of bins preceding the gene start
defined $command{parseGenomeUnique} or $command{parseGenomeUnique} = ['singleThread', '2:00:00', 1000, 0];
defined $command{loadGenomeUnique} or $command{loadGenomeUnique} = ['singleThread', '2:00:00', 1000, 0];
defined $command{clearGenomeUniqueFiles} or $command{clearGenomeUniqueFiles} = ['singleThread', '4:00:00', 1000, 0];
defined $command{parseBlocksGenomeUnique} or $command{parseBlocksGenomeUnique} = ['singleThread', '24:00:00', 2000, 0];

#defined $command{fixUnqTable} or $command{fixUnqTable} = ['singleThread', '24:00:00', 2000, 0];


my %wallTimes = (parse => '8:00:00', load => '24:00:00');
my %mbRams = (parse => '3000', load => '10000');
my (%segs, %counts, $maxBufferSegs);
my $rootLength = 6;
my ($segDir, $jobsDir, $logsDir, $chromDir);

sub parseGenomeUnique {
    getPathsGenomeUnique();
    status("queueing parse jobs by chromosome...\n"); 
    foreach my $chrom(1..nChrom()){
        status("  $chrom\n");
        queueJobGenomeUnique('parseSegsGenomeUnique', $chrom, $wallTimes{parse}, $mbRams{parse});
    } 
}
sub loadGenomeUnique { #only run this after all parseGenomeUnique jobs are finished!
    getPathsGenomeUnique();
    my $unqTable = newTable('Unq', "$param{refSeq}_$param{readLength}");
    status("queueing load jobs...\n"); 
    my @bases = ('A','C','G','T');
    foreach my $base1(@bases){
        foreach my $base2(@bases){
            my $baseKey = "$base1$base2";
            status("  $baseKey\n");
            queueJobGenomeUnique('loadSegsGenomeUnique', $baseKey, $wallTimes{load}, $mbRams{load}); 
        }
    }
}
sub clearGenomeUniqueFiles {
    foreach my $chrom(1..nChrom()){
        getPathsGenomeUnique($chrom);
        system("rm $chromDir/*"); #clear any previous results
    } 
}

sub queueJobGenomeUnique {
    my ($command, $param, $wallTime, $mbRam) = @_;
    my $perlTarget = "$param{vampPath}/bin/ExecuteInstruction.pl";
    my $vampCommand = "perl $perlTarget $paramCat $command $param";
    my $jobFile = "$jobsDir/$command\_$param";
    open my $jobFileH , ">", $jobFile;
    print $jobFileH ":\n";  #required first line
    print $jobFileH "#PBS -d $logsDir\n";  #where the job runs from and log files are deposited
    print $jobFileH "#PBS -j oe\n";  #merge STDERR and STDOUT into single log file
    print $jobFileH "#PBS -l walltime=$wallTime,mem=$mbRam"."M\n";  #CPU usage limits to assist with job queueing 
    print $jobFileH "$vampCommand\n";  #send the propagated information to ExecuteInstruction.pl
    close $jobFileH;   
    my $jobID = qx/qsub $jobFile/; 
}

sub getPathsGenomeUnique {
    my ($chrom) = @_;
    $segDir = "$param{refSeqPath}/$param{refSeq}";
    -d $segDir or mkdir $segDir;
    $segDir .= '/segs';
    -d $segDir or mkdir $segDir;
    $jobsDir = "$segDir/jobs";
    -d $jobsDir or mkdir $jobsDir;
    $logsDir = "$segDir/logs";
    -d $logsDir or mkdir $logsDir;   
    if ($chrom) {
        $chromDir = "$segDir/$chrom";
        -d $chromDir or mkdir $chromDir;   
    }   
}

sub parseSegsGenomeUnique {
    my ($chrom) = @_;
    getMaxBufferSegs();
    getPathsGenomeUnique($chrom);
    system("rm $chromDir/*"); #clear any previous results
    status("loading $param{refSeq} chromosome $chrom...\n"); 
    my $lineSize = getLineSize(); 
    my $chromLines = loadChromLines(undef, $chrom);     
    my $increment = $param{readLength} - 1;
    status("writing genome x-mers to tmp files...\n"); #subs found in findHomozygous.pl
    my $pos = 1;
    while(my $seg = getDNASegment($chrom, $pos, $pos + $increment, $lineSize, $chromLines)){
        $seg = "\U$seg"; #force uppercase 
        unless($seg  =~ m/N/) { #ignore any segment not purely ACGT
            ($seg) = sort {$a cmp $b} ($seg, reverseComplement($seg)); #order segment and its RC to allow dups to RC
            $seg =~ m/^(.{$rootLength})(.*)/;
            my ($binSeq, $seg) = ($1, $2);
            commitSegGenomeUnique($binSeq, $seg, $chrom, $pos); 
        }
        $pos++;
    }
    status("  maxPos $pos\n");
    commitAllBuffersGenomeUnique();;
}

sub commitSegGenomeUnique {
    my ($binSeq, $seg, $chrom, $pos) = @_;
    $segs{$binSeq} .= "$seg,$chrom,$pos\n";
    $counts{$binSeq}++;
    $counts{$binSeq} > $maxBufferSegs and commitBufferGenomeUnique($binSeq);
}
sub commitAllBuffersGenomeUnique {
    foreach my $binSeq(keys %segs){
         commitBufferGenomeUnique($binSeq);
    }
}
sub commitBufferGenomeUnique {
    my ($binSeq) = @_;
    my $binFile = "$chromDir/$binSeq.csv";
    open my $binFileH, ">>", $binFile or die "could not open $binFile: $!\n";
    $segs{$binSeq} and print $binFileH $segs{$binSeq};
    close $binFileH;
    $segs{$binSeq} = undef;
    $counts{$binSeq} = 0;
}

sub getMaxBufferSegs {
    my $bytesRam = 1E9; #allow up to 1 Gb RAM usage for all buffers
    my $nBuffers = 4 ** $rootLength;
    my $bytesPerSeq = $param{readLength} - $rootLength +    #the segment
                      length(nChrom()) +                    #the longest chromosome
                      length(1E8) +                         #the longest position
                      3;                                    #two commas and newline
    $maxBufferSegs = int($bytesRam / ($bytesPerSeq * $nBuffers));
    status("maxBufferSegs = $maxBufferSegs\n");
}

sub loadSegsGenomeUnique {
    my ($baseKey) = @_;
    getPathsGenomeUnique();
    my $unqTable = getTableName('Unq', "$param{refSeq}_$param{readLength}");    
    foreach my $file(<$segDir/1/*$baseKey.csv>){
        $file =~ m/$segDir\/1\/(.*)$baseKey\.csv/;
        my $binSeq = "$1$baseKey";
        status("$binSeq\n");        
        my (%good, %bad);   
        foreach my $chrom(1..nChrom()){
            getPathsGenomeUnique($chrom);
            my $binFile = "$chromDir/$binSeq.csv";
            -e $binFile or next;
            open my $binFileH, "<", $binFile or die "could not open $binFile\n";        
            while(my $seqLine = <$binFileH>){
                chomp $seqLine;
                my ($seg,$chrom,$pos) = split(",", $seqLine);
                $bad{$seg} and next; #third encounter of segment
                if($good{$seg}){ #second encounter of segment
                    $bad{$seg} = 1;
                    delete $good{$seg};
                } else { #first encoutnter of segment
                    $good{$seg} = [$chrom,$pos];
                }
            }    
            close $binFileH;  
            #unlink $binFile;            
        }   
        my $unqFile = "$unqTable\_$binSeq.csv";
        open my $unqFileH, ">", $unqFile or die "could not open $unqFile\n";        
        foreach my $seg(keys %good){
            my ($chrom, $pos) = @{$good{$seg}};
            my $bin25MB = int($pos/25E6);
            print $unqFileH "$chrom,$pos,$bin25MB\n";          
        } 
        close $unqFileH;
        loadData($unqFile, $unqTable, ",", "CHROMOSOME, POSITION, BIN25MB");             
    } 
}

sub parseBlocksGenomeUnique {
    my $unqTable = getTableName('Unq', "$param{refSeq}_$param{readLength}");
    my $unqBlockTable =newTable('UnqB', "$param{refSeq}_$param{readLength}");
    my $outFile = "$unqBlockTable.csv";
    open my $outFileH, ">", $outFile or die "could not open $outFile\n";
    foreach my $chrom(1..nChrom()){
        runSQL("SELECT position
                FROM $unqTable
                WHERE chromosome = $chrom
                ORDER BY position",
                \my$pos);
        fetchRow();
        $pos or next;
        $pos == 1 or commitBlockGenomeUnique($outFileH, $chrom, 1, $pos - 1, 0);  
        my ($prevPos, $goodStart) = ($pos, $pos);
        while (fetchRow()){
            unless($pos == $prevPos + 1){
                commitBlockGenomeUnique($outFileH, $chrom, $goodStart, $prevPos, 1);  
                commitBlockGenomeUnique($outFileH, $chrom, $prevPos + 1, $pos - 1, 0);  
                $goodStart = $pos;
            }
            $prevPos = $pos;
        }   
        commitBlockGenomeUnique($outFileH, $chrom, $goodStart, $prevPos, 1);       
    }  
    close $outFileH;
    loadData($outFile, $unqBlockTable, ",", "CHROMOSOME, START_, END_, UNIQUE");  
    
}
sub commitBlockGenomeUnique{
    my ($outFileH, $chrom, $blockStart, $blockEnd, $unique) = @_;
    print $outFileH "$chrom,$blockStart,$blockEnd,$unique\n";
}



#sub fixUnqTable {
#    my $unqTable = newTable('Unq', "$param{refSeq}_$param{readLength}");
#    my $inTable = "$unqTable\_tmp";
#    updateTable('Unq', $unqTable, "SELECT * FROM $inTable");
#    #dropTable($inTable);   
#}


1;


