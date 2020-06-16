#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %bowtieFields));

my %qseqColumns = (machine=>0,run=>1,lane=>2,
                    tile=>3,x=>4,y=>5,
                    index=>6,read=>7,
                    sequence=>8,quality=>9,filter=>10);
                    
sub prepareSRA{

    my @runs = (
        {sample => '090',
         library => '090a',
         run => '090a'
        },
        {sample => '090',
         library => '090b',
         run => '090b'
        },
        {sample => 'D2',
         library => 'D2a',
         run => 'D2a'
        },
        {sample => 'D2',
         library => 'D2b',
         run => 'D2b'
        },
        {sample => 'A1A1',
         library => 'A1A1a',
         run => 'A1A1a'
        },
        {sample => 'A1A1',
         library => 'A1A1b',
         run => 'A1A1b'
        },
        {sample => 'A3A2',
         library => 'A3A2ab',
         run => 'A3A2a'
        },
        {sample => 'A3A2',
         library => 'A3A2ab',
         run => 'A3A2b'
        },
        {sample => 'A3A2',
         library => 'A3A2c',
         run => 'A3A2c'
        }

#        {sample => 'vac6wt',
#         library => 'vac6wt',
#         run => 'vac6wt'
#        },
#        {sample => 'vac6mut',
#         library => 'vac6mut',
#         run => 'vac6mut'
#        },
#        {sample => 'vac22wt',
#         library => 'vac22wt',
#         run => 'vac22wt'
#        },
#        {sample => 'vac22mut',
#         library => 'vac22mut',
#         run => 'vac22mut'
#        }

    );
    
    my $outDir = "/home/wilsonte/sra/090_files";
    -d $outDir or system("mkdir -p $outDir");
    -d $outDir or die "could not find or create $outDir";      
    my $reportFile = "$outDir/md5sum_report.csv";
    open my $reportFileH, ">", $reportFile;

    foreach my $runRef (@runs) {
        print "$$runRef{run}\n";
        foreach my $read (1..2){
            my $inFile = "$$runRef{run}_$read\_sequence.txt";
            my $inPath = "$param{inputPath}/$$runRef{run}/$read/$inFile";
            my ($read, $readLength) = getReadLength($inPath);
            my $md5sum = calculateMD5Sum($inPath);
            copyAndCompressFile($inPath, $inFile, $outDir);
            my $report = "$inPath,$inFile,$read,$readLength,$md5sum\n";
            print $report;
            print $reportFileH $report;
#            push @inFiles, $inPath;      
        }
#        my $tarFile = "$outDir/$$runRef{run}.tar";
#        system("tar -cvf $tarFile @inFiles");      
    }
    close $reportFileH;
}

sub getReadLength{
    my ($inPath) = @_;
    print "  getting read length...\n";
    open my $inPathH, "<", $inPath or die "could not open $inPath";
    my $discard = <$inPathH>;
    my $read = <$inPathH>;
    chomp $read;
    $read =~s/\r//g;
    close $inPathH;
    return ($read, length($read));
}

sub calculateMD5Sum{
    my ($inPath) = @_;
    print "  calculating md5sum...\n";
    my $retVal = qx/md5sum $inPath/;
    my @retVals = split(" ", $retVal);
    #b4fb137af5fe2a3bb4d09877d466543e  /home/wilsonte/vamp/data/090a/1/090a_1_sequence.txt
    return $retVals[0];
}

sub copyAndCompressFile{
    my ($inPath, $inFile, $outDir) = @_;
    my $outPath = "$outDir/$inFile";
    print "  copying file...\n";
    system("cp $inPath $outPath");  
    print "  compressing file...\n";
    system("gzip $outPath");        
}

sub qseq2fastq {

    my @samples = ("vac6wt", "vac6mut");
    
    foreach my $sample (@samples) {
        foreach my $read (1..2) {
            my $fastqDir = "/home/wilsonte/vamp/data/$sample/$read";
            my $fastqFile = "$fastqDir/$sample\_$read\_sequence.txt";
            my $qseqDir = "$fastqDir/qseq";
            open my $fastqFileH, ">", $fastqFile or die "could not open $fastqFile";
            foreach my $qseqFile (<$qseqDir/*qseq.txt>){ #read files from qseq subdir
                open my $qseqFileH, "<", $qseqFile or die "could not open $qseqFile";
                while (<$qseqFileH>){
                    chomp $_;
                    my @in = split(/\t/, $_);                
                    if ($in[$qseqColumns{filter}]){                
                        my $readID = join(":", ($in[$qseqColumns{machine}], 
                                                @in[$qseqColumns{lane}..$qseqColumns{y}]) );
                        $readID .= "\#$in[$qseqColumns{index}]/$in[$qseqColumns{read}]";
                        print $fastqFileH '@'."$readID\n";
                        print $fastqFileH "$in[$qseqColumns{sequence}]\n";
                        print $fastqFileH '+'."$readID\n";
                        print $fastqFileH "$in[$qseqColumns{quality}]\n";     
                    }
                }
                close $qseqFileH;
            }
            close $fastqFileH;
        }
    }
}

1;



