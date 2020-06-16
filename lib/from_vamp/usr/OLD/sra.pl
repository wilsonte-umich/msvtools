#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs %bowtieFields));

sub prepareSRA{

    my @runs = (
#        {sample => '090',
#         library => '090a',
#         run => '090a'
#        },
#        {sample => '090',
#         library => '090b',
#         run => '090b'
#        },
#        {sample => 'D2',
#         library => 'D2a',
#         run => 'D2a'
#        },
#        {sample => 'D2',
#         library => 'D2b',
#         run => 'D2b'
#        },
#        {sample => 'A1A1',
#         library => 'A1A1a',
#         run => 'A1A1a'
#        },
#        {sample => 'A1A1',
#         library => 'A1A1b',
#         run => 'A1A1b'
#        },
#        {sample => 'A3A2',
#         library => 'A3A2ab',
#         run => 'A3A2a'
#        },
        {sample => 'A3A2',
         library => 'A3A2ab',
         run => 'A3A2b'
        },
        {sample => 'A3A2',
         library => 'A3A2c',
         run => 'A3A2c'
        }
    );

    foreach my $runRef (@runs) {
        print "$$runRef{run}\n";
        print "  archiving...\n";  
        
        my $outDir = "/home/wilsonte/sra/$$runRef{sample}/$$runRef{library}/$$runRef{run}";
        -d $outDir or system("mkdir -p $outDir");
        -d $outDir or die "could not find or create $outDir";   
        
        my @inFiles;
        foreach my $read (1..2){
            my $inFile = "$$runRef{run}_$read\_sequence.txt";
            my $inPath = "$param{inputPath}/$$runRef{run}/$read/$inFile";
            -e $inPath or die "could not find $inPath";
            push @inFiles, $inPath;
        }

        my $tarFile = "$outDir/$$runRef{run}.tar";
        system("tar -cvf $tarFile @inFiles");

        print "  calculating md5sum...\n";
        my $md5sum = qx/md5sum $tarFile/;
        print "    md5sum = $md5sum"; 

        print "  compressing...\n";
        system("gzip $tarFile");        
    }
}


1;



