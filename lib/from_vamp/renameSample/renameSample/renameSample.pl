#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %command %types %fields %refSeqs %reverseRefSeqs %gffFields %bowtieFields));

#callable parameters and commands in this script
#defined $param{xxx} or $param{xxx} = xxx; 
defined $param{renameTestOnly} or $param{renameTestOnly} = 0; 
defined $command{renameSample} or $command{renameSample} = ['singleThread', '1:00:00', 1000, 0];

#renames ALL files and database tables from one sample name to another
sub renameSample {
    my ($sampleOut, $sampleIn) = @_;
    ($sampleOut and $sampleIn) or die "usage:  renameSample sampleOut sampleIn\n";
    my $startPath = "$param{inputPath}/$sampleIn";
    -d $startPath and renameFiles($startPath, $param{inputPath}, $sampleOut, $sampleIn);
    renameTables($sampleOut, $sampleIn);  
}
sub renameFiles {
   my ($inFile, $filePath, $sampleOut, $sampleIn) = @_; 
   if ($inFile =~ m/^$filePath\/$sampleIn(.*)/){  
       my $outFile = "$filePath/$sampleOut$1";
       my $mvCommand = "mv $inFile $outFile";
       status("$mvCommand\n");
       unless($param{renameTestOnly}){
           system($mvCommand);
           $inFile = $outFile;          
       }
   }
   if (-d $inFile) {
       foreach my $file(<$inFile/*>){
           renameFiles($file, $inFile, $sampleOut, $sampleIn);
       } 
   }
}
sub renameTables {
    my ($sampleOut, $sampleIn) = @_;
    $sampleIn = "\U$sampleIn";
    $sampleOut = "\U$sampleOut";
    runSQL("SELECT table_name FROM user_tables ORDER BY table_name", \my($tableIn));
    my @tables;
    while (fetchRow()){ push @tables, $tableIn } 
    foreach my $tableIn (@tables){
        if($tableIn =~ m/(.+)_$sampleIn(.*)/){
            my $tableOut = "$1_$sampleOut$2";
            my $renameSQL = "ALTER TABLE $tableIn RENAME TO $tableOut";
            status("$renameSQL\n");
            unless($param{renameTestOnly}){
               runSQL($renameSQL);   
            } 
        }
    }
}

1;
