#!/usr/bin/perl
use strict;
use warnings;

#########################################################################
#converts BED files to CNV or Array tables in Oracle
#########################################################################

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs));  #common parameters used by most analysis scripts

################################################
#CNV tables
#-----------------------------------------------
sub importCNVs{
    #expects input file = vampPath/CNVs/sample/sample_CNVs.refSeq
    #in format: name,type,chrom,start,end 
    #where chrom = UCSC chrom field, e.g. 'chr13'
    #and type = numerical VAMP Frag/Set type
    #no header line
    my ($sample, $refSeq, $liftTo) = @_;
    my $filePath = "$param{vampPath}/CNVs/$sample";
    my $bedFile = writeCNVsToBED($filePath, $sample, $refSeq);
    importCNVs_BED($bedFile, 0, $sample, $refSeq, $liftTo);    
}
#-----------------------------------------------
sub importCNVs_BED{
    #expects pre-existent BED file = vampPath/CNVs/sample/sample_CNVs.refSeq.BED
    #BED chrom = UCSC chrom field, e.g. 'chr13'
    #BED name = type::name where type = numerical VAMP Frag/Set type
    my ($bedFile, $stripHeader, $sample, $refSeq, $liftTo) = @_;
    $stripHeader and $bedFile = stripBEDHeader($bedFile);    
    importBEDFile_CNV($bedFile, $sample, $refSeq);
    if ($liftTo){
        $bedFile = liftOverBED($bedFile, $refSeq, $liftTo);
        importBEDFile_CNV($bedFile, $sample, $liftTo);
    }
}
#-----------------------------------------------
sub writeCNVsToBED{ 
    my ($filePath, $sample, $refSeq) = @_;
    my $CNVsFile = "$filePath/$sample\_CNVs.$refSeq";      
    print "$CNVsFile\n";    
    open my $CNVsFileH, "<", $CNVsFile or die "could not open $CNVsFile"; 
    my $bedFile = "$CNVsFile.BED";      
    open my $bedFileH, ">", $bedFile;
    print "  $bedFile\n";      
    while (<$CNVsFileH>){
        chomp $_;  
        my ($name, $type, $chrom, $start, $end) = split(",", $_);  
        $end =~ m/^(\d*)/;  #some file formats had some weird character at the end, this gets rid of it!!!
        $end = $1;
        print $bedFileH "$chrom $start $end $type\::$name\n";
    }
    close $CNVsFileH;                
    close $bedFileH;
    return $bedFile;
}
#-----------------------------------------------
sub importBEDFile_CNV{  
    my ($bedFile, $sample, $refSeq) = @_;    
    my $CNVsTable = newTable('CNVs', "$sample\_$refSeq");    
    my $CNVsFile = "$CNVsTable.csv";
    open my $bedFileH, "<", $bedFile or die "could not open $bedFile";     
    open my $CNVsFileH, ">", $CNVsFile;    
    while (<$bedFileH>){
        chomp $_;      
        my ($chr, $start, $end, $typeName) = split(" ", $_);    
        my $chrom = $refSeqs{$refSeq}{$chr};
        my ($type, $name) = split("::", $typeName);
        print $CNVsFileH "$name,$type,$chrom,$start,$end\n";
    }
    close $bedFileH;
    close $CNVsFileH; 
    loadData($CNVsFile, $CNVsTable, ',', "CNVNAME, CNVTYPE, CHROMOSOME, START_, END_")
}
################################################


#################################################
sub importBEDFile_Array{  
    #expects BED chrom = UCSC chrom field, e.g. 'chr13'
    #expects BED name = array ratio::bAlleleFreq::cnvValue
    my ($sample, $arrayType, $bedFile, $refSeq) = @_;
    print "importing $bedFile\n"; 
    my $name = arrayTableName($sample, $arrayType, $refSeq);
    my $statsTable = newTable('Stats', $name);  
    updateStat($statsTable, 'maxHeterozygous', 0.6);
    updateStat($statsTable, 'minHomozygous', 0.9);    
    my $arrayTable = newTable('Array', $name);     
    my $arrayFile = "$arrayTable.csv";    
    open my $bedFileH, "<", $bedFile or die "could not open $bedFile";      
    open my $arrayFileH, ">", $arrayFile;
    while (<$bedFileH>){
        chomp $_;
        my ($chr, $start, $end, $name) = split(" ", $_);
        my $chrom = $refSeqs{$refSeq}{$chr};           
        my ($ratio, $bFreq, $cnvValue) = split("::", $name);
        $ratio = 2 ** $ratio;
        my $zygosity = -1;
        $bFreq > -1 and $zygosity = abs($bFreq - 0.5) + 0.5;
        print $arrayFileH "$chrom,$start,$ratio,$ratio,$bFreq,$zygosity,$zygosity,$cnvValue\n";
    }
    close $bedFileH;
    close $arrayFileH; 
    loadData($arrayFile, $arrayTable, ',', "CHROMOSOME, POSITION, RATIO, NORMALIZEDRATIO, 
                                            BFREQUENCY, ZYGOSITY, NORMALIZEDZYGOSITY, CNVVALUE")
}
#-----------------------------------------------
sub normalizeSNPs{
    my ($testSample, $refSample, $arrayType, $refSeq) = @_;
    if($testSample eq $refSample or !$refSample){normalizeSNPsToSelf($testSample, $arrayType, $refSeq); return}    
    my $testName = arrayTableName($testSample, $arrayType, $refSeq);
    my $refName = arrayTableName($refSample, $arrayType, $refSeq);
    my $testArrayTable = getTableName('Array', $testName);    
    my $refArrayTable = getTableName('Array', $refName);    
    my $refStatsTable = getTableName('Stats', $refName);
    getStatistics($refStatsTable, \my%refStats);
    runSQL("SELECT t.CHROMOSOME, t.POSITION, t.RATIO, r.RATIO refRatio, 
                   t.BFREQUENCY, t.ZYGOSITY, r.ZYGOSITY refZygosity, t.CNVVALUE
             FROM $testArrayTable t, $refArrayTable r 
             WHERE t.CHROMOSOME || 'x' || t.POSITION = r.CHROMOSOME || 'x' || r.POSITION",
             \my($chrom, $pos, $ratio, $refRatio, $bFreq, $zygosity, $refZygosity, $cnvValue) );
    my $arrayFile = "$testArrayTable.csv";       
    open my $arrayFileH, ">", $arrayFile;
    while (fetchRow()){
        my $normZygosity = $zygosity / $refZygosity;        
        $refZygosity >= $refStats{minHomozygous} and $normZygosity = -1; 
        my $normRatio = $ratio;
        $refRatio != 0 and $normRatio = $ratio / $refRatio;
        print $arrayFileH "$chrom,$pos,$ratio,$normRatio,$bFreq,$zygosity,$normZygosity,$cnvValue\n";   
    } 
    close $arrayFileH; 
    $testArrayTable = newTable('Array', $testName);  
    loadData($arrayFile, $testArrayTable, ',', "CHROMOSOME, POSITION, RATIO, NORMALIZEDRATIO, 
                                                BFREQUENCY, ZYGOSITY, NORMALIZEDZYGOSITY, CNVVALUE")
}
sub normalizeSNPsToSelf{
    my ($sample, $arrayType, $refSeq) = @_;
    my $name = arrayTableName($sample, $arrayType, $refSeq);
    my $arrayTable = getTableName('Array', $name);      
    my $statsTable = getTableName('Stats', $name);
    getStatistics($statsTable, \my%stats);
    runSQL("UPDATE $arrayTable 
            SET NORMALIZEDZYGOSITY = -1 
            WHERE ZYGOSITY >= $stats{minHomozygous}");
    runSQL("UPDATE $arrayTable 
            SET NORMALIZEDZYGOSITY = ZYGOSITY / 0.5 
            WHERE ZYGOSITY< $stats{minHomozygous}"); 
}
################################################


################################################
#common
#-----------------------------------------------
sub stripBEDHeader{
    my ($bedFile) = @_;
    my $bedFileStripped = "$bedFile.stripped";    
    open my $bedFileH, "<", $bedFile or die "could not open $bedFile";     
    open my $bedFileSH, ">", $bedFileStripped;
    my $line1 = <$bedFileH>;
    while (<$bedFileH>){print $bedFileSH $_}
    close $bedFileH;
    close $bedFileSH;    
    return $bedFileStripped;
}
#-----------------------------------------------
sub liftOverBED{
    my ($bedFile, $liftFrom, $liftTo) = @_;
    print "lifting $bedFile from $liftFrom to $liftTo\n"; 
    my $lODir = "$param{vampPath}/liftOver";    
    my $bedFileOut = $bedFile.".$liftTo";
    system("$lODir/liftOver $bedFile $lODir/$liftFrom"."To\u$liftTo.over.chain $bedFileOut $bedFileOut.unmapped");
    return $bedFileOut;
}

################################################

1;

