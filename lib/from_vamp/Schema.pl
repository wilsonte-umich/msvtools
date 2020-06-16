#!/usr/bin/perl -w
use strict;
use warnings;

#########################################################################
#Schema.pl is a storehouse for script requires file naming conventions and db schema
#########################################################################

use vars(qw(%param %types %fields %fieldNames %partitions %archiveTables));

#########################################################################
#load vamp required scripts
#------------------------------------------------------------------------
-d "$param{vampPath}/bin" or die "\n$param{vampPath} is not a valid VAMP path\n";
foreach my $script(<$param{vampPath}/bin/*.pl>){ $script =~ m/ExecuteInstruction.pl/ or require $script }
requireFolders();  #require user scripts
#########################################################################


#########################################################################
#load user-supplied scripts
#------------------------------------------------------------------------
#load any scripts found in 'XXX/vamp/bin/usr' where XXX is either the
#home directory or vampPath, with the home directory taking precedence
#if a script of the same name is found in both locations
#Any subs in usr script must be accessed by command 'run subName parameterList'
#------------------------------------------------------------------------
my $path = "$ENV{HOME}/vamp/bin/usr";
if(-d $path) {
    foreach my $script(<$path/*.pl>){ require $script }
}
#########################################################################


#########################################################################
#subs that determine naming conventions
#------------------------------------------------------------------------
sub getDirectories{
    my ($sample, $dirsRef) = @_;
    $$dirsRef{sample} = "$param{inputPath}/$sample"; 
    $$dirsRef{read1} = "$$dirsRef{sample}/1";
    $$dirsRef{read2} = "$$dirsRef{sample}/2";
    $$dirsRef{read3} = "$$dirsRef{sample}/3";
    $$dirsRef{qseq1} = "$$dirsRef{read1}/qseq";
    $$dirsRef{qseq2} = "$$dirsRef{read2}/qseq";
    $$dirsRef{qseq3} = "$$dirsRef{read3}/qseq";
    $$dirsRef{prepared_read1} = "$$dirsRef{read1}/prepared";###########
    $$dirsRef{prepared_read2} = "$$dirsRef{read2}/prepared";########### 
    $$dirsRef{prepared_read3} = "$$dirsRef{read3}/prepared";###########
}
sub mkdirBasic {
    my ($dirsRef) = @_;
    -d $$dirsRef{sample} or mkdir $$dirsRef{sample};
    -d $$dirsRef{read1} or mkdir $$dirsRef{read1};
    -d $$dirsRef{read2} or mkdir $$dirsRef{read2};  
}
sub getReadFiles{
    my ($sample, $readType, $filesRef) = @_;
    my %dirs;
    getDirectories($sample, \%dirs);
    if ($readType eq 'solexa_fastq'){
        $$filesRef{read1} = "$dirs{read1}/$sample\_1_sequence.txt";
        $$filesRef{read2} = "$dirs{read2}/$sample\_2_sequence.txt"; 
        $$filesRef{read3} = "$dirs{read3}/$sample\_3_sequence.txt";
    } elsif ($readType eq 'fastq'){
        $$filesRef{read1} = "$dirs{read1}/$sample\_1.fq";
        $$filesRef{read2} = "$dirs{read2}/$sample\_2.fq";
        $$filesRef{read3} = "$dirs{read3}/$sample\_3.fq";
    } elsif ($readType eq 'fasta') {
        $$filesRef{read1} = "$dirs{read1}/$sample\_1.fa";
        $$filesRef{read2} = "$dirs{read2}/$sample\_2.fa";
        $$filesRef{read3} = "$dirs{read3}/$sample\_3.fa";
    } elsif ($readType eq 'purgedFasta') {
        $$filesRef{read1} = "$dirs{read1}/$sample\_1.fa.p";
        $$filesRef{read2} = "$dirs{read2}/$sample\_2.fa.p";
        $$filesRef{read3} = "$dirs{read3}/$sample\_3.fa.p";
    } elsif ($readType eq 'sam') {
        $$filesRef{bam} = "$dirs{sample}/$sample.bam";    
        $$filesRef{sam} = "$dirs{sample}/$sample.sam";           
    } else {
        die "'$readType' is not a recognized readType"
    }
}
sub getMapFiles{
    my ($sample, $filesRef) = @_;
    getDirectories($sample, \my%dirs);
    my ($sampleSuffix, $readSuffix);
    if ($param{mapType} eq 'pass') {
        $readSuffix = 'gff';
    } elsif ($param{mapType} eq 'bowtie') {
        $readSuffix = 'bowtie';
    } elsif($param{mapType} eq 'sam') {
        $readSuffix = 'sam'; 
    } elsif($param{mapType} eq 'pairedSam') {
        $sampleSuffix = 'sam'; 
    } elsif($param{mapType} eq 'bam') {
        $readSuffix = 'bam'; 
    } elsif($param{mapType} eq 'pairedBam') {
        $sampleSuffix = 'bam'; 
    } elsif($param{mapType} eq 'bwa') {
        $readSuffix = 'sai'; 
        $sampleSuffix = 'bam';
    } else {
        die "$param{mapType} is not a recognized mapType";
    }
    $sampleSuffix and $$filesRef{sample} = "$dirs{sample}/$sample.$sampleSuffix";
    $readSuffix and $$filesRef{read1} = "$dirs{read1}/$sample\_1.$readSuffix";
    $readSuffix and $$filesRef{read2} = "$dirs{read2}/$sample\_2.$readSuffix";
    $readSuffix and $$filesRef{read3} = "$dirs{read3}/$sample\_3.$readSuffix";
}
sub getTableName{
    my ($tableType, $sample) = @_;
    return "\U$tableType\_$sample"; #force lane name to upper case for Oracle
}
sub splitSampleName {
    my ($sample, $dontDefaultOut) = @_;
    my ($sampleIn, $sampleOut);
    if ($sample =~ m/::/){
        ($sampleIn, $sampleOut) = split('::', $sample);
    } elsif ($sample =~ m/=/){
        ($sampleIn, $sampleOut) = split('=', $sample);
    } else {
        ($sampleIn, $sampleOut) = ($sample, $sample);
        $dontDefaultOut and $sampleOut = undef;
    }
    return ($sampleIn, $sampleOut);
}
#########################################################################


#########################################################################
#type enumerations
#------------------------------------------------------------------------
%types = (Pairs => {Convergent=>0, Divergent=>1, Colinear=>2, DiffChrom=>3, SingleRead=>4},
          Discs => {None=>0, A=>1, C=>2, G=>3, T=>4, '.' => 5, '?'=>5, N=>5, X=>5, Missing=>6, Extra=>7},
          Frags => {Discard => -1, Normal=>0, Deletion=>1, Insertion=>2, Inversion=>4, Duplication=>8,
                    DiffChrom=>16, ReverseNormal=>32, TooSmall=>64, Rooted=>128, SingleRead=>256, Anchored=>512}, 
          Consequences => {Silent=>0, Missense=>1, Nonsense=>2, InDel=>3},
          RevDiscs => {1=>'A', 2=>'C', 3=>'G', 4=>'T', 5=>'N'}, 
          RevDiscsFull => {0=>'-', 1=>'A', 2=>'C', 3=>'G', 4=>'T', 5=>'N', 6=>'Missing', 7=>'Extra'},
          Genes => {All=>0, Unq=>1, Ovr=>2, Nvr=>3} );
$types{Sets} = $types{Frags};
$types{Sets}{Z} = 3;
#########################################################################


#########################################################################
#table field definitions
#------------------------------------------------------------------------
$fields{Reads} = "PAIRID NUMBER, MAPS VARCHAR2(4000)";
$fields{Pairs} = "PAIRTYPE NUMBER,
                    CHROMOSOME1 NUMBER, POSITION1 NUMBER, LENGTH1 NUMBER, STRAND1 NUMBER, DISCREPANCIES1 NUMBER, NHITS1 NUMBER,
                    CHROMOSOME2 NUMBER, POSITION2 NUMBER, LENGTH2 NUMBER, STRAND2 NUMBER, DISCREPANCIES2 NUMBER, NHITS2 NUMBER,
                    FRAGMENTSIZE NUMBER, PAIRID NUMBER, NFRAGS NUMBER";
$fields{Discs} = "DISCREPANCYID NUMBER, 
                    CHROMOSOME NUMBER, POSITION NUMBER, DISCREPANCYTYPE NUMBER, EXTRA NUMBER,
                    S1_SAMECOUNT NUMBER, S1_CONSISTENTCOUNT NUMBER, S1_READCOUNT NUMBER, 
                    S2_SAMECOUNT NUMBER, S2_CONSISTENTCOUNT NUMBER, S2_READCOUNT NUMBER,
                    SAMELOD NUMBER(*,3), CONSISTENTLOD NUMBER(*,3), CONSEQUENCE NUMBER"; 
$fields{Stats} = "STATNAME VARCHAR2(255), STATVALUE NUMBER(*,3)";
$fields{Frags} = "FRAGMENTID NUMBER, FRAGMENTTYPE NUMBER,
                    CHROMOSOME1 NUMBER, POSITION1 NUMBER, LENGTH1 NUMBER, STRAND1 NUMBER, DISCREPANCIES1 NUMBER, NHITS1 NUMBER,
                    CHROMOSOME2 NUMBER, POSITION2 NUMBER, LENGTH2 NUMBER, STRAND2 NUMBER, DISCREPANCIES2 NUMBER, NHITS2 NUMBER,
                    FRAGMENTSIZE NUMBER, PAIRID NUMBER, NFRAGS NUMBER,
                    EVENTSIZE NUMBER, STDEVNORMAL NUMBER, ENDTOLERANCE NUMBER, 
                    NSETSFRAG NUMBER, NSETSPAIR NUMBER";    
$fields{FMap} = "CHROMOSOME NUMBER, POSITION NUMBER, COVERAGE NUMBER, NORMALIZEDCOVERAGE NUMBER(*,3)";
$fields{TMap} = "CHROMOSOME NUMBER, POSITION NUMBER, STRAND NUMBER, NONUVCOVERAGE NUMBER, UVCOVERAGE NUMBER";
$fields{RMap} = $fields{FMap};  
$fields{HMap} = "CHROMOSOME NUMBER, POSITION NUMBER, COVERAGE NUMBER, DENSITY NUMBER(*,5), 
                 NORMALIZEDCOVERAGE NUMBER(*,5), NORMALIZEDDENSITY NUMBER(*,5)";
$fields{GMap} = "NAME2 VARCHAR2(255), COVERAGE NUMBER, DENSITY NUMBER(*,5), 
                 NORMALIZEDCOVERAGE NUMBER(*,5), NORMALIZEDDENSITY NUMBER(*,5)";
$fields{IMap} = "NAME2 VARCHAR2(255), CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER,
                 COVERAGE NUMBER, DENSITY NUMBER(*,5),
                 NORMALIZEDCOVERAGE NUMBER(*,5), NORMALIZEDDENSITY NUMBER(*,5)";
$fields{EMap} = $fields{IMap};       
$fields{ZMap} = "CHROMOSOME NUMBER, POSITION NUMBER, COVERAGE NUMBER, MEAN NUMBER, Z NUMBER(*,1)";                         
#$fields{FRatio} = "CHROMOSOME NUMBER, POSITION NUMBER,
#                    COVERAGERATIO NUMBER(*,3), COVERAGERATIO1E4 NUMBER(*,3), COVERAGERATIO1E5 NUMBER(*,3)";                            
$fields{Sets} = "SETID NUMBER, SETTYPE NUMBER,
                    CHROMOSOME1 NUMBER, CHROMOSOME2 NUMBER,
                    SPANSTART NUMBER, SPANEND NUMBER, OVERLAPSTART NUMBER, OVERLAPEND NUMBER,
                    FRAGMENTMEAN NUMBER, EVENTMEAN NUMBER, STDEV NUMBER, 
                    MEANDELTA NUMBER(*,3), FRACTIONDELTA2 NUMBER(*,3), FRACTIONOVERLAPLEFT NUMBER(*,3), FRACTIONOVERLAPRIGHT NUMBER(*,3),
                    STRAND1 NUMBER, STRAND2 NUMBER,
                    NFRAGSTOTAL NUMBER, NFRAGSSAMPLE NUMBER, FRACBADFRAGS NUMBER(*,3),
                    MARK NUMBER, DESCRIPTION VARCHAR2(255)";                
$fields{IDs} = "SETID NUMBER, ID_ NUMBER";
$fields{CDS} = "CDSID NUMBER, NAME1 VARCHAR2(255), NAME2 VARCHAR2(255),
                CHROMOSOME NUMBER, STRAND NUMBER, EXON NUMBER,
                START_ NUMBER, CARRYOVER NUMBER, END_ NUMBER";
$fields{NCDS} = "NCDSID NUMBER, NAME1 VARCHAR2(255), NAME2 VARCHAR2(255),
                CHROMOSOME NUMBER, STRAND NUMBER, EXON NUMBER, 
                START_ NUMBER, END_ NUMBER";
$fields{h} = "SERIES NUMBER, X NUMBER, Y NUMBER"; #h = histogram, kept short since name appended to parent table name
$fields{rmsk} = "CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER, STRAND NUMBER, 
                 REPEATNAME VARCHAR(255), REPEATCLASS VARCHAR(255), REPEATFAMILY VARCHAR(255)";    
$fields{CNVs} = "CNVNAME VARCHAR2(255), CNVTYPE NUMBER, CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER";   
$fields{SNPs} = "STRAIN VARCHAR2(255), CHROMOSOME NUMBER, POSITION NUMBER, DISCREPANCYTYPE NUMBER";   
$fields{Array} = "CHROMOSOME NUMBER, POSITION NUMBER, RATIO NUMBER(*,5), NORMALIZEDRATIO NUMBER(*,5),
                  BFREQUENCY NUMBER(*,5), ZYGOSITY NUMBER(*,2), NORMALIZEDZYGOSITY NUMBER(*,2), CNVVALUE NUMBER(*,5),
                  GCSCORE NUMBER(*,3), X NUMBER(*,3), Y NUMBER(*,3), THETA NUMBER(*,3), R NUMBER(*,3), BF NUMBER(*,3)";    
$fieldNames{Arrays} = "CHROMOSOME VARCHAR(255), POSITION NUMBER, SAMPLE VARCHAR(255), GROUP VARCHAR(255), 
                       RAW NUMBER(*,5), SEXED NUMBER(*,5), CENTROID NUMBER(*,5), STATE NUMBER(*,5), SEGMENT NUMBER(*,5)";                    
$fields{NMA} = "NMAID NUMBER, NMATYPE VARCHAR(255), NPROBES NUMBER, NSTDDEV NUMBER(*,1), 
                CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER, CNVSIZE NUMBER, NPROBESINSET NUMBER, 
                TESTMEAN NUMBER(*,3), TESTSTDEV NUMBER(*,3), REFERENCEMEAN NUMBER(*,3), REFERENCESTDEV NUMBER(*,3), 
                NORMALIZEDRATIO NUMBER(*,3), NORMALIZEDZYGOSITY NUMBER(*,3), ZSCORE NUMBER(*,1),
                COPYNUMBER NUMBER, DESCRIPTION VARCHAR2(255) ";  
$fields{LD} = "CHROMOSOME NUMBER, REFERENCEPOSITION NUMBER, IMPUTEDPOSITION NUMBER, 
                DPRIME NUMBER(*,3), RSQUARED NUMBER(*,3), LOD NUMBER(*,3)"; 
$fields{AF} = "CHROMOSOME NUMBER, POSITION NUMBER, UNINFORMATIVEFREQ NUMBER(*,3)"; 
$fields{Events} = "EVENTID NUMBER, EVENTTYPE NUMBER,
                    CHROMOSOME NUMBER, SPANSTART NUMBER, SPANEND NUMBER, 
                    NSETS NUMBER ";  
$fields{CNCs} = "CNCID NUMBER, EVENTID NUMBER, CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER, 
                 COPYNUMBER NUMBER, NORMALIZEDCOVERAGE NUMBER(*,3) ";  
$fields{Marks} = "SETID NUMBER, USERINITIALS VARCHAR2(5), MARK NUMBER, DESCRIPTION VARCHAR2(255)"; 
$fields{Hits} = "CHROMOSOME NUMBER, POSITION NUMBER, STRAND NUMBER, COUNT_ NUMBER";
$fields{HB} = "HBID NUMBER, THRESHOLD NUMBER, BINSIZE NUMBER, CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER, 
               NBINS NUMBER, NORMALIZEDDENSITY NUMBER(*,5), INGENE NUMBER";
$fields{Unq} = "CHROMOSOME NUMBER, POSITION NUMBER, BIN25MB NUMBER";
$fields{UnqB} = "CHROMOSOME NUMBER, START_ NUMBER, END_ NUMBER, UNIQUE NUMBER";
#------------------------------------------------------------------------
$fieldNames{Pairs} = "PAIRTYPE,
                      CHROMOSOME1, POSITION1, LENGTH1, STRAND1, DISCREPANCIES1, NHITS1,
                      CHROMOSOME2, POSITION2, LENGTH2, STRAND2, DISCREPANCIES2, NHITS2,
                      FRAGMENTSIZE, PAIRID, NFRAGS";
$fieldNames{Frags} = "FRAGMENTID, FRAGMENTTYPE,
                      CHROMOSOME1, POSITION1, LENGTH1, STRAND1, DISCREPANCIES1, NHITS1,
                      CHROMOSOME2, POSITION2, LENGTH2, STRAND2, DISCREPANCIES2, NHITS2,
                      FRAGMENTSIZE, PAIRID, NFRAGS,
                      EVENTSIZE, STDEVNORMAL, ENDTOLERANCE, 
                      NSETSFRAG, NSETSPAIR";                
$fieldNames{Sets} = "SETID, SETTYPE,
                    CHROMOSOME1, CHROMOSOME2,
                    SPANSTART, SPANEND, OVERLAPSTART, OVERLAPEND,
                    FRAGMENTMEAN, EVENTMEAN, STDEV, 
                    MEANDELTA, FRACTIONDELTA2, FRACTIONOVERLAPLEFT, FRACTIONOVERLAPRIGHT,
                    STRAND1, STRAND2,
                    NFRAGSTOTAL, NFRAGSSAMPLE, FRACBADFRAGS,
                    MARK, DESCRIPTION"; 
$fieldNames{CNVs} = "CNVNAME, CNVTYPE, CHROMOSOME, START_, END_"; 
$fieldNames{SNPs} = "STRAIN, CHROMOSOME, POSITION, DISCREPANCYTYPE";   
$fieldNames{Array} = "CHROMOSOME, POSITION, RATIO, NORMALIZEDRATIO, BFREQUENCY, ZYGOSITY, NORMALIZEDZYGOSITY, 
                      CNVVALUE, GCSCORE, X, Y, THETA, R, BF";       
$fieldNames{Arrays} = "CHROMOSOME, POSITION, SAMPLE, GROUP, RAW, SEXED, CENTROID, STATE, SEGMENT";                 
$fieldNames{NMA} = "NMAID, NMATYPE, NPROBES, NSTDDEV, 
                    CHROMOSOME, START_, END_, CNVSIZE, NPROBESINSET, 
                    TESTMEAN, TESTSTDEV, REFERENCEMEAN, REFERENCESTDEV, 
                    NORMALIZEDRATIO, NORMALIZEDZYGOSITY, ZSCORE,
                    COPYNUMBER, DESCRIPTION";
$fieldNames{LD} = "CHROMOSOME, REFERENCEPOSITION, IMPUTEDPOSITION, DPRIME, RSQUARED, LOD"; 
$fieldNames{AF} = "CHROMOSOME, POSITION, UNINFORMATIVEFREQ";    
$fieldNames{Events} = "EVENTID, EVENTTYPE, CHROMOSOME, SPANSTART, SPANEND, NSETS";  
$fieldNames{CNCs} = "CNCID, EVENTID, CHROMOSOME, START_, END_, COPYNUMBER, NORMALIZEDCOVERAGE"; 
$fieldNames{Marks} = "SETID, USERINITIALS, MARK, DESCRIPTION"; 
$fieldNames{Hits} = "CHROMOSOME, POSITION, STRAND, COUNT_"; 
$fieldNames{FMap} = "CHROMOSOME, POSITION, COVERAGE, NORMALIZEDCOVERAGE";
$fieldNames{TMap} = "CHROMOSOME, POSITION, STRAND, NONUVCOVERAGE, UVCOVERAGE";
$fieldNames{RMap} = $fieldNames{FMap};  
$fieldNames{HMap} = "CHROMOSOME, POSITION, COVERAGE, DENSITY, NORMALIZEDCOVERAGE, NORMALIZEDDENSITY";  
$fieldNames{GMap} = "NAME2, COVERAGE, DENSITY, NORMALIZEDCOVERAGE, NORMALIZEDDENSITY";
$fieldNames{IMap} = "NAME2, CHROMOSOME, START_, END_, COVERAGE, DENSITY, NORMALIZEDCOVERAGE, NORMALIZEDDENSITY";
$fieldNames{EMap} = $fieldNames{IMap};
$fieldNames{HB} = "HBID, THRESHOLD, BINSIZE, CHROMOSOME, START_, END_, NBINS, NORMALIZEDDENSITY, INGENE";
$fieldNames{Wig} = "CHROMOSOME, POSITION, VALUE";
$fieldNames{Unq} = "CHROMOSOME, POSITION, BIN25MB";
$fieldNames{UnqB} = "CHROMOSOME, START_, END_, UNIQUE";
#########################################################################
  
               
#########################################################################
#table partition definitions
# = [$chromPartField, $typeSubPartField, $typeSubPartsRef, $partByPosField]
#------------------------------------------------------------------------     
$partitions{Reads} = []; #not partitioned
$partitions{LU} = [];
$partitions{Pairs} = ["CHROMOSOME1", "PAIRTYPE", $types{Pairs}];
#must purge ? and . out of Disc Types for subpartition, since will cause Oracle name error
my %discSubParts = (None=>0, A=>1, C=>2, G=>3, T=>4, N=>5, Missing=>6, Extra=>7);
$partitions{Discs} = ["CHROMOSOME", "DISCREPANCYTYPE", \%discSubParts];
$partitions{Stats} = [];
$partitions{Frags} = ["CHROMOSOME1", "FRAGMENTTYPE", $types{Frags}];
$partitions{FMap} = ["CHROMOSOME"];
$partitions{TMap} = $partitions{FMap};
$partitions{RMap} = $partitions{FMap};
$partitions{HMap} = $partitions{FMap};
$partitions{GMap} = [];
$partitions{IMap} = ["CHROMOSOME"];
$partitions{EMap} = $partitions{IMap};
$partitions{ZMap} = $partitions{FMap};
#$partitions{FRatio} = ["CHROMOSOME"];
$partitions{Sets} = ["CHROMOSOME1", "SETTYPE", $types{Sets}];
$partitions{IDs} = [];
$partitions{CDS} = [];
$partitions{NCDS} = [];
$partitions{h} = [];
$partitions{rmsk} = ["CHROMOSOME"];
$partitions{CNVs} = ["CHROMOSOME"];
$partitions{SNPs} = ["CHROMOSOME"];
$partitions{Array} = ["CHROMOSOME"];
$partitions{Arrays} = ["CHROMOSOME"];
$partitions{NMA} = ["CHROMOSOME"];
$partitions{LD} = ["CHROMOSOME"];
$partitions{AF} = ["CHROMOSOME"];
$partitions{Events} = [];
$partitions{CNCs} = [];
$partitions{Marks} = [];

$partitions{Hits} = ["CHROMOSOME"]; #, undef, undef, "BIN25MB"

$partitions{HB} = ["CHROMOSOME"];
$partitions{Wig} = ["CHROMOSOME"];

$partitions{Unq} = ["CHROMOSOME", undef, undef, "BIN25MB"];
$partitions{UnqB} = ["CHROMOSOME", ];
#########################################################################

#########################################################################
#archiving information
#------------------------------------------------------------------------ 
$archiveTables{sequence} = ['Reads', 'Pairs', 'Discs', 'Stats', 'Frags', 'CDS', 'NCDS',
                             'FMap', 'RMap', 'ZMap', 'Sets', 'Events', 'CNCs', 'Marks','Hits'];
#########################################################################

#########################################################################
#genetic code
#------------------------------------------------------------------------ 
our %geneticCode = (
        # * - Stop
        'TAA'=>'x','TAG'=>'x','TGA'=>'x',
        # A - Alanine
        'GCT'=>'A','GCC'=>'A','GCA'=>'A','GCG'=>'A',
        # C - Cysteine
        'TGT'=>'C','TGC'=>'C',
        # D - Aspartic Acid
        'GAT'=>'D','GAC'=>'D',
        # E - Glutamic Acid
        'GAA'=>'E','GAG'=>'E',
        # F - Phenylalanine
        'TTT'=>'F','TTC'=>'F',
        # G - Glycine
        'GGT'=>'G','GGC'=>'G','GGA'=>'G','GGG'=>'G',
        # H - Histidine
        'CAT'=>'H','CAC'=>'H',
        # I - Isoleucine
        'ATT'=>'I','ATC'=>'I','ATA'=>'I',
        # K - Lysine
        'AAA'=>'K','AAG'=>'K',
        # L - Leucine
        'CTT'=>'L','CTC'=>'L','CTA'=>'L','CTG'=>'L',
        'TTA'=>'L','TTG'=>'L',
        # M - Methionine
        'ATG'=>'M',
        # N - Asparagine
        'AAT'=>'N','AAC'=>'N',
        # P - Proline
        'CCT'=>'P','CCC'=>'P','CCA'=>'P','CCG'=>'P',
        # Q - Glutamine
        'CAA'=>'Q','CAG'=>'Q',
        # R - Arginine
        'CGT'=>'R','CGC'=>'R','CGA'=>'R','CGG'=>'R',
        'AGA'=>'R','AGG'=>'R',
        # S - Serine
        'TCT'=>'S','TCC'=>'S','TCA'=>'S','TCG'=>'S',
        'AGT'=>'S','AGC'=>'S',
        # T - Threonine
        'ACT'=>'T','ACC'=>'T','ACA'=>'T','ACG'=>'T',
        # V - Valine
        'GTT'=>'V','GTC'=>'V','GTA'=>'V','GTG'=>'V',
        # W - Tryptophan
        'TGG'=>'W',
        # Y - Tyrosine
        'TAT'=>'Y','TAC'=>'Y',
);
#########################################################################


1;

