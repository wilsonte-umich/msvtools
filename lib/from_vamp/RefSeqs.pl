#!/usr/bin/perl -w
use strict;
use warnings;

#######################################################
#refSeqs.pl contains chromosome definitions for refSeqs,
#noting that VAMP converts _everything_ to numbers only.
#Definitions are contained in a hash entry as
# {nChrom => number,
# chrName => number, etc.}
#A hash entry MUST exist in this file for any specified refSeqBase (or refSeq if refSeqBase not named)!
#VAMP expects to find sequential entries for 1 to nChrom, but there may be more
#see hg19 example below

#TODO: this is cumbersome, should really be replaced by
#an ability to read refSeq chromosomes from an external file

#######################################################

use vars(qw(%param %types %fields %refSeqs));

sub nChrom{return $refSeqs{$param{refSeqBase}}{nChrom}}

#when the hashes below are called, the primary key always corresponds to parameter refSeqBase
#noting that refSeqBase defaults to refSeq if refSeqBase not explicity defined

$refSeqs{hg19} = {
    nChrom => 24,
    #standard chromosomes
    chr1 => 1,
    chr2 => 2,
    chr3 => 3,
    chr4 => 4,
    chr5 => 5,
    chr6 => 6,
    chr7 => 7,
    chr8 => 8,
    chr9 => 9,
    chr10 => 10,
    chr11 => 11,
    chr12 => 12,
    chr13 => 13,
    chr14 => 14,
    chr15 => 15,
    chr16 => 16,
    chr17 => 17,
    chr18 => 18,
    chr19 => 19,
    chr20 => 20,
    chr21 => 21,
    chr22 => 22,
    chrX => 23,
    chrY => 24,
    #mitochondrial genome
    chrM => 25,
    #known alternative haplotypes
    chr4_ctg9_hap1 => 1041,
    chr6_apd_hap1 => 1061,
    chr6_cox_hap2 => 1062,
    chr6_dbb_hap3 => 1063,
    chr6_mann_hap4 => 1064,
    chr6_mcf_hap5 => 1065,
    chr6_qbl_hap6 => 1066,
    chr6_ssto_hap7 => 1067,
    chr17_ctg5_hap1 => 1171,
    #unplaced sequences, known chromosomes
    chr1_gl000191_random => 2011,
    chr1_gl000192_random => 2012,
    chr4_gl000193_random => 2041,
    chr4_gl000194_random => 2042,
    chr7_gl000195_random => 2071,
    chr8_gl000196_random => 2081,
    chr8_gl000197_random => 2082,
    chr9_gl000198_random => 2091,
    chr9_gl000199_random => 2092,
    chr9_gl000200_random => 2093,
    chr9_gl000201_random => 2094,
    chr11_gl000202_random => 2111,
    chr17_gl000203_random => 2171,
    chr17_gl000204_random => 2172,
    chr17_gl000205_random => 2173,
    chr17_gl000206_random => 2174,
    chr18_gl000207_random => 2181,
    chr19_gl000208_random => 2191,
    chr19_gl000209_random => 2192,
    chr21_gl000210_random => 2211,
    #unplaced sequences, unknown chromosomes
    chrUn_gl000211 => 301,
    chrUn_gl000212 => 302,
    chrUn_gl000213 => 303,
    chrUn_gl000214 => 304,
    chrUn_gl000215 => 305,
    chrUn_gl000216 => 306,
    chrUn_gl000217 => 307,
    chrUn_gl000218 => 308,
    chrUn_gl000219 => 309,
    chrUn_gl000220 => 310,
    chrUn_gl000221 => 311,
    chrUn_gl000222 => 312,
    chrUn_gl000223 => 313,
    chrUn_gl000224 => 314,
    chrUn_gl000225 => 315,
    chrUn_gl000226 => 316,
    chrUn_gl000227 => 317,
    chrUn_gl000228 => 318,
    chrUn_gl000229 => 319,
    chrUn_gl000230 => 320,
    chrUn_gl000231 => 321,
    chrUn_gl000232 => 322,
    chrUn_gl000233 => 323,
    chrUn_gl000234 => 324,
    chrUn_gl000235 => 325,
    chrUn_gl000236 => 326,
    chrUn_gl000237 => 327,
    chrUn_gl000238 => 328,
    chrUn_gl000239 => 329,
    chrUn_gl000240 => 330,
    chrUn_gl000241 => 331,
    chrUn_gl000242 => 332,
    chrUn_gl000243 => 333,
    chrUn_gl000244 => 334,
    chrUn_gl000245 => 335,
    chrUn_gl000246 => 336,
    chrUn_gl000247 => 337,
    chrUn_gl000248 => 338,
    chrUn_gl000249 => 339
};
$refSeqs{hg18} = $refSeqs{hg19};

$refSeqs{SC10012003} = {
    nChrom => 16,
    #standard chromosomes
    chr1 => 1,
    chr2 => 2,
    chr3 => 3,
    chr4 => 4,
    chr5 => 5,
    chr6 => 6,
    chr7 => 7,
    chr8 => 8,
    chr9 => 9,
    chr10 => 10,
    chr11 => 11,
    chr12 => 12,
    chr13 => 13,
    chr14 => 14,
    chr15 => 15,
    chr16 => 16,
    #mitochondrial genome
    chrM => 25
};

$refSeqs{sacCer2} = {
    nChrom => 16,
    #standard chromosomes
    chrI => 1,
    chrII => 2,
    chrIII => 3,
    chrIV => 4,
    chrV => 5,
    chrVI => 6,
    chrVII => 7,
    chrVIII => 8,
    chrIX => 9,
    chrX => 10,
    chrXI => 11,
    chrXII => 12,
    chrXIII => 13,
    chrXIV => 14,
    chrXV => 15,
    chrXVI => 16,
    CHR1 => 1,
    CHR2 => 2,
    CHR3 => 3,
    CHR4 => 4,
    CHR5 => 5,
    CHR6 => 6,
    CHR7 => 7,
    CHR8 => 8,
    CHR9 => 9,
    CHR10 => 10,
    CHR11 => 11,
    CHR12 => 12,
    CHR13 => 13,
    CHR14 => 14,
    CHR15 => 15,
    CHR16 => 16,
    CHRM => 25,
    chr1 => 1,
    chr2 => 2,
    chr3 => 3,
    chr4 => 4,
    chr5 => 5,
    chr6 => 6,
    chr7 => 7,
    chr8 => 8,
    chr9 => 9,
    chr10 => 10,
    chr11 => 11,
    chr12 => 12,
    chr13 => 13,
    chr14 => 14,
    chr15 => 15,
    chr16 => 16,
    #mitochondrial genome
    chrM => 25,
    #2 micron plasmid
    '2micron' => 26
};
$refSeqs{sacCerSNP} = $refSeqs{sacCer2};

our %reverseRefSeqs;
$reverseRefSeqs{sacCer2} = {
    1 => 'chrI',
    2 => 'chrII',
    3 => 'chrIII',
    4 => 'chrIV',
    5 => 'chrV',
    6 => 'chrVI',
    7 => 'chrVII',
    8 => 'chrVIII',
    9 => 'chrIX',
    10 => 'chrX',
    11 => 'chrXI',
    12 => 'chrXII',
    13 => 'chrXIII',
    14 => 'chrXIV',
    15 => 'chrXV',
    16 => 'chrXVI',
    25 => 'chrM',
    26 => '2micron'
};
$reverseRefSeqs{sacCerSNP} = $reverseRefSeqs{sacCer2};

$reverseRefSeqs{hg19} = {
    1 => 'chr1',
    2 => 'chr2',
    3 => 'chr3',
    4 => 'chr4',
    5 => 'chr5',
    6 => 'chr6',
    7 => 'chr7',
    8 => 'chr8',
    9 => 'chr9',
    10 => 'chr10',
    11 => 'chr11',
    12 => 'chr12',
    13 => 'chr13',
    14 => 'chr14',
    15 => 'chr15',
    16 => 'chr16',
    17 => 'chr17',
    18 => 'chr18',
    19 => 'chr19',
    20 => 'chr20',
    21 => 'chr21',
    22 => 'chr22',
    23 => 'chrX',
    24 => 'chrY'
};
$reverseRefSeqs{hg18} = $reverseRefSeqs{hg19};

$refSeqs{mm9} = {
    nChrom => 21,
    #standard chromosomes
    chr1 => 1,
    chr2 => 2,
    chr3 => 3,
    chr4 => 4,
    chr5 => 5,
    chr6 => 6,
    chr7 => 7,
    chr8 => 8,
    chr9 => 9,
    chr10 => 10,
    chr11 => 11,
    chr12 => 12,
    chr13 => 13,
    chr14 => 14,
    chr15 => 15,
    chr16 => 16,
    chr17 => 17,
    chr18 => 18,
    chr19 => 19,
    chrX => 20,
    chrY => 21,
    #mitochondrial genome
    chrM => 22,
    #unplaced sequences
    chr1_random => 2011,    
    chr3_random => 2031,
    chr4_random => 2041,
    chr5_random => 2051,
    chr7_random => 2071,
    chr8_random => 2081,
    chr9_random => 2091,   
    chr13_random => 2131,
    chr16_random => 2161,
    chr17_random => 2171,
    chrX_random => 2201,
    chrY_random => 2211,
    chrUn_random => 2001,
};

$reverseRefSeqs{mm9} = {
    1 => 'chr1',
    2 => 'chr2',
    3 => 'chr3',
    4 => 'chr4',
    5 => 'chr5',
    6 => 'chr6',
    7 => 'chr7',
    8 => 'chr8',
    9 => 'chr9',
    10 => 'chr10',
    11 => 'chr11',
    12 => 'chr12',
    13 => 'chr13',
    14 => 'chr14',
    15 => 'chr15',
    16 => 'chr16',
    17 => 'chr17',
    18 => 'chr18',
    19 => 'chr19',
    20 => 'chrX',
    21 => 'chrY'
};

$refSeqs{Ty1} = {
    nChrom => 1,
    YHRCTy1_1 => 1
};
$reverseRefSeqs{Ty1} = {
    1=>'YHRCTy1_1'
};
$refSeqs{chr2Ty1} = {
    nChrom => 1,
    chr2Ty1 => 1
};
$reverseRefSeqs{chr2Ty1} = {
    1=>'chr2Ty1'
};

$refSeqs{YP} = {
    nChrom => 2,
    PPP1CB => 1,
    YPEL5 => 2
};
$reverseRefSeqs{YP} = {
    1 => 'PPP1CB',
    2 => 'YPEL5'
};

$refSeqs{YP2} = $refSeqs{YP};
$reverseRefSeqs{YP2} =$reverseRefSeqs{YP};

$refSeqs{LR} = {
    nChrom => 2,
    left => 1,
    right => 2
};
$reverseRefSeqs{LR} = {
    1 => 'left',
    2 => 'right'
};

$refSeqs{phiX} = {
    nChrom => 1,
    phiX => 1
};
$reverseRefSeqs{phiX} = {
    1 => 'phiX'
};

1;

