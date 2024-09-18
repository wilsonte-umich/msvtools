use strict;
use warnings;

#========================================================================
# 'formats.pl' has file column formats
#========================================================================

my $i = 0;
our %prbCol = map { $_ => $i++ } (
    'CHROM','START','POS','PRB_NAME','IS_NA','STRAND',
    'LRR','BAF',                # refine=median, segment=rescaled and adjusted values used in HMM
    'AS_IN','CN_IN','INF_IN',   # input state model
    'LRR_RAW','BAF_RAW',        # original values prior to refine or segment actions
    # 'ot','os',
    'hs',                       # HMM indices
    'AS_OUT','CN_OUT','INF_OUT',# output state model
    'CNC','LOH',                # SV call metadata
    # FRAC_GC                   # in genome model files
);

$i = 0;
our %mdlCol = map { $_ => $i++ } (
    'CHROM','START','END','AS_OUT','CN_OUT','STRAND',
    'N_PROBES','LRR_MEAN','LRR_SD','ZYG_MEAN','ZYG_SD', # aggregrate probe raw data over probes in segment
    # relative likelihood of each allelic state
);

$i = 0;
our %svCol = map { $_ => $i++ } (
    'CHROM','START','END','SAMPLE','SV_ID','STRAND',
    'LRR_MEAN','LRR_SD','ZYG_MEAN','ZYG_SD',
    'AS_IN','CN_IN',
    'AS_OUT','CN_OUT', 
    'CNC', 'LOH',
    'N_PROBES', 'N_INF', 'SPAN',
    'MDL_LL','SV_LL','REL_LL',
    'CN_CLC', 'CN_P',
    'FLANK_START', 'FLANK_END'
);

our @RLL_COL = ('REGION_ID',
                'QRLL',   'MIN_TRLL',   'MAX_TRLL');
our @RLL_SMP_COL = ('MIN_TRLL_SMP',  'MAX_TRLL_SMP', 'OVLP_SMP');
our %cmpCol = (%svCol, map { $_ => $i++ } (@RLL_COL, @RLL_SMP_COL));
   
1;

#our @RLL_COL = ('REGION_ID',
#                'QRLL',   'MIN_TRLL',   'MAX_TRLL',   'MIN_BIAS',   'MAX_BIAS',   
#                'QERLL',  'MIN_TERLL',  'MAX_TERLL',  'MIN_EBIAS',  'MAX_EBIAS');
#our @RLL_SMP_COL = ('MIN_TRLL_SMP',  'MAX_TRLL_SMP',  'MIN_BIAS_SMP',  'MAX_BIAS_SMP',
#                    'MIN_TERLL_SMP', 'MAX_TERLL_SMP', 'MIN_EBIAS_SMP', 'MAX_EBIAS_SMP');
#our %cmpCol = (%svCol, map { $_ => $i++ } (@RLL_COL, @RLL_SMP_COL));

#our @RLL_COL = ('REGION_ID',
#                'QRLL',   'MIN_TRLL',   'MAX_TRLL',   'MIN_RATIO',   'MAX_RATIO',   
#                'QERLL',  'MIN_TERLL',  'MAX_TERLL',  'MIN_ERATIO',  'MAX_ERATIO');
#our @RLL_SMP_COL = ('MIN_TRLL_SMP',  'MAX_TRLL_SMP',
#                    'MIN_TERLL_SMP', 'MAX_TERLL_SMP');
#our %cmpCol = (%svCol, map { $_ => $i++ } (@RLL_COL, @RLL_SMP_COL));