#!/usr/bin/perl
use strict;
use warnings;

use vars(qw(%param %types %fields %refSeqs %reverseRefSeqs));

my @cnvs = (
['A3A2',7,61379057,61440507],
['A3A2',11,6476303,6564800],
['A3A2',15,69469199,69499819],
['A3A2',18,7630270,7864089],
['A3A2',3,117988930,118056525],
['A3A2',3,118056525,118673931],
['A3A2',3,118673931,118751078],
['A3A2',1,172809126,172859441],
['A3A2',9,130306938,130696029],
['A3A2',10,111820163,114253039],
['A3A2',4,123973232,123999958],
['A3A2',4,124044944,124068699],
['A1A1',7,100415834,100435634],
['A1A1',13,93393314,93545246],
['A1A1',1,29270423,29420653],
['A1A1',16,81847145,81882509]
);

sub pollEvents{
    foreach my $cnv(@cnvs){
        my ($sample, $chrom, $start, $end) = @$cnv;
        my ($nLog2R, $iLN, $iLog2R, $iBFN, $iBF) = ('','','','','','','','');
        my $arrayTable = "ARRAY_$sample"."N_NIM_HG18";
        my $nmaTable = "NMA_$sample"."N";
#select N, to_char(log2R,'FM0.00') ||  ' (' || to_char(p,'FM0.0EEEE') || ')' log2R
        runSQL("
select to_char(log2R,'FM0.00') || ' (' || N || ') ' ||
       case when p < 0.00001 then '****'
            when p < 0.0001 then '***'
            when p < 0.001 then '**'
            when p < 0.01 then '*' end 
from (
select count(*) N,  Avg(log2R) log2R, StdDev(log2R) sd, STATS_T_TEST_ONE(log2R, (select testmean from $nmaTable where nmatype = 'RATIO_STATS')) p
from (
select log(2,ratio) log2R
from $arrayTable 
where chromosome = $chrom 
and position >= $start
and position <= $end))", 
\($nLog2R));
        fetchRow();
        if($sample eq 'A3A2'){ 
        
#select N, to_char(log2R,'FM0.00') || ' +/- ' || to_char(sd,'FM0.00') || ' (' || to_char(p,'FM0.0EEEE') || ')' log2R             
        runSQL("
select to_char(log2R,'FM0.00') || ' (' || N || ') ' ||
       case when p < 0.00001 then '****'
            when p < 0.0001 then '***'
            when p < 0.001 then '**'
            when p < 0.01 then '*' end 
from (
select count(*) N,  Avg(nlog2R) log2R, StdDev(nlog2R) sd, STATS_T_TEST_ONE(nlog2R, (select testmean from nma_7974_7975 where nmatype = 'RATIO_STATS_$chrom')) p
from (
select log(2,nullif(t.ratio,0)) tlog2R, log(2,nullif(t.r,0)/nullif(r.r,0)) nlog2R
from array_7975_ill_hg18 t, array_7974_ill_hg18 r
where t.chromosome = r.chromosome
and t.position = r.position
and t.chromosome = $chrom 
and t.position >= $start
and t.position <= $end
AND r.GCSCORE >= 0.15 AND r.R > 0.1 AND r.R < 3.0))", 
\($iLog2R));
        fetchRow();

        
#select N, to_char(bf,'FM0.00') || ' +/- ' || to_char(sd,'FM0.00') || ' (' || to_char(p,'FM0.0EEEE') || ')' bf
        runSQL("
select to_char(bf,'FM0.00') || ' (' || N || ') ' ||
       case when p < 0.00001 then '****'
            when p < 0.0001 then '***'
            when p < 0.001 then '**'
            when p < 0.01 then '*' end 
from (
select count(*) N,  Avg(zygosity) bf, StdDev(zygosity) sd, STATS_T_TEST_ONE(nbf, (select testmean from nma_7974_7975 where nmatype = 'BFREQUENCY_STATS_$chrom')) p
from (
select t.zygosity, log(2,abs((nullif(t.bf,0)/nullif(r.bf,0)) - 1) + 1) nbf
from array_7975_ill_hg18 t, array_7974_ill_hg18 r
where t.chromosome = r.chromosome
and t.position = r.position
and t.chromosome = $chrom 
and t.position >= $start
and t.position <= $end
AND r.GCSCORE >= 0.15 AND r.R > 0.1 AND r.R < 3.0
and r.zygosity < 0.75
))", 
\($iBF));
        fetchRow();
        }
    
        print "$sample\t$chrom\t$start\t$end\t$nLog2R\t$iLog2R\t$iBF\n";
        
    }

}




1;



