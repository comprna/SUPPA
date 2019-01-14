#!/usr/bin/perl -w

# This program is distributed under the MIT license
# eduardo.eyras@upf.edu

use strict;
#ID      GeneID  geneSymbol      chr     strand  longExonStart_0base     longExonEnd     shortES shortEE flankingES      flankingEE      ID      IC_SAMPLE_1     SC_SAMPLE_1     IC_SAMPLE_2     SC_SAMPLE_2     IncFormLen      SkipFormLen     PValue  FDR     IncLevel1       IncLevel2       IncLevelDifference
#6342    "ENSG00000028203"       "VEZT"  chr12   +       95692869        95694948        95693940        95694948        95689826        95690034        6342    94,75,73        37,41,37        3,3,7   62,80,81        1143    86      0.0     0.0     0.16,0.121,0.129        0.004,0.003,0.006       0.132

my $count = 0;
while(<>){
    chomp;
    next if /GeneID/;
    $count++;
    my ($ID, $GeneID, $geneSymbol, $chr, $strand, 
	$longExonStart_0base, $longExonEnd, $shortES, $shortEE, $flankingES, $flankingEE, 
	$ID2, $IC_SAMPLE_1, $SC_SAMPLE_1, $IC_SAMPLE_2, $SC_SAMPLE_2, $IncFormLen, $SkipFormLen, $PValue, $FDR, $IncLevel1, $IncLevel2, $IncLevelDifference) = split;
    
    # we convert it to SUPPA coordindates:
    # For A5 and A3 the PSI and dPSI is for the long form
    # for forward we give first the fixed boundary and then the long form of the exon
    #ENSG00000002586;A3:chrX:2644414-2655922:2644414-2656241:+   0.0006483856        0.3746253746   
    # for reverse we give first the long form and then the fixed postion
    #ENSG00000002016;A3:chr12:1025701-1025805:1025698-1025805:-  nan       1.0
    
    $GeneID=~/\"(ENSG.*)\"/;
    my $gene_id = $1;

    my $suppa_id;
    if ($strand eq "+"){
	$suppa_id = $gene_id.";A3:".$chr.":".($flankingEE)."-".($longExonStart_0base+1).":".($flankingEE)."-".($shortES+1).":".$strand;
    }
    else{
	$suppa_id = $gene_id.";A3:".$chr.":".$longExonEnd."-".($flankingES+1).":".$shortEE."-". ($flankingES+1).":".$strand;
    }
    
    my $s = join "\t", ($suppa_id, $IncLevelDifference, $FDR); 
    print $s."\n";
}
    
#print $count."\n";
