#!/usr/bin/perl -w

# This program is distributed under the MIT license
# eduardo.eyras@upf.edu

use strict;

#ID      GeneID  geneSymbol      chr     strand  exonStart_0base exonEnd upstreamES      upstreamEE      downstreamES    downstreamEE    ID      IC_SAMPLE_1     SC_SAMPLE_1     IC_SAMPLE_2     SC_SAMPLE_2     IncFormLen      SkipFormLen     PValue  FDR     IncLevel1       IncLevel2       IncLevelDifference
#10054   "ENSG00000067836"       "ROGDI" chr16   -       4851267 4851322 4850498 4850579 4851503 4851586 10054   857,112,161     77,7,14 772,101,142     1,0,0   77      36      0.0     0.0     0.839,0.882,0.843       0.997,1.0,1.0   -0.144

my $count = 0;
while(<>){
    chomp;
    next if /GeneID/;
    $count++;
    my ($ID, $GeneID, $geneSymbol, $chr, $strand, 
	$exonStart_0base, $exonEnd, $upstreamES, $upstreamEE, $downstreamES, $downstreamEE, 
	$ID2, $IC_SAMPLE_1, $SC_SAMPLE_1, $IC_SAMPLE_2, $SC_SAMPLE_2, $IncFormLen, $SkipFormLen, $PValue, $FDR, $IncLevel1, $IncLevel2, $IncLevelDifference) = split;

    # we convert it to SUPPA coordindates:
    # ENSG00000103126;SE:chr16:339607-341190:341297-343488:-  0.2102558902    0.027972028     0.28
    
    $GeneID=~/\"(ENSG.*)\"/;
    my $gene_id = $1;

    my $suppa_id = $gene_id.";SE:".$chr.":".$upstreamEE."-".($exonStart_0base+1).":".$exonEnd."-".($downstreamES+1).":".$strand;

    my $s = join "\t", ($suppa_id, $IncLevelDifference, $FDR); 
    print $s."\n";
}
    
#print $count."\n";
