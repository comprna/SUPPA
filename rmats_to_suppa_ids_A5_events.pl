#!/usr/bin/perl -w

# This program is distributed under the MIT license
# eduardo.eyras@upf.edu

use strict;
#ID      GeneID  geneSymbol      chr     strand  longExonStart_0base     longExonEnd     shortES shortEE flankingES      flankingEE      ID      IC_SAMPLE_1     SC_SAMPLE_1     IC_SAMPLE_2     SC_SAMPLE_2     IncFormLen      SkipFormLen     PValue  FDR     IncLevel1       IncLevel2       IncLevelDifference
#4071    "ENSG00000187514"       "PTMA"  chr2    +       232576057       232576693       232576057       232576129       232577152       232577226       4071    8670,9224,9315  609,821,830     10645,11732,13398       23,25,41        596     46      0.0     0.0     0.524,0.464,0.464       0.973,0.973,0.962       -0.485

#ID      GeneID  geneSymbol      chr     strand  exonStart_0base exonEnd upstreamES      upstreamEE      downstreamES    downstreamEE    ID      IC_SAMPLE_1     SC_SAMPLE_1     IC_SAMPLE_2     SC_SAMPLE_2     IncFormLen      SkipFormLen     PValue  FDR     IncLevel1       IncLevel2       IncLevelDifference
#10054   "ENSG00000067836"       "ROGDI" chr16   -       4851267 4851322 4850498 4850579 4851503 4851586 10054   857,112,161     77,7,14 772,101,142     1,0,0   77      36      0.0     0.0     0.839,0.882,0.843       0.997,1.0,1.0   -0.144

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
    # for forward we give first the long form of the exon
    # ENSG00000004766;A5:chr7:92869263-92881966:92869247-92881966:+         0.0       0.3541458541
    # for reverse we give first the long form
    # ENSG00000004700;A5:chr12:21654306-21654397:21654306-21654415:-        0.0157526072        0.2602397602
    
    $GeneID=~/\"(ENSG.*)\"/;
    my $gene_id = $1;

    my $suppa_id;
    if ($strand eq "+"){
	$suppa_id = $gene_id.";A5:".$chr.":".$longExonEnd."-".($flankingES + 1).":".$shortEE."-".($flankingES + 1).":".$strand;
    }
    else{
	$suppa_id = $gene_id.";A5:".$chr.":".$flankingEE."-".($longExonStart_0base+1).":".$flankingEE."-". ($shortES+1).":".$strand;
    }
    
    my $s = join "\t", ($suppa_id, $IncLevelDifference, $FDR); 
    print $s."\n";
}
    
#print $count."\n";
