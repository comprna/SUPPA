#!/usr/bin/perl -w

use strict;

# This program is distributed under the MIT license
# eduardo.eyras@upf.edu

#perl rmats_to_suppa_ids_SE_events.pl  SE.MATS.JunctionCountOnly.txt              > rmats_se_junction.txt
#perl rmats_to_suppa_ids_SE_events.pl  SE.MATS.ReadsOnTargetAndJunctionCounts.txt > rmats_se_junction_and_reads.txt

my ($file) = @ARGV;
unless($file){
    print STDERR "$0 file\n";
    print STDERR "where file: SE.MATS.JunctionCountOnly.txt  or SE.MATS.ReadsOnTargetAndJunctionCounts.txt\n";
    exit(0);
}

while(<>){
    chomp;
#ID      GeneID  geneSymbol      chr     strand  exonStart_0base exonEnd upstreamES      upstreamEE      downstreamES    downstreamEE    ID      IJC_SAMPLE_1    SJC_SAMPLE_1    IJC_SAMPLE_2    SJC_SAMPLE_2    IncFormLen      SkipFormLen     PValue  FDR     IncLevel1       IncLevel2       IncLevelDifference 10012   "Dync1i2"       NA      chr2    +       71233648        71233708        71228599        71228690        71235914        71236030        10012   353,700,245     60,125,47       0,0,0   111,121,129     104     84      0.0     0.0     0.826,0.819,0.808       0.0,0.0,0.0     0.818

    next if /GeneID/;
    my ($ID, $GeneID, $geneSymbol, $chr, $strand, 
	$exonStart_0base, $exonEnd, $upstreamES, $upstreamEE, $downstreamES, $downstreamEE, 
	$ID2, $IC_SAMPLE_1, $SC_SAMPLE_1, $IC_SAMPLE_2, $SC_SAMPLE_2, $IncFormLen, $SkipFormLen, $PValue, $FDR, $IncLevel1, $IncLevel2, $IncLevelDifference) = split;

    # we convert it to SUPPA coordindates:
    # ENSG00000103126;SE:chr16:339607-341190:341297-343488:-  0.2102558902    0.027972028     0.28
    
    $GeneID=~/\"(\S+)\"/;
    my $gene_id = $1;
    
    # convert the info to a suppa ID (in closed-coordinates = 1-based)
    my $suppa_id = $gene_id.";SE:".$chr.":".$upstreamEE."-".($exonStart_0base+1).":".$exonEnd."-".($downstreamES+1).":".$strand;
    
    my $s = join "\t", ($suppa_id, $IncLevelDifference, $PValue, $FDR); 
    print $s."\n";
}
    
