# -*- coding: utf-8 -*-
"""
Created on Wed Aug 06 17:51:05 2014

@author: Gael P Alamancos
@email: gael.perez[at]upf.edu
"""

import fileMerger as joinFiles
import psiPerGene as psiPerIsoform
import psiCalculator as psiPerEvent
import eventGenerator as generateEvents
import eventClusterer as clusterAnalysis
import significanceCalculator as diffSplice
import logging
import argparse 
import sys


description = "Description:\n\n" + \
              "SUPPA allows you to generate all the possible Alternative Splicing events from an annotation file, \n" \
              "calculate the PSI values per event, calculate differential splicing across multiple conditions \n" \
              "with replicates, and cluster events across conditions \n" \
              "For further information, see the help of each subcommand."

parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
subparsers = parser.add_subparsers()

# EventGenerator parser
eventGeneratorSubparser = subparsers.add_parser(
    "generateEvents", parents=[generateEvents.parser],
    help="Calculates the Alternative Splicing events from an annotation file.")
eventGeneratorSubparser.set_defaults(which="generateEvents")

# psiCalculator parser
psiCalculatorSubparser = subparsers.add_parser(
    "psiPerEvent", parents=[psiPerEvent.parser],
    help="Calculates the PSI value for each event previously generated.")
psiCalculatorSubparser.set_defaults(which="psiPerEvent")

# psiPerGene parser
psiPerGeneSubparser = subparsers.add_parser(
    "psiPerIsoform", parents=[psiPerIsoform.parser],
    help="Calculates the PSI value for each isoform.")
psiPerGeneSubparser.set_defaults(which="psiPerIsoform")

# significanceCalculator parser
significanceCalculatorSubparser = subparsers.add_parser(
    "diffSplice", parents=[diffSplice.parser],
    help="Calculates differentially spliced events across multiple conditions.")
significanceCalculatorSubparser.set_defaults(which="diffSplice")

# eventClusterer parser
eventClustererSubparser = subparsers.add_parser(
    "clusterEvents", parents=[clusterAnalysis.parser],
    help="Calculates clusters of events across conditions.")
eventClustererSubparser.set_defaults(which="clusterEvents")

# fileMerger parser
fileMergerSubparser = subparsers.add_parser(
    "joinFiles", parents=[joinFiles.parser],
    help="Join multiple tab separated files into a single file.")
fileMergerSubparser.set_defaults(which="joinFiles")



# Setting logging preferences
logger = logging.getLogger(__name__)


def main():

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    try:
        args = parser.parse_args()
        if args.which == "generateEvents":
            generateEvents.parser = parser  # Setting the module aparser
            generateEvents.main()
        elif args.which == "psiPerEvent":
            psiPerEvent.parser = parser  # Setting the module parser
            psiPerEvent.main()
        elif args.which == "psiPerIsoform":
            psiPerIsoform.parser = parser  # Setting the module parser
            psiPerIsoform.main()
        elif args.which == "diffSplice":
            diffSplice.parser = parser  # Setting the module parser
            diffSplice.main()
        elif args.which == "clusterEvents":
            clusterAnalysis.parser = parser  # Setting the module parser
            clusterAnalysis.main()
        elif args.which == "joinFiles":
            joinFiles.parser = parser  # Setting the module parser
            joinFiles.main()
    except Exception:
        logger.error("Unknown error: {}".format(sys.exc_info()))
        sys.exit(1)
        
if __name__ == '__main__':
    main()
