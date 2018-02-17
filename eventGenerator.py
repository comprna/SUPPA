# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 14:57:05 2014

@authors: Gael P Alamancos and Miha Skalic
@email: gael.perez[at]upf.edu, miha.skalic[at]gmail.com
"""

import sys
import logging
from argparse import ArgumentParser, RawTextHelpFormatter
from lib.tools import *
from lib.gtf_store import *
from lib.event import *

# Setting argument parser
# parser = argparse.ArgumentParser()
description = \
    "Description:\n\n" + \
    "This tool reads an annotation file and generates different alternative\n" + \
    "splicing(AS) events depending on the user's choice.\n" + \
    "It outputs a \"gtf\" file to load on genome browser to visualize the \n" + \
    "events as well as an \"ioe\" file to use for the further calculation of\n" + \
    "the PSI values of each possible event conformation."
parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=False)
parser.add_argument("-i", "--input-file", help="specify input file",
                    required=True)
parser.add_argument("-o", "--output-file", help="specify output path and" +
                    " name without extension", required=True)

parser.add_argument("-e", "--event-type", nargs='+', required=False, choices=["SE", "SS", "MX", "RI", "FL"],
                    help="list of events to analyze. "
                    "(space separated)\n\n"
                    "Options:\n"
                    "\tSE -- Skipping Exon\n"
                    "\tSS -- Alternative Splice Site (5'/3')\n"
                    "\tMX -- Mutually Exclusive Exon\n"
                    "\tRI -- Retained Intron\n"
                    "\tFL -- Alternative First/Last Exon\n")
parser.add_argument("-b", "--boundary", choices=["S", "V"], default="S", help="Boundary type."
                    "Options:\n"
                    "\tS -- Strict (Default)\n"
                    "\tV -- Variable\n")
parser.add_argument("-t", "--threshold", default=10, type=int,
                    help="Variability treshold. In case of strict boundaries this argument is ignored" +
                    "(Default: 10nt).")
parser.add_argument("-p", "--pool-genes", action="store_true",
                    help="pool together overlapping genes.")
parser.add_argument("-l", "--exon-length", default=100, type=int,
                    help="length of the exons for its visualization. " +
                    "(Default: 100nt)")
parser.add_argument("-m", "--mode", default="INFO",
                    help="to choose from DEBUG, INFO, WARNING, ERROR and CRITICAL")
parser.add_argument('-f', '--format',
                    dest="format",
                    action="store",
                    required=True,
                    choices=['ioe', 'ioi'],
                    help="Format of the annotation file. Required.")


def main():

    args = parser.parse_args()

    # Parsing arguments
    mode = "logging." + args.mode

    # Setting logging preferences
    logger = logging.getLogger(__name__)
    logger.setLevel(eval(mode))

    # Setting the level of the loggers in lib
    setToolsLoggerLevel(mode)

    exon_length = int(args.exon_length)

    my_genome = Genome()

    # Check if event list exist only if format is IOE
    if args.format == "ioe" and args.event_type == None:
        logger.info("No event type found. It is necessary to determine at least a type of event to create an IOE file.")
        exit(1)

    logger.info("Reading input data.")
    fetched_exons = gtf_reader(args.input_file, logger)

    # Check for empy sequence
    if len(fetched_exons) == 0:
        logger.info("No exons found. Check format and content of your GTF file.")
        exit(1)

    for exon_meta in fetched_exons:
        my_genome.add_to_genes(exon_meta)

    my_genome.sort_transcripts()

    if args.pool_genes:
        my_genome.split_genes()
        logger.info("Pooling genes")
        my_genome.poll_genes()

    if args.format == "ioi":
        ioi_writer(my_genome, args.output_file)
    else:
        make_events(args.event_type, my_genome, args.input_file, args.output_file, exon_length, logger,
                b_type=args.boundary, th=args.threshold)
    logger.info("Done")
    
if __name__ == '__main__':
    main()
