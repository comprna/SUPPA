# -*- coding: utf-8 -*-
"""
Created on Fri May 23 10:17:33 2014

@author: Miha Skalic
@email: miha.skalic[at]gmail.com
"""


import sys
import logging
from argparse import ArgumentParser, RawTextHelpFormatter
from lib.tools import *
from lib.gtf_store import *


description = \
    "Description:\n\n" + \
    "This tool calculates the PSI (Percentatge Splice In) for the different\n" + \
    "transcripts of a gene.\n" + \
    "It reads a gtf to get transcript-gene relationship and an expression file\n" + \
    "of the different transcripts\n" 
    
parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=False)
parser.add_argument("-g", "--gtf-file", help="Input gtf file",
                    required=True)
parser.add_argument("-e", "--expression-file", required=True,
                    help="Input expression file")
parser.add_argument("-o", "--output-file", required=True,
                    help="Path and name of the ouput file")
parser.add_argument("-m", "--mode", default="INFO",
                    help="to choose from DEBUG, INFO, WARNING, ERROR and CRITICAL")


def expression_reader(exp_file):
    """
    Reads in expression file and returns dict
    of transcript expressions and first line.
    """
    if not os.path.isfile(exp_file):
        sys.stderr.write("Expression file does not exist. Quiting\n")
        exit(1)

    expressions = {}
    with open(exp_file, 'r') as handle:
        first_line = nextel(handle).strip()
        for line in handle:
            line = line.strip().split('\t')
            expressions[line[0]] = [float(xp) for xp in line[1:]]
    return expressions, first_line


def expression_writer(genomeinfo, expressions, firstline, output_file):
    """
    Function to write perIsoform inclusion
    """
    output_file += '_isoform.psi'
    entriesnumber = len(expressions[nextel(expressions.__iter__())])
    with open(output_file, 'w') as handle:
        handle.write(firstline + '\n')
        for gene, _, _ in genomeinfo:
            expr_sum = [0 for _ in range(entriesnumber)]

            # collect expression
            for transcript in gene.sortedTranscripts:
                if transcript not in expressions:
                    logger.info(('Expression for transcript "{}" not found. '
                                 'Ignoring it in calculation.').format(transcript))
                else:
                    expr_sum = list(map(lambda exp_pair: exp_pair[0] + exp_pair[1], zip(expr_sum, expressions[transcript])))

            # calculate expression
            if 0 in expr_sum:
                logger.debug('Gene "{}" has at least one replicate with 0 expression.'.format(gene.name))
                expr_sum = [y if y else float('NaN') for y in expr_sum]

            for transcript in gene.sortedTranscripts:
                if transcript not in expressions:
                    continue
                t_exp = map(lambda exp_pair: exp_pair[1] / exp_pair[0], zip(expr_sum, expressions[transcript]))
                handle.write('{};{}\t{}\n'.format(gene.name, transcript,
                                                  '\t'.join([str(exp_val) for exp_val in t_exp])))


def main(): 
    args = parser.parse_args()

    #Parsing arguments
    mode = "logging." + args.mode

    #Setting logging preferences
    logger = logging.getLogger(__name__)
    logger.setLevel(eval(mode))

    #Setting the level of the loggers in lib
    setToolsLoggerLevel(mode)

    #PREPAIRING GTF
    my_genome = Genome()
    logger.info("Reading GTF data.")
    fetched_exons = gtf_reader(args.gtf_file, logger)

    # Check for empy sequences
    if len(fetched_exons) == 0:
        logger.info("No exons found. Check format and content of your GTF file.")
        exit(1)

    for exon_meta in fetched_exons:
        my_genome.add_to_genes(exon_meta)

    # split non overlapping genes
    my_genome.sort_transcripts()
    my_genome.split_genes()

    logger.info("Reading Expression data.")
    trans_expres, sample_names = expression_reader(args.expression_file)
    if not trans_expres:
        logger.info("No expressions found. Check format and content of your expression file.")
        exit(1)

    # Calculate and write output
    logger.info("Calculating inclusion and generating output.")
    expression_writer(my_genome, trans_expres, sample_names, args.output_file)


if __name__ == '__main__':
    main()
