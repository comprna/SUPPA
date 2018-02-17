# -*- coding: utf-8 -*-
"""
Created on Tue May  6 12:19:45 2014

@author: Gael P Alamancos
@email: gael.perez[at]upf.edu
"""

import sys
import logging
import numpy as np
from argparse import ArgumentParser, RawTextHelpFormatter
from lib.tools import *


description = \
    "Description:\n\n" + \
    "This tool reads an ioe file and a transcript expression file and calculates\n" + \
    "the Percentage of Sliced In (PSI)\n"
    
parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=False)
parser.add_argument("-i", "--ioe-file", help="Input file with the event-transcripts equivalence (.ioe format).", required=True)
parser.add_argument("-e", "--expression-file", required=True,
                    help="Input transcript expression file.")
parser.add_argument("-o", "--output-file", required=True,
                    help="Output psi file.")
parser.add_argument("-f", "--total-filter", type=float, default=0,
                    help="Minimum total expression of the transcripts involved in the event.")
parser.add_argument('-s', '--save_tpm_events',
                    action="store_true",
                    default=False,
                    help="Boolean. If True, save the TPM of the events in an external file (Default: False).")
parser.add_argument("-m", "--mode", default="INFO",
                    help="to choose from DEBUG, INFO, WARNING, ERROR and CRITICAL")


def main():

    args = parser.parse_args()

    #Parsing arguments
    mode = "logging." + args.mode

    #Setting logging preferences
    logger = logging.getLogger(__name__)
    logger.setLevel(eval(mode))

    # Setting the level of the loggers in lib
    setToolsLoggerLevel(mode)

    output_file = args.output_file
    total_filter = args.total_filter
    expression_file = args.expression_file

    expression_dictionary = {}  # to store transcript expression [transcript_id] = {[colId] = expression}
    psi_dictionary = {}  # to store the psi values calculated [event_id] = {[colId] = expression}
    expression_dictionary_events = {}  # to store the expression values calculated [event_id] = {[colId] = expression}
    col_ids = []  # to store all the column_id for the expression fields
    try:
        #Buffering EXPRESSION FILE
        factory = FactoryReader()
        r = factory.getReader("expression")
        r.openFile(expression_file)
        logger.info("Buffering transcript expression levels.")
        line = r.readLine()
        try:
            while True:
                try:
                    arguments = nextel(line)
                    col_ids = arguments["colIds"]
                    expression_dictionary[arguments["transcript_id"]] = {}
                    # Assign per each colId the expression [colId] = expression
                    for key, value in arguments["expression"].items():
                        expression_dictionary[arguments["transcript_id"]][key] = float(value)
                except ValueError:
                    logger.error("%s expression is not a float. Skipping..." % arguments["transcript_id"])
                    continue
        except StopIteration:
            if not expression_dictionary:
                logger.error("No expression values have been buffered.")
                sys.exit(1)
                
            #Buffering IOE and calculating PSI
            factory = FactoryReader()
            r = factory.getReader("ioe")
            # r.openFile(args.ioe_file)
            r.openFile(args.ioe_file)
            logger.info("Calculating PSI from the ioe file.") 
            line = r.readLine()
            try:
                while True:
                    alternative_transcripts = {}
                    total_transcripts = {}
                    arguments = nextel(line)
                    if psi_dictionary.get(arguments["event_id"]):
                        logger.error("Duplicated event %s. Skipping line..." %
                                     arguments["event_id"])
                        continue
                    skip = False  # to avoid checking total_iso if event == "NA"
                    psi_dictionary[arguments["event_id"]] = {}
                    expression_dictionary_events[arguments["event_id"]] = {}
                    for x in col_ids:
                        #set to 0 all the cumulative expression values
                        alternative_transcripts[x] = 0
                        total_transcripts[x] = 0
                    #Add all the expression values of the alternative transcripts
                    for tr in arguments["alt_iso"].rstrip("\n").split(","):
                        try:
                            for x in col_ids:
                                alternative_transcripts[x] += \
                                    expression_dictionary[tr][x]
                        except KeyError:
                            logger.error(
                                "transcript %s not found in the \"expression file\"." %
                                tr)
                            logger.error("PSI not calculated for event %s." %
                                         arguments["event_id"])
                            for x in col_ids:
                                psi_dictionary[arguments["event_id"]][x] = "NA"
                                expression_dictionary_events[arguments["event_id"]][x] = "NA"
                            skip = True
                    # Add all the expression values of the total transcripts
                    if not skip:
                        for tr in arguments["total_iso"].rstrip("\n").split(","):
                            try:
                                for x in col_ids:
                                    total_transcripts[x] += \
                                        expression_dictionary[tr][x]
                            except KeyError:
                                logger.error(
                                    "transcript %s not found in the \"expression file\"." %
                                    tr)
                                logger.error("PSI not calculated for event %s." %
                                             arguments["event_id"])
                                for x in col_ids:
                                    psi_dictionary[arguments["event_id"]][x] = "NA"
                                    expression_dictionary_events[arguments["event_id"]][x] = "NA"
                                skip = True
                    for x in col_ids:
                        #If it passes the filter and skip == False
                        if total_transcripts[x] >= total_filter and not skip:
                            try:
                                psi_dictionary[arguments["event_id"]][x] = (
                                    alternative_transcripts[x]/total_transcripts[x])
                                expression_dictionary_events[arguments["event_id"]][x] = total_transcripts[x]
                            except ZeroDivisionError:
                                logger.debug("Zero division for event %s.(psi= nan)." % arguments["event_id"])
                                psi_dictionary[arguments["event_id"]][x] = np.nan
                                expression_dictionary_events[arguments["event_id"]][x] = total_transcripts[x]
                        #If it passes the filter but skip == True
                        elif not skip:   
                            #for x in col_ids:
                            psi_dictionary[arguments["event_id"]][x] = np.nan
                            expression_dictionary_events[arguments["event_id"]][x] = total_transcripts[x]

            except StopIteration:
                writer = Writer.getWriter("PSI")
                logger.info("Generating output %s" % (output_file + ".psi"))
                writer.openFile(output_file)
                writer.writeLine("\t".join(col_ids), False)
                for key, value in sorted(psi_dictionary.items()):
                    logger.debug("Calculating psi for %s" % key)
                    psi_line = PsiWriter.lineGenerator(key, value, col_ids)
                    writer.writeLine("\t".join(psi_line))
                writer.closeFile()

                if (args.save_tpm_events == True):
                    #Save the expression of the events
                    writer = Writer.getWriter("TPM")
                    logger.info("Generating output %s" % (output_file + ".tpm"))
                    writer.openFile(output_file)
                    writer.writeLine("\t".join(col_ids), False)
                    for key, value in sorted(expression_dictionary_events.items()):
                        logger.debug("Calculating tpm for %s" % key)
                        tpm_line = TpmWriter.lineGenerator(key, value, col_ids)
                        writer.writeLine("\t".join(tpm_line))
                    writer.closeFile()

    except BaseException:
        logger.error("Unknown error: %s" % sys.exc_info()[1].args[0])
        sys.exit(1)
    logger.info("Done")

if __name__ == '__main__':
    main()