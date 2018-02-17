# -*- coding: utf-8 -*-
"""
Created on Wed May 25 04:20:00 CEST 2016

@authors: Juan C Entizne
@email: juancarlos.entizne01[at]estudiant.upf.edu

Modified by Juan L. Trincado
@email: juanluis.trincado[at].upf.edu

"""

import os
import logging
from lib.diff_tools import multiple_conditions_analysis
from argparse import *


description = \
    "Description:\n" + \
    "This tool calculates the significance to the change in mean PSI values between conditions, across multiple conditions.\n" \
    "The conditions are tested in a sequential order specified as input.\n"

parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=False)

parser.add_argument('-m', '--method',
                    dest="method",
                    action="store",
                    required=True,
                    choices=['empirical', 'classical'],
                    help="Method to test significance. Required.")

parser.add_argument('-p', '--psi',
                    dest="conds",
                    action="store",
                    nargs="+",
                    help="Path of the PSI files. PSI files and the transcript expression (TPM) files "
                         "must have the same order."
                         "The conditions files and the tpm files must have the same order.")

parser.add_argument('-e', '--tpm',
                    dest="tpms",
                    action="store",
                    nargs="+",
                    help="Path of the transcript expression (TPM) files. Conditions files and the transcript expression "
                         "(TPM) files must have the same order."
                         "The conditions files and the tpm files must have the same order.")

parser.add_argument('-i', '--input',
                    dest="iox",
                    action="store",
                    nargs=1,
                    default=None,
                    help="Input file with the event-transcripts equivalence (.ioe or .ioi format)")

parser.add_argument('-a', '--area',
                    dest="area",
                    action="store",
                    nargs=1,
                    type=int,
                    default=[1000],
                    help="Number indicating the number of points in the local area distribution. (default: 1000)")

parser.add_argument('-l', '--lower-bound',
                    dest="lower_bound",
                    action="store",
                    nargs=1,
                    type=float,
                    default=[0],
                    help="Lower-bound for the absolute delta PSI value to test for significance. (Default: 0).")

parser.add_argument('-pa', '--paired',
                    dest="paired",
                    action="store_true",
                    default=False,
                    help="Boolean. Indicates if replicates in conditions are paired. (Default: False).")

parser.add_argument('-gc', '--gene-correction',
                    dest="gene_cor",
                    action="store_true",
                    default=False,
                    help="Boolean. If True, SUPPA correct the p-values by gene. (Default: False).")

parser.add_argument('-al', '--alpha',
                    dest="alpha",
                    action="store",
                    nargs=1,
                    type=float,
                    default=[0.05],
                    help="Family-wise error rate to use for the multiple test correction. (Default: 0.05).")

parser.add_argument('-s', '--save_tpm_events',
                    action="store_true",
                    default=False,
                    help="Boolean. If True, the average log TPM of the events will be saved in an external file (Default: False).")

parser.add_argument('-c', '--combination',
                    action="store_true",
                    dest="seq",
                    default=False,
                    help="Boolean. If True, SUPPA perform the analysis between all the possible combinations of conditions (Default: False).")

parser.add_argument('-me', '--median',
                    dest="median",
                    action="store_true",
                    default=False,
                    help="Boolean. If True, SUPPA use the median to calculate the Delta PSI. (Default: False).")

parser.add_argument('-th', '--tpm-threshold',
                    dest="tpm_th",
                    action="store",
                    nargs=1,
                    type=float,
                    default=[0.0],
                    help="Minimum transcript average TPM value within-replicates and between-conditions to be included in the analysis. (Default: 1.0).")

def nan_threshold_type(x):
    x = float(x)
    if x < 0.0 and x > 1.0:
        raise ArgumentTypeError("nan_threshold should be a float number between 0 and 1")
    return x

parser.add_argument('-nan', '--nan-threshold',
                    dest="nan_th",
                    action="store",
                    nargs=1,
                    type=nan_threshold_type,
                    default=[0.0],
                    help="Percentage allowed of samples per condition with nan values for returning a DeltaPSI (Default: 0, no missing values allowed).")

parser.add_argument('-o', '--output',
                    dest="output",
                    action="store",
                    default=None,
                    help="Name of the output files.")

parser.add_argument("-mo", "--mode", default="INFO",
                    help="to choose from DEBUG, INFO, WARNING, ERROR and CRITICAL")


def create_path(lst):

    temp_lst = []
    for fl in lst:
        if not os.path.isabs(fl):
            fl_path = os.getcwd()+"/"+fl
            temp_lst.append(fl_path)
        else:
            temp_lst.append(fl)

    return temp_lst


def main():

    args = parser.parse_args()

    # Parsing arguments
    mode = "logging." + args.mode

    # Setting logging preferences
    logger = logging.getLogger(__name__)
    logger.setLevel(eval(mode))

    # Setting the level of the loggers in lib
    # setToolsLoggerLevel(mode)

    # Check if path is absolute, if not the program use the current working path
    cond_files = create_path(args.conds)
    expr_files = create_path(args.tpms)
    ioe_fl = create_path(args.iox)

    # Check extension of input file
    id_type = ioe_fl[0].split(".")[-1].strip("\n")
    if id_type != "ioe" and id_type != "ioi":
        logger.info("Invalid input file. Input file has to be either IOE or IOI format "
                    "it must present the appropriate suffix.")
        exit(1)

    #multiple_conditions_analysis(args.method, cond_files, expr_files, ioe_fl[0], args.area[0], args.lower_bound[0],
     #                            args.paired, args.gene_cor, args.alpha[0], args.output)

    multiple_conditions_analysis(args.method, cond_files, expr_files, ioe_fl[0], args.area[0],
                                 args.lower_bound[0], args.paired, args.gene_cor, args.alpha[0],
                                 args.save_tpm_events, args.seq, args.median, args.tpm_th[0],
                                 args.nan_th[0],args.output)


if __name__ == "__main__":
    main()
