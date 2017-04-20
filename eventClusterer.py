# -*- coding: utf-8 -*-
"""
Created on Wed May 25 04:20:00 CEST 2016

@authors: Juan C Entizne, Juan L. Trincado
@email: juancarlos.entizne01[at]estudiant.upf.edu,
        juanluis.trincado[at]upf.edu
"""

import os
import logging
from lib.cluster_tools import cluster_analysis
from argparse import ArgumentParser, RawTextHelpFormatter


description = \
    "Description:\n\n" + \
    "This tool cluster events that change significantly in at least one pair of conditions, across multiple conditions.\n" + \
    "This tool takes as input the .dpsi and .psivec files generate by SUPPA differentialAnalysis method\n" + \
    "and generates a .clustvec on which the events has been tagged according to their cluster membership\n"


parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=False)

parser.add_argument('-d', '--dpsi',
                    dest="dpsi",
                    nargs=1,
                    action="store",
                    help="Input file of delta-PSI values (.dpsi format)")

parser.add_argument('-p', '--psivec',
                    dest="psivec",
                    nargs=1,
                    action="store",
                    help="Input file with PSI values (.psivec format)")

parser.add_argument('-st', '--sig-threshold',
                    dest="sig_threshold",
                    action="store",
                    type=float,
                    default=0.05,
                    help="P-value cut-off for significant events. (Default: 0.05).")

parser.add_argument('-dt', '--dpsi-threshold',
                    dest="dpsi_threshold",
                    action="store",
                    type=float,
                    default=0.05,
                    help="Lower-bound for the absolute delta PSI value to cluster. (Default: 0.05).")

parser.add_argument('-e', '--eps',
                    dest="eps",
                    action="store",
                    type=float,
                    default=0.05,
                    help="Maximum (Euclidean) distance (between 0 and 1) to consider two events as members of "
                         "the same cluster. (Default: 0.05).")

parser.add_argument('-s', '--separation',
                    action="store",
                    type=float,
                    default=0,
                    help="Minimum separation for considering two points in different clusters. (Default: 0).")

parser.add_argument('-n', '--min-pts',
                    dest="minpts",
                    action="store",
                    type=int,
                    default=20,
                    help="Minimum number of events required per cluster. (Default: 20).")

parser.add_argument("-m", "--metric", dest="metric", choices=["euclidean", "manhattan", "cosine"],
                    default="euclidean", help="Distance function to be used."
                                                                              "Options:\n"
                                                                              "\teuclidean (Default),\n"
                                                                              "\tmanhattan,\n"
                                                                              "\tcosine.\n")
parser.add_argument("-c", "--clustering", choices=["OPTICS", "DBSCAN"],
                    default="DBSCAN", help="Clustering method to use."
                                                                              "Options:\n"
                                                                              "\tOPTICS ,\n"
                                                                              "\tDBSCAN (Default).\n")
parser.add_argument('-g', '--groups',
                    dest="indexes",
                    action="store",
                    required=True,
                    type=str,
                    nargs="*",
                    help="Ranges of column numbers specifying the replicates per condition. "
                         "Column numbers have to be continuous, with no overlapping or missing columns between them. "
                         "Ex: 1-3,4-6")

parser.add_argument('-o', '--output',
                    dest="output",
                    action="store",
                    help="Name of the output file.")


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

    # Check if path is absolute, if not the program use the current working path
    dpsi_file = create_path(args.dpsi)
    psivec_file = create_path(args.psivec)

    cluster_analysis(dpsi_file[0], psivec_file[0], args.sig_threshold, args.dpsi_threshold, args.eps, args.minpts,
                     args.metric,args.indexes[0], args.clustering, args.separation, args.output)

if __name__ == "__main__":
    main()
