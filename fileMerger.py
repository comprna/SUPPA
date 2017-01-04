# -*- coding: utf-8 -*-
"""
Created on Wed May 25 04:20:00 CEST 2016

@authors: Juan C Entizne
@email: juancarlos.entizne01[at]estudiant.upf.edu
"""

import os
import argparse
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter


description = \
    "Description:\n" + \
    "This tool joins multiple psi or expression files together.\n" \
    "This tool assume that the first field (column) of the files are in common.\n"

parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=False)

parser.add_argument('-i', '--input-files',
                    dest="input",
                    action="store",
                    nargs="+",
                    help="Space separated list of the files to be joined. "
                         "If the absolute path is not indicate the program use the current working directory instead.")

parser.add_argument('-f', '--file-extension',
                    dest="ext",
                    action="store",
                    required=True,
                    choices=['psi', 'tpm'],
                    help="Extension of the output file. Required.")

parser.add_argument('-o', '--output',
                    dest="output",
                    action="store",
                    default=None,
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


def merge_files(fl_lst, output, ext):


    df_lst = []
    for fl in fl_lst:
        df = pd.read_table(fl, sep='\t', index_col=0, header=0)

        old_header = df.columns.values
        new_header = [os.path.basename(fl).split(".")[0]+"_"+col_id for col_id in old_header]
        df.rename(columns=dict(zip(old_header, new_header)), inplace=True)

        df_lst.append(df)

    merged_dfs = pd.concat(df_lst, axis=1)

    header = merged_dfs.columns.values

    with open("%s.%s" % (output, ext), "w+") as fh:
            ln = "\t".join(header)
            fh.write(ln+"\n")

    with open("%s.%s" % (output, ext), "a") as fh:
            merged_dfs.to_csv(fh, sep="\t", na_rep="nan", header=False)


def main():

    args = parser.parse_args()

    input_files = args.input
    outname = args.output
    fl_ext = args.ext

    fl_lst = create_path(input_files)

    merge_files(fl_lst, outname, fl_ext)

if __name__ == "__main__":
    main()
