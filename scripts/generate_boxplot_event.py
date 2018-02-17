# The next script will format a phenotype table (junctions, events, trasncripts...)
# for runnning FastQTL analysis

#This version is for formatting the SCLC phenotype

"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

generate_boxplot_event.py: Generates a boxplot with the PSI values, given which samples are in which conditions
"""

import sys
import logging
import matplotlib.pyplot as plt
import numpy as np
import re



from argparse import ArgumentParser, RawTextHelpFormatter

description = \
    "Description:\n\n" + \
    "This script accept a phenotype table (junctions, events, transcripts...)\n" + \
    "and a genotype table (mutations associated to K-mers or SMRs) and returns a formatted table\n" + \
    "for using with FastQTL"

parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=True)
parser.add_argument("-i", "--input", required=True,
                    help="Input file")
parser.add_argument("-e", "--event", required=True, type=str,
                    help="Event to plot")
parser.add_argument('-g', '--groups',
                    action="store",
                    required=True,
                    type=str,
                    nargs="*",
                    help="Ranges of column numbers specifying the replicates per condition. "
                         "Column numbers have to be continuous, with no overlapping or missing columns between them. "
                         "Ex: 1-3,4-6")
parser.add_argument('-c', '--conds',
                    action="store",
                    required=False,
                    default="0",
                    type=str,
                    nargs="*",
                    help="Name of each one of the conditions. Ex: Mutated,Non_mutated")
parser.add_argument("-o", "--output", required=True,
                    help="Output path")

# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# create console handler and set level to info
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)

def main():

    args = parser.parse_args()

    input_file = args.input
    event = args.event
    groups = re.findall(r"[\w]+", args.groups[0])
    output_path = args.output

    # input_file = "/home/juanluis/Desktop/Work/Master_class/events.psi"
    # event = "ENSG00000149554;SE:chr11:125496728-125497502:125497725-125499127:+"
    # groups = ['1','3','4','6']
    # output_path = "/home/juanluis/Desktop/Work/Master_class/"

    try:

        logger.info("Reading input file...")

        dict_PSI = {}
        cond = 1
        success = False
        file = open(input_file)
        for line in file:
            tokens = line.rstrip().split("\t")
            if (tokens[0]==event):
                success = True
                for i,x in enumerate(groups):
                    if(i%2==1):
                        continue
                    PSI = []
                    samples = range(int(groups[i]),int(groups[i+1])+1)
                    #Get the PSI of this group of samples
                    for j in samples:
                        PSI.append(tokens[j])
                    dict_PSI[cond] = PSI
                    cond = cond + 1
                break

        if(success):
            #Create the boxplot
            data_to_plot = []
            for key in dict_PSI.keys():
                data_to_plot.append(list(map(float,dict_PSI[key])))
            # Create a figure instance
            fig = plt.figure(figsize=(9, 6))
            # Create an axes instance
            ax = fig.add_subplot(111)
            # Create the boxplot
            bp = ax.boxplot(data_to_plot, patch_artist=True, sym='')
            # change the style of fliers and their fill
            for flier in bp['fliers']:
                flier.set(marker='.', color='#000000', alpha=0.7)
            # Assign different colors
            colors = ['lightblue', 'pink']
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
            for j in range(len(data_to_plot)):
                y = data_to_plot[j]
                x = np.random.normal(1 + j, 0.02, size=len(y))
                plt.plot(x, y, 'ko', alpha=0.5)
            # Custom x-axis labels if the user has input conditions
            if (args.conds != "0"):
                conditions = re.findall(r"[\w]+", args.conds[0])
                ax.set_xticklabels(conditions)
            # Leave just ticks in the bottom
            ax.get_xaxis().tick_bottom()
            ax.set_ylabel('PSI')
            # Set the title
            title = "Event: " + event
            ax.set_title(title, fontsize=10)
            # Add a horizontal grid to the plot,
            ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
            # Set the limits for the y axes
            ax.set_ylim([-0.05, 1.05])
            # Save the figure
            output_path = output_path + "/" + event + ".png"
            logger.info("Created " + output_path)
            fig.savefig(output_path, bbox_inches='tight')
        else:
            logger.info("Event not found.")

        logger.info("Done.")

        exit(0)

    except Exception as error:
        logger.error(repr(error))
        logger.error("Aborting execution")
        sys.exit(1)


if __name__ == '__main__':
    main()