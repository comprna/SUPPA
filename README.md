Here we provide an extensive explanation about how SUPPA works and all the running options.
Read our [short SUPPA tutorial on an example dataset](https://github.com/comprna/SUPPA/wiki).

----------------------------
# Table of Contents
----------------------------

   * [Overview](#overview)
   * [Set up](#set-up)
      * [Requirements](#requirements)
   * [Command and subcommand structure](#command-and-subcommand-structure)
   * [Event generation](#event-generation)
      * [Input files](#input-files)
      * [Command and options](#command-and-options)
      * [Output files](#output-files)
         * [ioe](#ioe)
         * [GTF](#gtf)
   * [PSI calculation for each event](#psi-calculation-for-each-event)
      * [Input files](#input-files-1)
      * [Command and options](#command-and-options-1)
      * [Output files](#output-files-1)
   * [PSI calculation for each isoform](#psi-calculation-for-each-isoform)
   * [Combining multiple expression files](#combining-multiple-expression-files)
   * [Differential splicing analysis](#differential-splicing-analysis)
      * [Input files](#input-files-2)
      * [Command and options](#command-and-options-2)
      * [Output files](#output-files-2)
         * [dpsi file](#dpsi-file)
         * [psivec file](#psivec-file)
   * [Cluster analysis](#cluster-analysis)
      * [Input files](#input-files-3)
      * [Command and options](#command-and-options-3)
      * [Output files](#output-files-3)
         * [clustvec](#clustvec)
   * [License](#license)

----------------------------
# Overview
----------------------------

SUPPA is a tool to study splicing across multiple conditions at high speed and accuracy. The tool has four modular operations that are run separately to:

- Generate events from an annotation
- Quantify event inclusion levels (PSIs) from multiple samples
- Calculate differential splicing across multiple conditions with replicates
- Cluster events according to PSI values across conditions

![Slide1.jpg](https://bitbucket.org/repo/4gEBMd/images/3120745382-Slide1.jpg)

 
**Fig 1.** SUPPA generates the different alternative splicing events from an input annotation file (GTF format). The method reads transcript and gene information solely from the "exon" lines in the GTF. It then generates the events and outputs an ioe file, which contains the for each event the transcripts that describe either form of the event. Specifically, it provides the transcripts that contribute to the numerator (one form of the event) and the denominator (both forms of the event) of the PSI calculation. For the generation of PSI values, SUPPA reads the ioe file generated in the previous step and a transcript expression file with the transcript abundances to calculate the PSI value for each of the events. 

![figure2_for_readme_bitbucket.jpg](https://bitbucket.org/repo/4gEBMd/images/612499447-figure2_for_readme_bitbucket.jpg)

**Fig 2.** (a) SUPPA calculates the magnitude of splicing change (ΔPSI) and their significance across multiple biological conditions, using two or more replicates per condition. Conditions are analyzed in a sequential order specified as input. (b) Statistical significance is calculated by comparing the observed ΔPSI between conditions with the distribution of the ΔPSI between replicates as a function of the gene expression (measured as the expression of the transcripts defining the events). Alternatively, when there is a large (>10) number of replicates per condition, SUPPA can apply a classical statistical test (Wilcoxon) per event (not shown in the figure). (c) Using the output from the differential splicing analysis, SUPPA can also cluster events with similar splicing patterns using a density-based clustering algorithm.
 
We provide below detailed information on how to install and run SUPPA. Please join the SUPPA google-group for sharing your thoughts and questions (suppa-users@googlegroups.com). 


----------------------------
# Set up
----------------------------

## Requirements

SUPPA has been developed in Python 3.4. 

If necessary, to install python3 we recommend to download from the official site ( https://www.python.org/downloads/ ) the corresponding version for your OS.

SUPPA uses the following modules:

- SciPy ( 0.15.1 )
- NumPy ( 1.11.0 )
- Pandas ( 0.18.0 )
- statsmodels ( 0.6.1 )
- scikit-learn ( 0.16.1 )

In case you do not have already these modules, we reccomend to use pip3 to install the specific version of the modules. To install pip3, please go to the official site ( https://pip.pypa.io/en/latest/installing/ ),
download the file "get-pip.py", and run the following command:
```
python3 get-pip.py
sudo pip3 install --upgrade pip3
```

Then, to install the modules:
```
sudo pip3 install scipy==0.15.1 numpy==1.11.0 pandas==0.18.0 statsmodels==0.6.1 scikit-learn==0.16.1
```

SUPPA is ready to use. Once downloaded, it can be used directly from the command line.

----------------------------
# Command and subcommand structure
----------------------------

SUPPA works with a command/subcommand structure:

```
suppa.py subcommand options
```
where the subcommand can be one of these five:

- **generateEvents**    : Generates events from an annotation.
- **psiPerEvent**       : Quantifies event inclusion levels (PSIs) from multiple samples.
- **psiPerIsoform**       : Quantifies isoform inclusion levels (PSIs) from multiple samples.
- **diffSplice**        : Calculate differential splicing across multiple conditions with replicates.
- **clusterEvents**     : Cluster events according to PSI values across conditions.

----------------------------
# Event generation
----------------------------

SUPPA generates the different AS events from an input annotation file (GTF format). The method reads transcript and gene information solely from the "exon" lines. It then generates the events and outputs an ioe file. The ioe file provides for each AS event in a gene, the transcripts that describe either form of the event. Specifically, it provides the transcripts that contribute to the numerator (one form of the event) and the denominator (both forms of the event) of the PSI calculation. 

Different event types generated by SUPPA:

- Skipping Exon (SE)
- Alternative 5'/3' Splice Sites (A5/A3)
- Mutually Exclusive Exons (MX)
- Retained Intron (RI)
- Alternative First/Last Exons (AF/AL)


![Events_legend1.jpg](https://cloud.githubusercontent.com/assets/23315833/22699555/19788560-ed58-11e6-8340-1ce83445bdb3.png)

![Events_legend2.jpg](https://cloud.githubusercontent.com/assets/23315833/22699557/199ac7c4-ed58-11e6-8512-2d3950001a8d.png)


**Fig 3.** The figure describes the nomenclature for events in the forward (upper panel) and reverse (lower panel) strands. Each event is identified uniquely by a set of coordinates: The start (s) and end (e) coordinates for the different exonic regions involved in the event. The external coordinates of the event are only used for the RI, AF and AL events. The form of the alternative splicing event that includes the region in black is the one for which the relative inclusion level (PSI) is given. The gray area denotes the other alternative form of the event. For instance, for RI the inclusion level is given for the form that retains the intron. **Important**: for the non-symmetrical events note that the meaning of the coordinates varies depending on the strand


![Slide2.jpg](https://bitbucket.org/repo/4gEBMd/images/1362804745-Slide2.jpg)

**Fig 4.** In case the option for variable boundaries is used **-b V** (see below), an user input variability (detault: 10nt) is allowed in some of the boundaries (indicated with grey arrows). For some even types there will be less coordinates defining the event, as some boundaries are not tested any more. The variable boundaries (grey arrows) allow incorporating other transcripts contributing to the event, and therefore mimicking more closely the PSI calculate from RT-PCR primers. The events are still reported with the indicated coordinates. Those events that are redundant in this new setting will have the same transcripts contributing to the inclusion and skipping forms, however, they will be still reported in separate lines. As before, the form of the alternative splicing event that includes the region in black is the one for which the relative inclusion level or PSI is given, whereas the gray area denotes the other alternative form of the event. 


## Input files

An annotation file in GTF format is required (see e.g. http://mblab.wustl.edu/GTF22.html):

```
chr14 Ensembl exon  73741918  73744001  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 
chr14 Ensembl exon  73749067  73749213  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1";  
chr14 Ensembl exon  73750789  73751082  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 
chr14 Ensembl exon  73753818  73754022  0.0 - . gene_id "ENSG00000000001"; transcript_id "ENST00000000001.1"; 

```

The *generateEvents* operation just uses the lines where the feature (column 3) is "exon". It then reads the different transcripts and genes. For that purpose "gene_id" and "transcript_id" tags are required in the attributes field (column 9). This command also generates a GTF file with the calculated events and with a track header ready to be uploaded into the UCSC browser for visualization (see below). 

## Command and options

To generate the events from the GTF file one has to run the following command:

```
suppa.py generateEvents [options]
```
List of options available:

- **-i**  | **--input-file**: a GTF format file containing at least "exon" lines

- **-o**  | **--output-file**: name of the output file without any extension

- **-e**  | **--event-type**: space separated list of events to generate from the following list:

  - **SE**: Skipping exon (SE)
  - **SS**: Alternative 5' (A5) or 3' (A3) splice sites (generates both)
  - **MX**: Mutually Exclusive (MX) exons
  - **RI**: Retained intron (RI)
  - **FL**: Alternative First (AF) and Last (AL) exons (generates both)

- **-b** | **--boundary**: [S,V]
        - Boundary type. Options: S -- Strict (Default) V -- Variable

- **-t** | **--threshold**: THRESHOLD
        - Variability treshold (Default: 10nt). In case of strict boundaries this argument is ignored

- **-l**  | **--exon-length**: defines the number of nucleotides to display in the output GTF. (Default: 100 nt)

- **-p**  | **--pool-genes**: pools together into a composite gene those transcripts overlapping in the same strand that share at least one splice site.

- **-h**  | **--help**: display the help message describing the different paramenters


The command line to generate the different output files will be of the form:

```
suppa.py generateEvents -i <input-file> -o <output-file> -e <list-of-events>
```

## Output files

The *generateEvents* operation outputs two files:

1. An *ioe* file that shows the relationship between each event and the trasnscripts that that particular event.

2. A *GTF* file (with a track header) to load into the [UCSC genome browser](http://genome-euro.ucsc.edu/cgi-bin/hgGateway?redirect=auto&source=genome.ucsc.edu) to visualize the different events.

The name of the output is generated as follows:

**&lt;output-file&gt;**.**&lt;event-type&gt;**.ioe/gtf
    
**&lt;output-file&gt;**: -o option specified when launching the program.

**&lt;event_type&gt;**: a two letter code referring to the event type (SE, A3, A5, MX, RI, AF, AL).


### ioe

*ioe* stands for Isoforms Overlapping Event and is a tab separated file with 5 different fields:

1. **seqname**: field 1 from the input GTF file of the generateEvents operation (generally the chromosome name)

2. **gene_id**: ID of the gene where the event appears. (comma separated list in case of using the --pool-genes option)

3. **event_id**: ID of the event.

4. **Inclusion transcripts**: IDs of the transcripts that define the form of the event for which we calculate the PSI (e.g. exon inclusion) and contribute to the numerator of the PSI formula. For more details see Figures above. 

5. **Total transcripts**: IDs of the all transcripts that define either of the two forms of the event (e.g. inclusion and skipping) and contribute to the denominator of the PSI formula. For more details see Figures above.


**Event ID:** The event_id is formatted as follows:
```
<gene_id>;<event-type>:<seqname>:<coordinates-of-the-event>:<strand>
```
where:

- &lt;gene_id&gt;: is the gene where the event takes place
- &lt;event-type&gt;: correspond to the two letter code of the event from the following list.
  - **SE**: Skipping Exon
  - **A5**: Alternative 5' Splice Site
  - **A3**: Alternative 3' Splice Site
  - **MX**: Mutually Exclusive Exon
  - **RI**: Retained Intron
  - **AF**: Alternative First Exon
  - **AL**: Alternative Last Exon  
- &lt;seqname&gt;: coordinate reference system (e.g. chr1) 
- &lt;coordinates-of-the-event&gt;: the coordinates of the event depends on the type of event (see above)
- &lt;strand&gt;: either '+' or '-'


**Inclusion transcripts**: The IDs must be the same as those provided in the expression (TPM) files. These transcripts define the inclusion form of the event. This form is chosen according to the convention described above in Figures 3 and 4: 

- **SE**: transcripts including the middle exon
- **A5/A3**: transcripts minimizing the intron length
- **MX**: transcripts containing the alternative exon with the smallest (left most) start coordinate.
- **RI**: transcripts that have the retain intron
- **AF/AL**: transcripts maximizing the intron length

**Total transcripts**: The IDs must be the same as those provided in the expression (TPM) files. Total transcripts are the transcripts that define both forms of the event and therefore contribute to the denominator of the PSI formula. 

An example of an *ioe* file is the following one:
```
seqname gene_id event_id  inclusion_transcripts total_transcripts
chr1  ENSG00000067606   ENSG00000067606;AF:chr1:1986880:1987001-1987923:1987118:1987315-1987923:+ ENST00000479566.1 ENST00000466651.1,ENST00000479566.1
chr1  ENSG00000067606   ENSG00000067606;AF:chr1:2013687:2013899-2066701:2043093:2043256-2066701:+ ENST00000481140.1 ENST00000470986.1,ENST00000481140.1
chr1  ENSG00000067606   ENSG00000067606;AF:chr1:2006294:2006572-2066701:2043093:2043256-2066701:+ ENST00000470596.1 ENST00000470986.1,ENST00000470596.1
```

### GTF

This output file is aimed for visualization and it is a regular *gtf* in half-open coordinate convention with a track header to be uploaded in UCSC. In this file, each of the two possible forms of an event is described as a separated *transcript_id*. Moreover, the *gene_id* is defined as the event_id.

The **-l** | **--exon-length** option of the *eventGenerator* operation defines the length used to display the constitutive exons in the GTF file and it is used only for visualization purposes. For instance, in a *Skipping exon* event the flanking exons will be shown with a length given by this parameter.


-------------------
# PSI calculation for each event
-------------------

SUPPA reads the ioe file generated in the previous step and a transcript expression file with the transcript abundances to calculate the PSI value per sample for each of event. 

## Input files

An **ioe** file and a "**transcript expression file**" are required as input. 

The transcript expression file is a tab separated file where each line provides the estimated abundance of each transcript (ideally in TPM units). This file might contain multiple columns with the expression values in different samples. The expression file must have a header with the naming of the different expression fields, i.e., the sample name of each expression value. 

An example of a transcript expression file for one single sample:
```
sample1
transcript1 <expression>
transcript1 <expression>
transcript1 <expression>
```
A transcript expression file with multiple samples:

```
sample1 sample2 sample3 sample4
transcript1 <expression>  <expression>  <expression>  <expression>
transcript2 <expression>  <expression>  <expression>  <expression>
transcript3 <expression>  <expression>  <expression>  <expression>
```

**Note:** these files have a header with only the sample names (1 less column)

## Command and options
To calculate the psi from the *ioe* and the *transcript expression file* one has to run the following command:

```
suppa.py psiPerEvent [options]
```
List of options available:

- **-e** | **--expression-file**: transcript expression file containing the abundances of all transcripts

- **-i** | **--ioe-file**: input file with the definition of the events and corresponding transcripts

- **-f** | **--total-filter**: Minimum total expression of the transcripts involved in the event. If used, it will filter out the events that do not reach this total expression value for the transcripts defining the event (the denominator of the PSI calculation).

- **-o** | **--output-file**: output file name

- **-h**  | **--help**: display the help message describing the different parameters

- **-m** | **--mode**: verbose mode to choose from DEBUG, INFO, WARNING, ERROR and CRITICAL

An example of the usage of the program is:

```
suppa.py psiPerEvent --ioe-file <ioe-file> --expression-file <expression-file> -o <output-file> 
```

## Output files

This operation generates a *psi* file per each expression file. This file follows the same exact format of the input expression file used in the *psiPerEvent* operation. The PSI values are normalized between 0 and 1 ([0,1]). There are two exceptions for the PSI value:

- **-1**: PSI = -1 if the expression is 0 for all the transcripts involved in the event or if the event does not pass any the transcript expression filter (see *Psi calculation, Command and options*).
- **NA**: PSI = NA if one or more transcripts of the event do not appear in the expression file.


-----------------------------------
# PSI calculation for each isoform
------------------------------------

SUPPA also allows to calculate the PSI per transcript isoform, which is simply defined as the normalized abundance of the transcript over the abundances of all transcripts in the same gene. In this case we only need the *transcript expression file* and the original annotation file in GTF as inputs. The command is as follows:
```
suppa.py psiPerIsoform [options]
```
List of options available:

- **-g** | **--gtf-file**: input file with the annotation

- **-e** | **--expression-file**: expression file containing the abundances of all transcripts (ideally in TPM units)

- **-o** | **--output-file**: output file

- **-h**  | **--help**: display the help message describing the different parameters

- **-m** | **--mode**: verbose mode to choose from DEBUG, INFO, WARNING, ERROR and CRITICAL

An example of the usage of the program is:

```
suppa.py psiPerIsoform -g <gtf-file> -e <expression-file> -o <output-file> 
```

The output file is similar as for events, but where each line has the PSI of each isoform in every given sample. 

---------------------------------------
# Combining multiple expression files
------------------------------------

Transcript expression files used with SUPPA typically come from multiple calculations with different samples. To facilitate the generation of a single file with all the transcript expression values for all samples, SUPPA distribution includes a program to combine multiple simple transcript expression files into one single file:

```
suppa.py joinFiles [options]
```


where the options are:

- **-i** | **--input-files**:    Space separated list of the files to be joined.

- **-f** | **--file-extension**  Extension of the output file. Required.

- **-o** | **--output-file**:   name of the output file.

- **-h** | **--help**            


We show below an example of the usage of the program for reading multiple output files from Sailfish to join together the 3rd column, given that all files have in the first column the transcript ids (which are kept for the output):

```
suppa.py joinFiles -f tpm -i sample1.tpm sample2.tpm sample3.tpm -o all_samples_tpms  
```

The output will look like an expression file with multiple files as described above.

**Important:** Right now this method assumes that the individual files have in the first column the IDs that are common across all files, and in the second column the PSI or TPM values to be joined. This will be changed in the future to a more flexible method so that generic tab separated files can be used directly. 

-------------------
# Differential splicing analysis
-------------------

SUPPA reads the PSI and transcript expression values from multiple samples, grouped by condition, and the ioe file, to calucate the events that are differentially spliced between a pair of conditions. SUPPA can test multiple conditions with a variabe number of samples per condition. 

## Input files

The **ioe** file, a **PSI** files and the **transcript expression** files are required.

The PSI and transcript expression files are both tab separated files, where each line provides, respectively, the PSI values and the estimated abundance of each transcript. These files have to contain multiple columns with the PSI values and expression values in different samples, keeping the same ordering of samples in both files.

The PSI and transcript expression files must have a header with the sample names. For instance, the PSI file could look like this: 
```
sample1 sample2 sample3 sample4
event1  <psi_value> <psi_value> <psi_value> <psi_value>
event2  <psi_value> <psi_value> <psi_value> <psi_value>
event3  <psi_value> <psi_value> <psi_value> <psi_value>
```

and the transcript expression file could be of the form: 

```
sample1 sample2 sample3 sample4
transcript1 <expression>  <expression>  <expression>  <expression>
transcript2 <expression>  <expression>  <expression>  <expression>
transcript3 <expression>  <expression>  <expression>  <expression>
```
where the expression values are ideally given in TPM units.

**Note:** these files have a header with only the sample names (1 less column)

**Important:** SUPPA will read one PSI file and one transcript expression file per condition. Each of these files will contain the multiple replicates (or individual samples) that are grouped into a given condition, and the should be specified in the same order within the files.

## Command and options
To calculate the dpsi from the *ioe*, *psi* and the *expression file* one has to run the following command:
```
suppa.py diffSplice [options]
```
List of options available:

- **-m** | **--method**: The method to use to calculate the significance (empirical/classical)

- **-p** | **--psi**:  PSI files. The psi and expression files must be given in the same order.

- **-e** | **--expression-files**: Transcript expression files. The expression files and PSI files must be given in the same order.

- **-i** | **--ioe**: File with the definition of the events

- **-a** | **--area**: Integer indicating the number of points in the local area of the delta PSI - average TPM distribution. (Default: 1000).

- **-l** | **--lower-bound**:  Lower-bound for the absolute delta PSI value to test for significance. Events with less than this delta PSI will not be tested. (Default: 0).

- **-pa** | **--paired**: Boolean. Indicates if replicates across conditions are paired. (Default: False).

- **-gc** | **--gene-correction**: Boolean. If True, SUPPA correct the p-values by gene. (Default: False).

- **-al** | **--alpha**: Family-wise error rate to use for the multiple test correction. (Default: 0.05).

- **-o** | **--output**: Name of the output

- **-h**  | **--help**: display the help message describing the different parameters

An example of the usage of the program is:

```
suppa.py diffSplice --method <empirical> --ioe <ioe-file> --psi <Cond1.psi> <Cond2.psi> --expression-file <Cond1_expression-file> <Cond2_expression-file> --area <1000> --lower-bound <0.05> -gc -o <output-file>
```

**Note:** ΔPSI values are calculated pairwise between each pair of adjacent conditions as provided, and calculating the PSI different between a given condition and the previous one. For instance, for three conditions 1,2,3, the ΔPSI values will be calculated for 3 - 2 and 2 - 1. 

## Output files

The differential splicing operation generates a *dpsi* file and a *psivec* file. 


### dpsi file

*dpsi* is a tab separated file, with the event ID in the first column, followed by a variable even number of fields, two for each pair of conditions being compared:


1. **Cond1_Cond2_dPSI**: Event PSI difference (ΔPSI) between Cond2 and Cond1.

2. **Cond1_Cond2_pvalue**: Significance of the difference of PSI between Cond2 and Cond1

An example of an *dpsi* file is the following one:
```
Cond1_Cond2_dPSI  Cond1_Cond2_p-val
ENSG00000000003;A5:chrX:99890743-99891188:99890743-99891605:- nan 1.0
ENSG00000000419;A3:chr20:49557492-49558568:49557470-49558568:-  0.1307455855  0.0484515485
ENSG00000000419;SE:chr20:49557470-49557642:49557746-49558568:-  0.0469449126  0.0954045953
ENSG00000000419;SE:chr20:49557470-49557666:49557746-49558568:-  0.0164051352  0.1518481518
```

### psivec file

*psivec* is a tab separated file, with the event ID in the first column, followed by a variable of fields, each giving the PSI value per replicate:


1. **sample PSI**: PSI value for a given replicate/condition.

An example of an *psivec* file is the following one:
```
sample1 sample2 sample1 sample2
ENSG00000000003;A5:chrX:99890743-99891188:99890743-99891605:- 0.14343855447180462 0.02929736320730957 0.12495621749266525 -1.0
ENSG00000000003;A5:chrX:99890743-99891605:99890743-99891790:- 1.0     1.0     1.0       -1.0
ENSG00000000419;A3:chr20:49557492-49557642:49557470-49557642:-  0.0     0.0     0.0       0.0
ENSG00000000419;A3:chr20:49557492-49558568:49557470-49558568:-  0.004508763698090278  0.0     0.028759955911389606 0.0638008250417565 
ENSG00000000419;A5:chr20:49557470-49557642:49557470-49557666:-  0.4495047957648805  0.490664195110865 0.6304018235083193   0.6781317189253738 
```

**Note:** these files have a header with only the sample names (1 less column)

-------------------
# Cluster analysis
-------------------

Using the PSI values in all samples, and the information of which events change significantly in at least one comparison, SUPPA calculates clusters of events using a density-based clustering. Density-based clustering has the advantage that it can put together into the same cluster two events that even though they might not have similar PSI values, their behaviour is similar across conditions and there are sufficient events between them with a similar behaviour. 

## Input files

A **dpsi** file and **psivec** file are required (see definitions above). The dpsi file contains the information about which events are signficantly differentially spliced in each pairwise comparison. For instance, a dpsi file with data from the comparison of three conditions will look like:
```
Cond1_Cond2_dPSI  Cond1_Cond2_p-val Cond2_Cond3_dPSI  Cond2_Cond3_p-val
event1    <dpsi_value>    <dpsi_value>    <dpsi_value>    <dpsi_value>
event2    <dpsi_value>    <dpsi_value>    <dpsi_value>    <dpsi_value>
event3    <dpsi_value>    <dpsi_value>    <dpsi_value>    <dpsi_value>
```

The psivec file contains the PSI values for all samples, either per replicate or the average PSI value per condition, averaging over the replicates. For instance:

```
sample1 sample2 sample3 sample1 sample2 sample3 
event1    <psi_value> <psi_value> <psi_value> <psi_value> <psi_value> <psi_value>
event2    <psi_value> <psi_value> <psi_value> <psi_value> <psi_value> <psi_value>
event3    <psi_value> <psi_value> <psi_value> <psi_value> <psi_value> <psi_value> 
```


## Command and options

SUPPA will use the psivec file to cluster events according to the PSI values across samples using those events that show significant change in at least one pairwise comparison using the dpsi file. Two methods are available: DBSCAN and OPTICS. Both methods require as input the minimum number of events in a cluster. OPTICS also requires as the maximum reachability distance (s), which represents the maximum distance in PSI space of an event to a cluster  To perform the clustering from the *dpsi* and the *psivec* one has to run the following command:

```
suppa.py clusterEvents [options]
```
List of options available:

- **-d** | **--dpsi**: Input .dpsi file generated from SUPPA diffSplice operation.

- **-p** | **--psivec**: Input .psivec file generated from SUPPA diffSplice operation.

- **-st** | **--sig-threshold**: p-value threshold to consider an event significant from the dpsi file.

- **-dt** | **--dpsi-threshold**: Lower-bound for the absolute delta PSI value to cluster. (Default: 0.05)

- **-e** | **--eps**: Maximum (Euclidean) distance (between 0 and 1) to consider two events as members of the same cluster. (Default: 0.05).

- **-s** | **--separation**: maximum distance in PSI space of an event to a cluster. Required for OPTICS method

- **-n** | **--min-pts**: Minimum number of events required per cluster. (Default: 20).

- **-g** | **--groups**:  Ranges of column numbers specifying the replicates per condition.
                        Column numbers have to be continuous, with no
                        overlapping or missing columns between them. Ex:1-3,4-6

- **-c** | **--clustering**: Clustering method to use (DBSCAN, OPTICS). (Default: DBSCAN)

- **-o** | **--output**: Name of the output file

- **-h**  | **--help**: display the help message describing the different parameters

An example of the usage of the program is:

```
suppa.py clusterEvents --dpsi <dpsi-file> --psivec <psivec-file> --sig-threshold <0.05> --eps <0.05> --min-pts <20> --groups <1-3,4-6> -o <output-file> 

```

**Note:** In the example above, the replicates 1,2,3 and 4,5,6 are grouped into 2 conditions. 

## Output files

SUPPA clusterEvents operation generates a *clustvec* file per each pair of dpsi and psivec files. 


### clustvec

*clustvec* is a tab separated file, with the event ID in the first column, the cluster ID as second column, and then followed by a variable even number of fields, one for each condition, giving the mean PSI value per condition:

1. **event_id**: Event ID.

2. **cluster_id**: Cluster ID. An integer starting from 0 that indicates to which cluster the events belongs to. Singletons (non-clustered events) are indicated with -1.

3. **cond_PSI_avg**: (variable number of columns): Average PSI value for each condition.

An example of a *clustvec* file is the following one:

**Note:** The clustvec file doesn't contain a header. The header is for illustraion purpouse only. 
```
event_id  cluster_id  cond1_PSI_avg cond2_PSI_avg
ENSG00000004897;SE:chr17:45247408-45249283:45249430-45258928:-  1 0.967582617827  0.885135296685
ENSG00000004961;A5:chrX:11129580-11130139:11129543-11130139:+ 1 0.950872555163  0.720740016826
ENSG00000005189;A3:chr16:20824624-20826249:20824624-20826307:+  0 0.632147193342  0.332088364504
ENSG00000005189;A5:chr16:20818139-20818274:20818027-20818274:+  0 0.344510537149  0.161862745525
ENSG00000005189;SE:chr16:20851790-20855256:20855348-20855951:+  -1  0.853283401943  0.548674896371
ENSG00000005194;SE:chr16:57463192-57464170:57464251-57466399:-  -1  0.971516032983  0.829712829609
```

----------------------------
# License
----------------------------

SUPPA is released under the MIT license.
