# SUPPA2

Fast, accurate, and uncertainty-aware differential splicing analysis across multiple conditions. Publications:

* Trincado JL, Entizne JC, Hysenaj G, Singh B, Skalic M, Elliott DJ, Eyras E. [SUPPA2: fast, accurate, and uncertainty-aware differential splicing analysis across multiple conditions](https://www.ncbi.nlm.nih.gov/pubmed/29571299). Genome Biol. 2018 Mar 23;19(1):40.
* Alamancos GP, Pagès A, Trincado JL, Bellora N, Eyras E. [Leveraging transcript quantification for fast computation of alternative splicing profiles](https://www.ncbi.nlm.nih.gov/pubmed/26179515). RNA. 2015 Sep;21(9):1521-31.

Here we provide an extensive explanation about how SUPPA works and all the running options.
Read our [short SUPPA tutorial on an example dataset](https://github.com/comprna/SUPPA/wiki).

----------------------------
# Table of Contents
----------------------------

   * [Overview](#overview)
   * [Installation](#installation)
   * [Command and subcommand structure](#command-and-subcommand-structure)
   * [Generation of transcript events and local alternative splicing events](#generation-of-transcript-events-and-local-alternative-splicing-events)
      * [Input files](#input-files)
      * [Command and options](#command-and-options)
      * [Output files](#output-files)
         * [ioi](#ioi)
         * [ioe](#ioe)
         * [gtf](#gtf-for-local-events)
   * [PSI calculation for Transcripts and Events](#psi-calculation-for-transcripts-and-events)
      * [Input files](#input-files-1)
         * [PSI per transcript Isoform](#psi-per-transcript-isoform)
         * [PSI per local event](#psi-per-local-event)
      * [Output files](#output-files-1)
   * [Combining multiple expression files](#combining-multiple-expression-files)
   * [Differential splicing analysis for transcripts and local events](#differential-splicing-analysis-for-transcripts-and-local-events)
      * [Input files](#input-files-2)
      * [Command and options](#command-and-options-2)
      * [Differential transcript usage](#differential-transcript-usage)
      * [Differential splicing with local events](#differential-splicing-with-local-events)
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

SUPPA is a flexible and powerful tool to study splicing at the transcript isoform or at the local alternative splicing event level, across multiple conditions, and at high speed and accuracy. SUPPA has various modular operations that can be run separately to:

- Generate transcript events and local alternative splicing (AS) events from an annotation
- Quantify transcript and local AS event inclusion levels (PSIs) from multiple samples
- Calculate differential splicing for AS events and differential transcript usage across multiple conditions with replicates
- Cluster events and transcripts according to PSI values across conditions

![Slide1.jpg](https://bitbucket.org/repo/4gEBMd/images/3120745382-Slide1.jpg)

 
**Fig 1.** SUPPA generates the alternative splicing events from an input annotation file (GTF format). The method reads transcript and gene information solely from the "exon" lines in the GTF. It then generates the events and outputs an ioe file, which contains the for each event the transcripts that describe either form of the event. Specifically, it provides the transcripts that contribute to the numerator (one form of the event) and the denominator (both forms of the event) of the PSI calculation. For the generation of PSI values, SUPPA reads the ioe file generated in the previous step and a transcript expression file with the transcript abundances to calculate the PSI value for each of the events. 

![figure2_for_readme_bitbucket.jpg](https://bitbucket.org/repo/4gEBMd/images/612499447-figure2_for_readme_bitbucket.jpg)

**Fig 2.** (a) SUPPA calculates the magnitude of splicing change (ΔPSI) (for events or transcripts) and its significance across multiple biological conditions, using two or more replicates per condition. Conditions are analyzed in a sequential order specified as input. (b) Statistical significance is calculated by comparing the observed ΔPSI between conditions with the distribution of the ΔPSI between replicates as a function of the expression of the transcripts defining the events (for events) or as a functon of the gene expression (for transcripts). When there is a large (>10) number of replicates per condition, you can also run SUPPA with a classical statistical test (Wilcoxon) per local event or per transcsript (not shown in the figure). (c) Using the output from the differential splicing analysis, SUPPA can cluster events or transcripts with similar splicing patterns across conditions using a density-based clustering algorithm.
 
We provide below detailed information on how to install and run SUPPA. Please join the SUPPA google-group for sharing your thoughts and questions (suppa-users@googlegroups.com). 


----------------------------
# Installation
----------------------------

SUPPA has been developed in Python 3.4. 

If necessary, to install python3 we recommend to download from the official site https://www.python.org/downloads/ the corresponding version for your OS.

A installation using pip is available using the next command:

```
pip install SUPPA==2.3 
```
By default SUPPA is installed into the Python package library directory. The following command can be executed to obtain the directory location:

```
pip show SUPPA 
```

SUPPA is ready to use. Once downloaded, it can be used directly from the command line by specifying the absolute path to the SUPPA executable (suppa.py).

Another option is via bioconda (thanks to Devon Ryan)

```
conda install -c bioconda suppa
```

Alternatively, you can obtain the code for the SUPPA2.3 release here:
https://github.com/comprna/SUPPA/releases/tag/v2.3

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


Note: Unless suppa was installed using a dependency manager (i.e. suppa would just be in the PATH), you can run suppa as:
```
python3.4 suppa.py subcommand options
```

----------------------------
**Generation of transcript events and local alternative splicing events**
==============

----------------------------

SUPPA can work with local alternative splicing events or with transcripts "events" per gene. The local alternative splicing events are
standard local splicing variations (see below), whereas a transcript event is an isoform-centric approach, where each isoform in a gene is described separately. 

SUPPA generates the AS events or transcript events from an input annotation file (GTF format). The method reads transcript and gene information solely from the "exon" lines. It then generates the events and outputs an event file: **.ioe** format for local AS events, and **.ioi** format for transcripts. 

The ioe file provides for each AS event in a gene, the transcripts that describe either form of the event. Specifically, it provides the transcripts that contribute to the numerator (one form of the event) and the denominator (both forms of the event) of the PSI calculation. 

The ioi file provides for each transcript in a gene, the set of all transcripts from that gene from which the transcript relative abundance is calculated. Transcript events are important as they can describe splicing variations that are complex and cannot be encapsulated in a simple event (see e.g. Sebestyen et al. 2015 https://www.ncbi.nlm.nih.gov/pubmed/25578962).

Different local event types generated by SUPPA:

- Skipping Exon (SE)
- Alternative 5'/3' Splice Sites (A5/A3) (generated together with the option SS)
- Mutually Exclusive Exons (MX)
- Retained Intron (RI)
- Alternative First/Last Exons (AF/AL) (generated together with the option FL)


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

The *generateEvents* operation uses the lines where the feature (column 3) is "exon". It then reads the different transcripts and genes. For that purpose "gene_id" and "transcript_id" tags are required in the attributes field (column 9). For local AS events, this command also generates a GTF file with the calculated events and with a track header ready to be uploaded into the UCSC browser for visualization (see below). 

## Command and options

To generate the events from the GTF file one has to run the following command:

```
python3.4 suppa.py generateEvents [options]
```
List of options available:

- **-i**  | **--input-file**: a GTF format file containing at least "exon" lines

- **-o**  | **--output-file**: name of the output file without any extension

- **-f**  | **--format**: [ioe,ioi]
  	  - Required. Format of the event annotation file: ioe for local events, ioi for transcript events.

- **-p**  | **--pool-genes**: 
  	  - Optional. Redefine genes by clustering together transcripts by genomic stranded overlap and sharing at least one exon.
            It is crucial when creating ioe/ioi from annotations that are not loci-based, e.g.: RefSeq and UCSC genes.

- **-e**  | **--event-type**: (only used for local AS events) space separated list of events to generate from the following list:

  - **SE**: Skipping exon (SE) events
  - **SS**: Alternative 5' (A5) and 3' (A3) splice sites (it generates both)
  - **MX**: Mutually Exclusive (MX) exons
  - **RI**: Retained intron (RI)
  - **FL**: Alternative first (AF) and last (AL) exons (it generates both)

- **-b** | **--boundary**: [S,V]
        - Boundary type (only used for local AS events). Options: S -- Strict (Default) V -- Variable

- **-t** | **--threshold**: THRESHOLD
        - Variability treshold (Default: 10nt. Only used for local AS events). In case of strict boundaries this argument is ignored

- **-l**  | **--exon-length**: (only used for local AS events). Defines the number of nucleotides to display in the output GTF. (Default: 100 nt)

- **-h**  | **--help**: display the help message describing the different paramenters


The command line to generate local AS events will be of the form:

```
python3.4 suppa.py generateEvents -i <input-file.gtf> -o <output-file> -f ioe -e <list-of-events>
```

The command to generate the transcript "events" would be of the form:

```
python3.4 suppa.py generateEvents -i <input-file.gtf> -o <output-file> -f ioi 
```

## Output files

For transcripts, the *generateEvents* operation outputs one single file: An *ioi* file that shows the transcript "events" in each gene. This is a tab separated file with the following fields

### **ioi**
-------

1. **seqname**: field 1 from the input GTF file of the generateEvents operation (generally the chromosome name)

2. **gene_id**: ID of the gene where the event appears taken from the GTF file. (comma separated list in case of using the --pool-genes option)

3. **event_id**: ID of the event, formatted as **gene_id**;**transcript_id**

4. **transcript_id**: ID of the transcript that defines the event, for which the relative inclusion (PSI) is calculated.

5. **Total transcripts**: IDs of the all transcripts in the gene (including the transcript in 4.) 

Note that we call a transcript "event", a transcript in a gene. 


An example of an *ioi* file is the following one:
```
seqname	gene_id	event_id	inclusion_transcripts	total_transcripts
chr14   ENSG00000133961.15      ENSG00000133961.15;ENST00000556772.1    ENST00000556772.1       ENST00000554546.1,ENST00000356296.4,ENST00000557597.1,ENST00000555238.1,ENST00000556772.1,ENST00000355058.3,ENST00000359560.3,ENST00000555394.1,ENST00000454166.4,ENST00000544991.3,ENST00000560335.1,ENST00000559312.1,ENST00000554521.2,ENST00000555738.2,ENST00000535282.1,ENST00000554014.2,ENST00000553997.1,ENST00000557486.1,ENST00000555859.1,ENST00000555307.1,ENST00000554394.1,ENST00000554315.1,ENST00000556989.1,ENST00000555987.1,ENST00000556112.1,ENST00000557774.1,ENST00000554818.1,ENST00000557031.1,ENST00000556700.1,ENST00000556600.1,ENST00000553415.1,ENST00000557577.1,ENST00000557581.1
chr14   ENSG00000133961.15      ENSG00000133961.15;ENST00000554818.1    ENST00000554818.1       ENST00000554546.1,ENST00000356296.4,ENST00000557597.1,ENST00000555238.1,ENST00000556772.1,ENST00000355058.3,ENST00000359560.3,ENST00000555394.1,ENST00000454166.4,ENST00000544991.3,ENST00000560335.1,ENST00000559312.1,ENST00000554521.2,ENST00000555738.2,ENST00000535282.1,ENST00000554014.2,ENST00000553997.1,ENST00000557486.1,ENST00000555859.1,ENST00000555307.1,ENST00000554394.1,ENST00000554315.1,ENST00000556989.1,ENST00000555987.1,ENST00000556112.1,ENST00000557774.1,ENST00000554818.1,ENST00000557031.1,ENST00000556700.1,ENST00000556600.1,ENST00000553415.1,ENST00000557577.1,ENST00000557581.1
chr14   ENSG00000133961.15      ENSG00000133961.15;ENST00000556600.1    ENST00000556600.1       ENST00000554546.1,ENST00000356296.4,ENST00000557597.1,ENST00000555238.1,ENST00000556772.1,ENST00000355058.3,ENST00000359560.3,ENST00000555394.1,ENST00000454166.4,ENST00000544991.3,ENST00000560335.1,ENST00000559312.1,ENST00000554521.2,ENST00000555738.2,ENST00000535282.1,ENST00000554014.2,ENST00000553997.1,ENST00000557486.1,ENST00000555859.1,ENST00000555307.1,ENST00000554394.1,ENST00000554315.1,ENST00000556989.1,ENST00000555987.1,ENST00000556112.1,ENST00000557774.1,ENST00000554818.1,ENST00000557031.1,ENST00000556700.1,ENST00000556600.1,ENST00000553415.1,ENST00000557577.1,ENST00000557581.1
chr14   ENSG00000133961.15      ENSG00000133961.15;ENST00000557774.1    ENST00000557774.1       ENST00000554546.1,ENST00000356296.4,ENST00000557597.1,ENST00000555238.1,ENST00000556772.1,ENST00000355058.3,ENST00000359560.3,ENST00000555394.1,ENST00000454166.4,ENST00000544991.3,ENST00000560335.1,ENST00000559312.1,ENST00000554521.2,ENST00000555738.2,ENST00000535282.1,ENST00000554014.2,ENST00000553997.1,ENST00000557486.1,ENST00000555859.1,ENST00000555307.1,ENST00000554394.1,ENST00000554315.1,ENST00000556989.1,ENST00000555987.1,ENST00000556112.1,ENST00000557774.1,ENST00000554818.1,ENST00000557031.1,ENST00000556700.1,ENST00000556600.1,ENST00000553415.1,ENST00000557577.1,ENST00000557581.1
```



For local AS events, the *generateEvents* operation outputs two files:

1. An *ioe* file that shows the relationship between each event and the transcripts that define that particular event.

2. A *GTF* file (with a track header) to load into the [UCSC genome browser](http://genome-euro.ucsc.edu/cgi-bin/hgGateway?redirect=auto&source=genome.ucsc.edu) to visualize the different local AS events.

The name of the output is generated as follows:

**&lt;output-file&gt;**.**&lt;event-type&gt;**.ioe/gtf
    
**&lt;output-file&gt;**: -o option specified when launching the program.

**&lt;event_type&gt;**: a two letter code referring to the event type (SE, A3, A5, MX, RI, AF, AL).


The ioe file has the following fields:

### **ioe**
-------


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
seqname	gene_id	event_id	inclusion_transcripts	total_transcripts
chr14	ENSG00000133961.15	ENSG00000133961.15;SE:chr14:73763993-73789838:73789937-73822334:-	ENST00000556112.1    ENST00000556112.1,ENST00000555859.1
chr14	ENSG00000133961.15	ENSG00000133961.15;SE:chr14:73876776-73925740:73926011-73929235:-	ENST00000553415.1    ENST00000553415.1,ENST00000556700.1
chr14	ENSG00000133961.15	ENSG00000133961.15;SE:chr14:73744001-73745989:73746132-73749067:-	ENST00000557597.1,ENST00000555238.1,ENST00000556772.1,ENST00000359560.3,ENST00000355058.3,ENST00000535282.1   ENST00000557597.1,ENST00000555238.1,ENST00000556772.1,ENST00000359560.3,ENST00000355058.3,ENST00000535282.1,ENST00000554546.1,ENST00000356296.4,ENST00000555394.1,ENST00000454166.4,ENST00000560335.1,ENST00000555738.2,ENST00000554014.2
chr14	ENSG00000133961.15	ENSG00000133961.15;SE:chr14:73749213-73750789:73751082-73753818:-	ENST00000554546.1,ENST00000356296.4,ENST00000557597.1,ENST00000555238.1,ENST00000556772.1,ENST00000359560.3,ENST00000355058.3,ENST00000555394.1,ENST00000535282.1,ENST00000553997.1,ENST00000557486.1    ENST00000554546.1,ENST00000356296.4,ENST00000557597.1,ENST00000555238.1,ENST00000556772.1,ENST00000359560.3,ENST00000355058.3,ENST00000555394.1,ENST00000535282.1,ENST00000553997.1,ENST00000557486.1,ENST00000454166.4,ENST00000560335.1,ENST00000555738.2

```


**Important: --pool-genes**
-------


This option is important when creating ioe/ioi from annotations that are not loci-based, e.g.: RefSeq and UCSC genes.
Unlike Ensembl or Gencode, which annotate gene loci, i.e. a set of transcripts will be uniquely be related to a gene at a locus, other annotations, like UCSC and Refseq 
dowloaded from UCSC, do not have this unequivocal link of transcripts to a genomic locus. 

We thus re-cluster transcripts according to their overlap in genomic extent on the same strand, and according to them sharing 
at least one splice-site position in the genome. This method borrows from the (pre-Gencode) Ensembl function to create
genes from transcripts (see Curwen et al. 2004 https://www.ncbi.nlm.nih.gov/pubmed/15123590). If transcripts are not appropriately clustered, 
transcript relative abundances or event inclusion levels may not make sense.

Using the **--pool-genes** option is also advisable to use with Ensembl and Gencode. The annotation contains genes with overlapping transcripts
that share a great deal of sequence, hence their relative contribution to alternative splicing events should be taken into account. 


**GTF for local events**
-------

This output file is aimed for visualization and it is a regular *gtf* in half-open coordinate convention with a track header to be directly uploaded in UCSC. In this file, each of the two possible forms of an event is described as a separated *transcript_id*. Moreover, the *gene_id* is defined as the event_id. Note that SUPPA generally operates with closed coordinates.

The **-l** | **--exon-length** option of the *eventGenerator* operation defines the length used to display the constitutive exons in the GTF file and it is used only for visualization purposes. For instance, in a *Skipping exon* event the flanking exons will be shown with a length given by this parameter.


-------------------
**PSI calculation for Transcripts and Events**
==============
-------------------


SUPPA reads the ioi or ioe file generated in the previous step and a transcript expression file with the transcript abundances (TPM units) to calculate the relative abundance (PSI) value per sample for each transcript or each local event. 

## Input files

To calculate the PSI for local events, the required input files are an **ioe** file and a **transcript expression file**.

However, to calculate the PSI per isoform, the required input files are an **GTF** file and a **transcript expression file**.

The transcript expression file is a tab separated file where each line provides the estimated abundance of each transcript (in TPM units). This file might contain multiple columns with the expression values in different samples. The expression file must have a header with the naming of the different expression fields, i.e., the sample name of each expression value. 

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

### **PSI per transcript Isoform** ###

At the moment the PSI per transcript isoform is calculated in the following way:

```
python3.4 suppa.py psiPerIsoform [options]
```
List of options available:

- **-g** | **--gtf-file**: input file with the annotation

- **-e** | **--expression-file**: expression file containing the abundances of all transcripts (ideally in TPM units)

- **-o** | **--output-file**: output file

- **-h**  | **--help**: display the help message describing the different parameters

- **-m** | **--mode**: verbose mode to choose from DEBUG, INFO, WARNING, ERROR and CRITICAL

An example of the usage of the program is:

```
python3.4 suppa.py psiPerIsoform -g <gtf-file> -e <expression-file> -o <output-file>
```

### **PSI per local event** ###

To calculate the PSI value for each event from the *ioe* and the *transcript expression file* one has to run the following command:

```
python3.4 suppa.py psiPerEvent [options]

```
List of options available:

- **-e** | **--expression-file**: transcript expression file containing the abundances of all transcripts

- **-i** | **--ioe-file**: input file with the definition of the events and corresponding transcripts

- **-f** | **--total-filter**: Minimum total expression of the transcripts involved in the event (Default is zero). If used, it will filter out the events that do not reach this total expression value for the transcripts defining the event (the denominator of the PSI calculation).

- **-o** | **--output-file**: output file name

- **-h**  | **--help**: display the help message describing the different parameters

- **-m** | **--mode**: verbose mode to choose from DEBUG, INFO, WARNING, ERROR and CRITICAL

An example of the usage of the program is:

```
python3.4 suppa.py psiPerEvent --ioe-file <ioe-file> --expression-file <expression-file> -o <output-file> 
```

### **Output files** ###

These operations generate a *psi* file per each expression file. This file follows a similar format to the input expression file used in the *psiPerEvent* operation, where in the first column you find the event ID (transcript event or local AS event ID), and the following columns you find the PSI values in each sample. The PSI values are normalized between 0 and 1 ([0,1]). There are two exceptions for the PSI value:

- **NA**: PSI = if the expression is 0 for all the transcripts involved in the event or if the event does not pass any the transcript expression filter (see *PSI calculation, Command and options*). Also if one or more transcripts of the event do not appear in the expression file, an NA will be returned


Example of a PSI file for transcript "events":
```
sample1 sample2 sample3
ENSG00000000003.10;ENST00000373020.4	0.9891352763622773      0.8702144635757746      0.9758296160419822
ENSG00000000003.10;ENST00000494424.1	0.012124454316514821    0.01863029006342797     0.014658155865899102
ENSG00000000003.10;ENST00000496771.1	0.9067534025774345      0.9116862902516033      0.911059790520273
```

Example of a PSI file for local AS events:
```
sample1 sample2 sample3
ENSG00000000003.14;SE:chrX:100630866-100632485:100632568-100633405:-    0.9420871643913276      0.9304369364607596	0.8999960720735679
ENSG00000000419.12;SE:chr20:50940933-50941105:50941209-50942031:-       0.02302214526941425     0.08920651168486651	0.01093999389147593
```



---------------------------------------
# Combining multiple expression files
------------------------------------

Transcript expression files used with SUPPA typically come from calculations with multiple samples. To facilitate the generation of a single file with all the transcript expression values for all samples, SUPPA distribution includes a program to combine multiple simple transcript expression files into one single file:

```
python3.4 suppa.py joinFiles [options]
```


where the options are:

- **-i** | **--input-files**:    Space separated list of the files to be joined.

- **-f** | **--file-extension**  Extension of the output file. Required.

- **-o** | **--output-file**:   name of the output file.

- **-h** | **--help**            


We show below an example of the usage of the program for reading multiple output files from Sailfish to join together the 3rd column, given that all files have in the first column the transcript ids (which are kept for the output):

```
python3.4 suppa.py joinFiles -f tpm -i sample1.tpm sample2.tpm sample3.tpm -o all_samples_tpms  
```

The output will look like an expression file with multiple files as described above.

**Important:** Right now this method assumes that the individual files have in the first column the IDs that are common across all files, and in the second column the PSI or TPM values to be joined. This will be changed in the future to a more flexible method so that generic tab separated files can be used directly. 


-------------------
**Differential splicing analysis for transcripts and local events**
==============
-------------------

SUPPA reads the PSI for the (transcripts or local) events and the transcript expression values from multiple samples, grouped by condition, and the ioe/ioi file, to calucate the events that are differentially spliced between a pair of conditions. SUPPA can test multiple conditions with a variabe number of samples per condition. 

## Input files

The **ioe/ioi** file, the **PSI** files and the **transcript expression** files are required.

The PSI and transcript expression files are both tab separated files, where each line provides the (transcript or local) event PSI values and the estimated abundance of each transcript, respectively. These files must contain multiple columns with the PSI values and expression values in the different samples, always keeping the same order of samples in both files.

The PSI and transcript expression files must have a header with the sample names. For instance, the PSI file could look like this: 
```
sample1 sample2 sample3 sample4
event1  <psi_value> <psi_value> <psi_value> <psi_value>
event2  <psi_value> <psi_value> <psi_value> <psi_value>
event3  <psi_value> <psi_value> <psi_value> <psi_value>
```

and the transcript expression file should be then of the form:

```
sample1 sample2 sample3 sample4
transcript1 <expression>  <expression>  <expression>  <expression>
transcript2 <expression>  <expression>  <expression>  <expression>
transcript3 <expression>  <expression>  <expression>  <expression>
```

where the expression values are given in TPM units.

**Note:** these files have a header with only the sample names (1 less column)

**Important:** SUPPA will read one PSI file and one transcript expression file per condition. Each of these files will contain the multiple replicates (or individual samples) that are grouped into a given condition, and the should be specified in the same order within the files.

### **Command and options** ###
To calculate the dpsi from the *ioe*, *psi* and the *expression file* one has to run the following command:
```
python3.4 suppa.py diffSplice [options]
```

List of options available:

- **-m** | **--method**: The method to use to calculate the significance (empirical/classical)

- **-p** | **--psi**:  PSI files. The psi and expression files per condition must be given in the same order.

- **-e** | **--tpm**: Transcript expression (TPM) files. The expression files and PSI files per condition must be in the same order.

- **-i** | **--input**: Input file with the local or transcript events, .ioe or .ioi format, respectively.

- **-a** | **--area**: Integer indicating the number of points in the local area of the delta PSI - average TPM distribution. (Default: 1000).

- **-l** | **--lower-bound**:  Lower-bound for the absolute delta PSI value to test for significance. Events with less than this delta PSI will not be tested. (Default: 0).

- **-pa** | **--paired**: Indicates if replicates across conditions are paired.

- **-gc** | **--gene-correction**: Correction of the p-values by gene.

- **-al** | **--alpha**: Family-wise error rate to use for the multiple test correction. (Default: 0.05).

-  **-s** | **--save_tpm_events**: The average log TPM of the events will be saved in an external file.

-  **-c** | **--combination**: SUPPA will perform the analysis between all the possible combinations of conditions.

-  **-me** | **--median**: SUPPA will use the median to calculate the Delta PSI, instead of the mean.

-  **-th** | **--tpm-threshold**: Minimum expression (calculated as average TPM value within-replicates and between-conditions) to be included in the analysis. (Default: 0).

-  **-nan** | **--nan-threshold**: Proportion of samples with nan values allowed per condition to calculate a DeltaPSI (Default: 0, i.e. no missing values allowed).

- **-o** | **--output**: Name of the output

- **-h** | **--help**: display the help message describing the different parameters


### **Differential transcript usage** ###

An example of the usage of the program with transcripts is, indicating that replicates are paired (-pa), to apply a multple testing correction (-gc) and perform pairwise comparison between all conditions (-c):

```
python3.4 suppa.py diffSplice --method <empirical> --input <ioi-file> --psi <Cond1.psi> <Cond2.psi> --tpm <Cond1_expression-file> <Cond2_expression-file> --area <1000> --lower-bound <0.05> -pa -gc -c -o <output-file>
```

### **Differential splicing with local events** ###

An example of the usage of the program with local events, applying a multple testing correction (-gc):

```
python3.4 suppa.py diffSplice --method <empirical> --input <ioe-file> --psi <Cond1.psi> <Cond2.psi> --tpm <Cond1_expression-file> <Cond2_expression-file> --area <1000> --lower-bound <0.05> -gc -o <output-file>
```

### **Output files** ###

The differential splicing operation generates a *dpsi* file and a *psivec* file for both cases, differential splicing of events and differential transcript usage.  


### dpsi file

*dpsi* is a tab separated file, with the (transcript or local) event ID in the first column, followed by a variable even number of fields, two for each pair of conditions being compared:

1. **Cond1_Cond2_dPSI**: Event PSI difference (ΔPSI) between Cond1 and Cond2 (ΔPSI = PSI_2 - PSI_1).

2. **Cond1_Cond2_pvalue**: Significance of the difference of PSI between Cond1 and Cond2

An example of a *dpsi* file for transcript events is the following one:
```
Cond1_Cond2_dPSI	Cond1_Cond2_p-val
ENSG00000000419.8;ENST00000413082.1	0.1204556785	0.0119388277
ENSG00000000419.8;ENST00000371583.5	0.0940010396	0.0897896221
ENSG00000000419.8;ENST00000494752.1	0.1904161332	0.2039492118
```

An example of a *dpsi* file for local events is the following one:
```
Cond1_Cond2_dPSI	Cond1_Cond2_p-val
ENSG00000000419;A3:chr20:49557492-49558568:49557470-49558568:-	0.1307455855	0.0484515485
ENSG00000000419;SE:chr20:49557470-49557642:49557746-49558568:-	0.0469449126	0.0954045953
ENSG00000000419;SE:chr20:49557470-49557666:49557746-49558568:-	0.0164051352	0.1518481518
```
By default, SUPPA calculates ΔPSI values pairwise between each pair of adjacent conditions as provided, and calculating the PSI difference between a given condition and the previous one. For instance, for three conditions 1,2,3, the ΔPSI values will be calculated for 2 - 1 (Cond1_Cond2_dPSI) and 3 - 2 (Cond2_Cond3_dPSI). 
On the other hand, if  the **--combination** flag is indicated, ΔPSI values are calculated pairwise between all the possible combinations of conditions. For the previous example with conditions provided as 1,2,3, the ΔPSI values will then be calculated for 2 - 1 (Cond1_Cond2_dPSI), 3 - 1 (Cond1_Cond3_dPSI) and 3 - 2 (Cond2_Cond3_dPSI).


**psivec file**
-------


*psivec* is a tab separated file, with the (transcript or local) event ID in the first column, followed by a variable number of fields, each giving the PSI value per replicate:


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

Using the relative abundances (PSI values) of (transcript or local) events in all samples, and the information of which events change significantly in at least one comparison, SUPPA calculates clusters of events using a density-based clustering. Density-based clustering has various advantages over other methods: it does not require to choose the number of clusters, as this is driven by the data; and it cluster together events that even though they might not have similar PSI values, they behave similarly across conditions as long as there are sufficient events between them. This function works in the same way for transcripts or for local alternative splicing events.


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
python3.4 suppa.py clusterEvents [options]

```
List of options available:

- **-d** | **--dpsi**: Input .dpsi file generated from SUPPA diffSplice operation.

- **-p** | **--psivec**: Input .psivec file generated from SUPPA diffSplice operation.

- **-st** | **--sig-threshold**: p-value threshold to consider an event significant from the dpsi file.

- **-dt** | **--dpsi-threshold**: Lower-bound for the absolute delta PSI value to cluster. (Default: 0.05)

- **-e** | **--eps**: Maximum distance (between 0 and 1) to consider two events as members of the same cluster. (Default: 0.05).

- **-m** | **--metric**: distance metric. Choices: euclidean, manhattan, cosine. (Default: euclidean).

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
python3.4 suppa.py clusterEvents --dpsi <dpsi-file> --psivec <psivec-file> --sig-threshold <0.05> --eps <0.05> --min-pts <20> --groups <1-3,4-6> -o <output-file> 

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

and similarly for transcript events:
**Note:** The clustvec file doesn't contain a header. The header is for illustraion purpouse only. 
```
event_id	cluster_id	cond1_PSI_avg	cond2_PSI_avg
ENSG00000133961.15;ENST00000556772.1 0    0.728937981178  0.223877291002  
ENSG00000133961.15;ENST00000554818.1 1    0.041992298773  0.508135647228
ENSG00000133961.15;ENST00000556600.1 -1   0.106752966827  0.127736627782
```


The cluster_id is 
- **-1**: For unclustered events
- **0,1,2,...**: For each of the clusters formed. Clusteres are numbered from zero. 



----------------------------
# License
----------------------------

SUPPA is released under the MIT license.
