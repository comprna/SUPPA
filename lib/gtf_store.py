"""
Functions and Classes to parse and store GTF file data


@author: Miha Skalic
@email: miha.skalic[at]gmail.com
"""

import sys
import os


class Genome():
    def __init__(self):
        self.data = {}

    def add_to_genes(self, t_dict):
        """
        Appends found exon to his place in genome
        """
        # and strand if and chr if not present
        if t_dict['chr'] not in self.data:
            self.data[t_dict['chr']] = {}
        if t_dict['strand'] not in self.data[t_dict['chr']]:
            self.data[t_dict['chr']][t_dict['strand']] = {}

        # append to gene:
        my_gene = self.data[t_dict['chr']][t_dict['strand']].get(t_dict['gene'], Gene(t_dict['gene']))
        my_gene.add_transcript(t_dict)
        self.data[t_dict['chr']][t_dict['strand']][t_dict['gene']] = my_gene

    def __iter__(self):
        """
        Generates tuples: Gene object, chr str, strand str
        """
        for chr_p in self.data:
            for strand in self.data[chr_p]:
                for gene in self.data[chr_p][strand]:
                    yield self.data[chr_p][strand][gene], chr_p, strand

    def sort_transcripts(self):
        """
        For each gene in genome sorts transcripts (for further
        precessing).
        """
        for gene, _, _ in self.__iter__():
            gene.sort_transcripts()

    def split_genes(self):
        """
        Checks Transkript overlap and split gene if transcripts do not overlap
        """
        # new list consists of pairs (chr, strand, gene_name, list of first transcript indencies)
        new_genes = []

        # construct new genes
        for gene, xchr, strand in self.__iter__():
            overlaps = gene.get_non_overlaps()
            if overlaps:
                new_gene_group = gene.split_gene(overlaps)
                new_genes.append((xchr, strand, gene.name, new_gene_group))

        # delete old gene and add new genes
        for xchr, strand, gname, xgenes in new_genes:
            self.data[xchr][strand].pop(gname)
            for ngene in xgenes:
                self.data[xchr][strand][ngene.name] = ngene

    def merge_two_genes(self, gene1, gene2, chr_p, strand):
        """
        deletes the old gene and creates new merged gene
        """
        # construct new gene
        merge_g = Gene('{}_and_{}'.format(gene1.name, gene2.name))
        merge_g.transcripts = gene1.transcripts
        merge_g.transcripts.update(gene2.transcripts)
        merge_g.sort_transcripts()

        # add new genes old gene
        self.data[chr_p][strand][merge_g.name] = merge_g

        # delte old genes
        self.data[chr_p][strand].pop(gene1.name)
        self.data[chr_p][strand].pop(gene2.name)
        span = merge_g.get_span()
        return span[0], span[1], merge_g

    def merge_overlaps(self, g_list, chr_p, strand):
        """
        Given a ordered list of genes creates new polled genes.
        """
        current_g = g_list.pop(0)
        while g_list:
            new_g = g_list.pop(0)
            # check if they overlap
            if new_g[0] <= current_g[1]:
                current_g = self.merge_two_genes(current_g[2], new_g[2], chr_p, strand)
            else:
                current_g = new_g

    def poll_genes(self):
        """
        Pools genes that overlap
        """

        for chr_p in self.data:
            for strand in self.data[chr_p]:
                # create list of genes
                g_list = []
                for gene in self.data[chr_p][strand]:
                    coords = self.data[chr_p][strand][gene].get_span()
                    g_list.append((coords[0], coords[1], self.data[chr_p][strand][gene]))
                g_list.sort()
                self.merge_overlaps(g_list, chr_p, strand)


class Gene():
    def __init__(self, name):
        self.name = name
        self.transcripts = {}
        self.sortedTranscripts = []

    def __lt__(self, other):
        return 0

    def add_transcript(self, t_dict):
        transcript = self.transcripts.get(t_dict['transcript'], Transcript())
        transcript.add_exon(t_dict['coordinates'][0], t_dict['coordinates'][1])
        self.transcripts[t_dict['transcript']] = transcript

    def sort_transcripts(self):
        self.sortedTranscripts = sorted(
            self.transcripts.keys(),
            key=lambda x: self.transcripts[x].span)

    def get_non_overlaps(self):
        """
        Returns indencies of transcripts that do not overlap with previous transcript
        """
        new_sites = []
        m_upper = self.transcripts[self.sortedTranscripts[0]].span[1]
        for i in range(1, len(self.sortedTranscripts)):
            if self.transcripts[self.sortedTranscripts[i]].span[0] > m_upper:
                    # self.transcripts[self.sortedTranscripts[i - 1]].span[1]:
                new_sites.append(i)
            m_upper = max([m_upper, self.transcripts[self.sortedTranscripts[i]].span[1]])
        return new_sites

    def split_gene(self, idxs):
        """
        Splits gene into two new genes at give incendies.
        Returns two new gene objects
        """
        new_genes = []
        idxs = [0] + idxs + [len(self.sortedTranscripts)]
        for i in range(1, len(idxs)):
            gname = '{}_locus{}'.format(self.name, i)
            new_gene = Gene(gname)
            new_gene.sortedTranscripts = self.sortedTranscripts[idxs[i-1]:idxs[i]]
            new_gene.transcripts = {xtranscript: self.transcripts[xtranscript]
                                    for xtranscript in new_gene.sortedTranscripts}
            new_genes.append(new_gene)
        return new_genes

    def get_span(self):
        """
        Fetches the gene range span and returns results in tuple
        """
        return min([self.transcripts[transc].span[0] for transc in self.transcripts]), \
            max([self.transcripts[transc].span[1] for transc in self.transcripts])


class Transcript():
    def __init__(self):
        self.exons = set()
        self.span = (float('inf'), float('-inf'))

    def add_exon(self, start, stop):
        """
        Adds exon to transkript
        """
        self.exons.add((start, stop))
        self.span = (min([self.span[0], start]), max([self.span[1], stop]))


def parse_gtf_line(line, l_count, logger):
    #  Comment line
    if line.startswith('#'):
        return
    line = line.strip().split('\t')

    # missformated line
    if len(line) != 9:
        logger.info('Missmatch in number of Fields. skipping line {}'.format(l_count))
        return
    # Skip if it not an exon
    if line[2] != 'exon':
        return
    # parsing the attributes
    attributes = [att.split(' ', 1) for att in line[-1].strip(';').split('; ')]
    #Some gtf has an extra space in column 9 (like GRCh37 v67 Ensembl). Remove this space if it's in there and reformat the attribute
    if(attributes[0][0]==""):
        attributes[0] = attributes[0][1].split(" ")
    att_dict = dict(map(lambda x: (x[0], x[1].strip('"')), attributes))
    if 'gene_id' not in att_dict or 'transcript_id' not in att_dict:
        logger.info('Missing gene_id or transcript_id field in line {}'.format(l_count))
        return

    return {'gene': att_dict['gene_id'],
            'transcript': att_dict['transcript_id'],
            'coordinates': (int(line[3]), int(line[4])),
            'strand': line[6],
            'chr': line[0]}


def gtf_reader(parse_file, logger):
    """
    Parses the input GTF file and creates exon
    representations for further processing
    """
    l_count = 0
    all_results = []
    # Check if file exists
    if not os.path.isfile(parse_file):
        sys.stderr.write("Input file does not exist. Quiting\n")
        exit(1)
    # logger.info("Reading input data.")
    with open(parse_file, 'r') as handle:
        for line in handle:
            l_count += 1
            results = parse_gtf_line(line, l_count, logger)
            if results:
                all_results.append(results)
    return all_results
