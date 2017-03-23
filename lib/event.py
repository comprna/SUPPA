"""
Functions and classes for event creation and outputting.

@author: Miha Skalic
@email: miha.skalic[at]gmail.com
"""
import os


class MGene:
    """
    Class holding gene information. For easy event construction
    """
    def __init__(self, gene, chrn, strand):
        self.name = gene.name
        self.chr = chrn
        self.metadata = gene
        self.strand = strand
        self.sortedTranscripts = []
        self.sortedExons = {}
        self.construct_transcripts()
        self.construct_exons()

    def construct_transcripts(self):
        """
        Adding ordered Transcripts to gene
        """
        if len(self.metadata.sortedTranscripts) < 2:
            return
        self.sortedTranscripts = self.metadata.sortedTranscripts

    def construct_exons(self):
        """
        Adding ordered Exons to transcript (gene)
        """
        for transcript in self.sortedTranscripts:
            self.sortedExons[transcript] = sorted(
                self.metadata.transcripts[transcript].exons)


class Event:
    file = ''

    def __init__(self, gene):
        self.gene = gene
        self.alt_search = {}
        self.positive_ids = {}
        self.negative_ids = {}
        self.etype = ''
        self.gtf_string = ('{}\t{}\texon\t{{}}\t{{}}\t.\t{{}}\t.\tgene_id "{{}}"; transcript_id "{{}}:{{}}";\n'
                           .format(self.gene.chr, Event.file))

    def construct_events(self, gene):
        for transcript in gene.sortedTranscripts:
            self.add_putative_events(gene.sortedExons[transcript], transcript)
        for transcript in gene.sortedTranscripts:
            self.add_real_events(gene.sortedExons[transcript], transcript)

    def export_events_ioe(self):
        """
        Returns a list of possible events
        """
        for event in self.positive_ids:
            pos_trans = ','.join(self.positive_ids[event])
            all_trans = ','.join(list(set(self.positive_ids[event] + self.negative_ids[event])))
            full_event = '{};{}:{}'.format(self.gene.name, self.etype, event)

            yield ('{}\t{}\t{}\t{}\t{}\n'.format(self.gene.chr, self.gene.name, full_event,
                                                 pos_trans, all_trans),
                   self.etype)

    def add_putative_events(self, *_):
        pass

    def add_real_events(self, *_):
        pass


class ARSS(Event):
    """
    Alternative 'Right' Splice Site
    """
    def __init__(self, gene, *_):
        Event.__init__(self, gene)
        if self.gene.strand == '-':
            self.etype = 'A3'
        else:
            self.etype = 'A5'
        self.positive_coords = {}
        self.negative_coords = {}

        self.construct_events(gene)

    def add_putative_events(self, exons, transcript):
        if len(exons) < 2:
            return
        for i in range(len(exons) - 1):
            firstexon = exons[i]
            secondexon = exons[i + 1]
            search_coords = (firstexon[1], secondexon[0])
            self.negative_coords[search_coords] = self.negative_coords.get(search_coords, []) + [transcript]

    def add_real_events(self, exons, transcript):
        if len(exons) < 2:
            return
        for i in range(len(exons) - 1):
            firstexon = exons[i]
            secondexon = exons[i + 1]
            self.create_arss(firstexon[0], firstexon[1], secondexon[0], transcript)

    def create_arss(self, coord1, coord2, coord3, transcript):
        gene = self.gene
        for neg1, neg2 in self.negative_coords:
            if neg2 == coord3 and (coord1 < neg1 < coord2):
                eventid = '{}:{}-{}:{}-{}:{}'.format(gene.chr, coord2, coord3, neg1, neg2, gene.strand)
                self.positive_ids[eventid] = self.positive_ids.get(eventid, []) + [transcript]
                self.negative_ids[eventid] = self.negative_coords[(neg1, neg2)]

    def export_events_gtf(self, edge):
        """
        Generator of GTF lines
        """
        strand = self.gene.strand
        for event in self.positive_ids:
            full_event = '{}:{}'.format(self.etype, event)
            e_vals = full_event.replace('-', ':').split(':')

            line1 = self.gtf_string.format(int(e_vals[4]) - edge, e_vals[2], strand, full_event, full_event,
                                           'alternative1')
            yield line1, self.etype

            line2 = self.gtf_string.format(e_vals[3], int(e_vals[3]) + edge, strand, full_event, full_event,
                                           'alternative1')
            yield line2, self.etype

            line3 = self.gtf_string.format(int(e_vals[4]) - edge, e_vals[4], strand, full_event, full_event,
                                           'alternative2')
            yield line3, self.etype

            line4 = self.gtf_string.format(e_vals[5], int(e_vals[5]) + edge, strand, full_event, full_event,
                                           'alternative2')
            yield line4, self.etype


class ALSS(Event):
    """
    Alternative 'Left' Splice Site
    """
    def __init__(self, gene, *_):
        Event.__init__(self, gene)
        if self.gene.strand == '-':
            self.etype = 'A5'
        else:
            self.etype = 'A3'
        self.positive_coords = {}
        self.negative_coords = {}

        self.construct_events(gene)

    def add_putative_events(self, exons, transcript):
        if len(exons) < 2:
            return
        for i in range(len(exons) - 1):
            firstexon = exons[i]
            secondexon = exons[i + 1]
            search_coords = (firstexon[1], secondexon[0])
            self.negative_coords[search_coords] = self.negative_coords.get(search_coords, []) + [transcript]

    def add_real_events(self, exons, transcript):
        if len(exons) < 2:
            return
        for i in range(len(exons) - 1):
            firstexon = exons[i]
            secondexon = exons[i + 1]
            self.create_alss(firstexon[1], secondexon[0], secondexon[1], transcript)

    def create_alss(self, coord1, coord2, coord3, transcript):
        gene = self.gene
        for neg1, neg2 in self.negative_coords:
            if neg1 == coord1 and (coord2 < neg2 < coord3):
                eventid = '{}:{}-{}:{}-{}:{}'.format(gene.chr, coord1, coord2, neg1, neg2, gene.strand)
                self.positive_ids[eventid] = self.positive_ids.get(eventid, []) + [transcript]
                self.negative_ids[eventid] = self.negative_coords[(neg1, neg2)]

    def export_events_gtf(self, edge):
        """
        Generator of GTF lines
        """
        strand = self.gene.strand
        for event in self.positive_ids:
            full_event = '{}:{}'.format(self.etype, event)
            e_vals = full_event.replace('-', ':').split(':')

            line1 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[2], strand, full_event, full_event,
                                           'alternative1')
            yield line1, self.etype

            line2 = self.gtf_string.format(e_vals[3], int(e_vals[5]) + edge, strand, full_event, full_event,
                                           'alternative1')
            yield line2, self.etype

            line3 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[2], strand, full_event, full_event,
                                           'alternative2')
            yield line3, self.etype

            line4 = self.gtf_string.format(e_vals[5], int(e_vals[5]) + edge, strand, full_event, full_event,
                                           'alternative2')
            yield line4, self.etype


class ALTL(Event):
    """
    Alternative on 'left' side of the gene
    """
    def __init__(self, gene, *_):
        Event.__init__(self, gene)
        if self.gene.strand == '-':
            self.etype = 'AL'
        else:
            self.etype = 'AF'
        self.positive_coords = {}
        self.negative_coords = {}

        self.construct_events(gene)

    def add_putative_events(self, exons, transcript):
        if len(exons) < 2:
            return
        firstexon = exons[0]
        secondexon = exons[1]
        search_coords = (firstexon[0], firstexon[1], secondexon[0])
        self.negative_coords[search_coords] = self.negative_coords.get(search_coords, []) + [transcript]

    def add_real_events(self, exons, transcript):
        if len(exons) < 2:
            return
        firstexon = exons[0]
        secondexon = exons[1]
        self.create_altl(firstexon[0], firstexon[1], secondexon[0], transcript)

    def create_altl(self, coord1, coord2, coord3, transcript):
        """
        Find and create alternative events
        """
        gene = self.gene
        for neg1, neg2, neg3 in self.negative_coords:
            if neg3 == coord3 and coord2 < neg1:
                eventid = '{}:{}:{}-{}:{}:{}-{}:{}'.format(gene.chr, coord1, coord2, coord3,
                                                           neg1, neg2, neg3, gene.strand)
                self.positive_ids[eventid] = self.positive_ids.get(eventid, []) + [transcript]
                self.negative_ids[eventid] = self.negative_coords[(neg1, neg2, neg3)]

    def export_events_gtf(self, edge):
        """
        Generator of GTF lines
        """
        strand = self.gene.strand
        for event in self.positive_ids:
            full_event = '{}:{}'.format(self.etype, event)
            e_vals = full_event.replace('-', ':').split(':')

            line1 = self.gtf_string.format(e_vals[2], e_vals[3], strand, full_event, full_event,
                                           'alternative1')
            yield line1, self.etype

            line2 = self.gtf_string.format(e_vals[4], int(e_vals[4]) + edge, strand, full_event, full_event,
                                           'alternative1')
            yield line2, self.etype

            line3 = self.gtf_string.format(e_vals[5], e_vals[6], strand, full_event, full_event,
                                           'alternative2')
            yield line3, self.etype

            line4 = self.gtf_string.format(int(e_vals[7]), int(e_vals[7]) + edge, strand, full_event, full_event,
                                           'alternative2')
            yield line4, self.etype


class ALTR(Event):
    """
    Alternative on 'Right' side of the gene
    """
    def __init__(self, gene, *_):
        Event.__init__(self, gene)
        if self.gene.strand == '-':
            self.etype = 'AF'
        else:
            self.etype = 'AL'
        self.positive_coords = {}
        self.negative_coords = {}

        self.construct_events(gene)

    def add_putative_events(self, exons, transcript):
        if len(exons) < 2:
            return
        firstexon = exons[-2]
        secondexon = exons[-1]
        search_coords = (firstexon[1], secondexon[0], secondexon[1])
        self.negative_coords[search_coords] = self.negative_coords.get(search_coords, []) + [transcript]

    def add_real_events(self, exons, transcript):
        if len(exons) < 2:
            return
        firstexon = exons[-2]
        secondexon = exons[-1]
        self.create_altr(firstexon[1], secondexon[0], secondexon[1], transcript)

    def create_altr(self, coord1, coord2, coord3, transcript):
        """
        Find and create alternative events
        """
        gene = self.gene
        for neg1, neg2, neg3 in self.negative_coords:
            if neg1 == coord1 and neg2 > coord3:
                eventid = '{}:{}-{}:{}:{}-{}:{}:{}'.format(gene.chr, coord1, coord2, coord3,
                                                           neg1, neg2, neg3, gene.strand)
                self.positive_ids[eventid] = self.positive_ids.get(eventid, []) + [transcript]
                self.negative_ids[eventid] = self.negative_coords[(neg1, neg2, neg3)]

    def export_events_gtf(self, edge):
        """
        Generator of GTF lines
        """
        strand = self.gene.strand
        for event in self.positive_ids:
            full_event = '{}:{}'.format(self.etype, event)
            e_vals = full_event.replace('-', ':').split(':')

            line1 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[2], strand, full_event, full_event,
                                           'alternative2')
            yield line1, self.etype

            line2 = self.gtf_string.format(e_vals[3], e_vals[4], strand, full_event, full_event,
                                           'alternative2')
            yield line2, self.etype

            line3 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[2], strand, full_event, full_event,
                                           'alternative1')
            yield line3, self.etype

            line4 = self.gtf_string.format(e_vals[6], e_vals[7], strand, full_event, full_event,
                                           'alternative1')
            yield line4, self.etype


class MXE(Event):
    """
    Mutually exclusive exons event
    """
    def __init__(self, gene, *_):
        Event.__init__(self, gene)
        self.etype = 'MX'
        self.positive_coords = {}
        self.negative_coords = {}

        self.construct_events(gene)

    def add_putative_events(self, exons, transcript):
        if len(exons) < 3:
            return
        for i in range(len(exons) - 2):
            firstexon = exons[i]
            midexon = exons[i + 1]
            lastexon = exons[i + 2]
            search_coords = firstexon[1], midexon[0], midexon[1], lastexon[0]
            self.negative_coords[search_coords] = self.negative_coords.get(search_coords, []) + [transcript]

    def add_real_events(self, exons, transcript):
        if len(exons) < 3:
            return
        for i in range((len(exons) - 2)):
            firstexon = exons[i]
            midexon = exons[i + 1]
            lastexon = exons[i + 2]
            search_coords = firstexon[1], midexon[0], midexon[1], lastexon[0]

            for eventid, neg_coord in self.coord_match(search_coords):
                self.negative_ids[eventid] = self.negative_ids.get(eventid, []) + [transcript]
                self.positive_ids[eventid] = self.negative_coords[neg_coord]

    def coord_match(self, new_coords):
        """
        Creates and returns event IDs that match tupled with
        alternative IDs
        """
        all_events = []
        gene = self.gene

        for orig_cord in self.negative_coords:
            compare_f = orig_cord[0] + orig_cord[-1]
            compare_t = new_coords[0] + new_coords[-1]
            mm_f = orig_cord[2]
            mm_t = new_coords[1]

            if compare_f == compare_t and mm_f < mm_t:
                eventid = ('{}:{}-{}:{}-{}:{}-{}:{}-{}:{}'
                           .format(gene.chr, orig_cord[0], orig_cord[1], orig_cord[2], orig_cord[3],
                                   new_coords[0], new_coords[1], new_coords[2], new_coords[3], gene.strand))
                all_events.append([eventid, orig_cord])
        return all_events

    def export_events_gtf(self, edge):
        """
        Generator of GTF lines
        """
        strand = self.gene.strand
        for event in self.positive_ids:
            full_event = '{}:{}'.format(self.etype, event)
            e_vals = full_event.replace('-', ':').split(':')

            line1 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[2], strand, full_event, full_event,
                                           'alternative1')
            yield line1, self.etype

            line2 = self.gtf_string.format(e_vals[3], e_vals[4], strand, full_event, full_event,
                                           'alternative1')
            yield line2, self.etype

            line3 = self.gtf_string.format(e_vals[9], int(e_vals[9]) + edge, strand, full_event, full_event,
                                           'alternative1')
            yield line3, self.etype

            line4 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[2], strand, full_event, full_event,
                                           'alternative2')
            yield line4, self.etype

            line5 = self.gtf_string.format(e_vals[7], e_vals[8], strand, full_event, full_event,
                                           'alternative2')
            yield line5, self.etype

            line6 = self.gtf_string.format(e_vals[9], int(e_vals[9]) + edge, strand, full_event, full_event,
                                           'alternative2')
            yield line6, self.etype


class RI(Event):
    """
    Retained intro event
    """
    def __init__(self, gene, *_):
        Event.__init__(self, gene)
        self.etype = 'RI'

        self.construct_events(gene)

    def add_putative_events(self, exons, transcript):
        gene = self.gene
        for i in range(len(exons) - 1):
            firstexon = exons[i]
            secondexon = exons[i + 1]
            eventid = '{}:{}:{}-{}:{}:{}'.format(gene.chr, firstexon[0], firstexon[1], secondexon[0],
                                                 secondexon[1], gene.strand)

            self.negative_ids[eventid] = self.negative_ids.get(eventid, []) + [transcript]
            search_coord = (firstexon[0], secondexon[1])
            search = self.alt_search.get(search_coord, set())
            search.add(eventid)
            self.alt_search[search_coord] = search

    def add_real_events(self, exons, transcript):
        for i in range(len(exons)):
            search_coord = exons[i]

            if search_coord in self.alt_search:
                eventids = self.alt_search[search_coord]
                for eventid in eventids:
                    self.positive_ids[eventid] = self.positive_ids.get(eventid, []) + [transcript]

    def export_events_gtf(self, *_):
        """
        Generator of GTF lines
        """
        strand = self.gene.strand
        for event in self.positive_ids:
            full_event = '{}:{}'.format(self.etype, event)
            e_vals = full_event.replace('-', ':').split(':')

            line1 = self.gtf_string.format(e_vals[2], e_vals[3], strand, full_event, full_event, 'alternative2')
            yield line1, self.etype

            line2 = self.gtf_string.format(e_vals[4], e_vals[5], strand, full_event, full_event, 'alternative2')
            yield line2, self.etype

            line3 = self.gtf_string.format(e_vals[2], e_vals[5], strand, full_event, full_event, 'alternative1')
            yield line3, self.etype


class SE(Event):
    """
    Skipped exon event
    """
    def __init__(self, gene, *_):
        Event.__init__(self, gene)
        self.etype = 'SE'

        self.construct_events(gene)
        self.clear_redundant()

    def add_putative_events(self, exons, transcript):
        if len(exons) < 3:
            return
        gene = self.gene
        for i in range(len(exons) - 2):
            firstexon = exons[i]
            midexon = exons[i + 1]
            lastexon = exons[i + 2]
            eventid = '{}:{}-{}:{}-{}:{}'.format(gene.chr, firstexon[1], midexon[0],
                                                 midexon[1], lastexon[0], gene.strand)

            self.positive_ids[eventid] = self.positive_ids.get(eventid, []) + [transcript]
            search_coord = (firstexon[1], lastexon[0])
            search = self.alt_search.get(search_coord, set())
            search.add(eventid)
            self.alt_search[search_coord] = search

    def add_real_events(self, exons, transcript):
        if len(exons) < 2:
            return
        for i in range((len(exons)-1)):
            firstexon = exons[i]
            secondexon = exons[i + 1]
            search_coord = (firstexon[1], secondexon[0])

            if search_coord in self.alt_search:
                eventids = self.alt_search[search_coord]
                for eventid in eventids:
                    self.negative_ids[eventid] = self.negative_ids.get(eventid, []) + [transcript]

    def export_events_gtf(self, edge):
        """
        Generator of GTF lines
        """
        strand = self.gene.strand
        for event in self.positive_ids:
            full_event = '{}:{}'.format(self.etype, event)
            e_vals = full_event.replace('-', ':').split(':')

            line1 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[2], strand, full_event, full_event,
                                           'alternative2')
            yield line1, self.etype

            line2 = self.gtf_string.format(e_vals[5], int(e_vals[5]) + edge, strand, full_event, full_event,
                                           'alternative2')
            yield line2, self.etype

            line3 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[2], strand, full_event, full_event,
                                           'alternative1')
            yield line3, self.etype

            line4 = self.gtf_string.format(e_vals[3], e_vals[4], strand, full_event, full_event, 'alternative1')
            yield line4, self.etype

            line5 = self.gtf_string.format(e_vals[5], int(e_vals[5]) + edge, strand, full_event, full_event,
                                           'alternative1')
            yield line5, self.etype

    def clear_redundant(self):
        for id_key in list(self.positive_ids.keys()):
            if id_key not in self.negative_ids:
                self.positive_ids.pop(id_key)


class EWriter:
    def __init__(self, all_events, output_name, etype, boundary):
        """
        Creates 'write to' buckets
        """
        self.wdict = {}
        for event in all_events:
            if event == 'SS':
                self.wdict['A3'] = open('{}_A3_{}.{}'.format(output_name, boundary, etype), 'w')
                self.wdict['A5'] = open('{}_A5_{}.{}'.format(output_name, boundary, etype), 'w')
            elif event == 'FL':
                self.wdict['AF'] = open('{}_AF_{}.{}'.format(output_name, boundary, etype), 'w')
                self.wdict['AL'] = open('{}_AL_{}.{}'.format(output_name, boundary, etype), 'w')
            else:
                self.wdict[event] = open('{}_{}_{}.{}'.format(output_name, event, boundary, etype), 'w')
        for bucket in self.wdict:
            if etype == 'ioe':
                self.wdict[bucket].write("seqname\tgene_id\tevent_id\talternative_transcripts\ttotal_transcripts\n")
            if etype == 'gtf':
                self.wdict[bucket].write("track name={} visibility=2\n".format(bucket))

    def write(self, line, bucket):
        """
        Given writes the line to the given bucket
        """
        self.wdict[bucket].write(line)

    def close(self):
        """
        Closes all files in writer
        """
        for idx in self.wdict:
            self.wdict[idx].close()


def process_events(my_gene, event, ioe_writer, gtf_writer, edge_len, th):
    """
    Generates specified event occurrences in gene and writes the GTF/IOE output
    """
    gene_event = event(my_gene, th)
    for event_repr, etype in gene_event.export_events_ioe():
        ioe_writer.write(event_repr, etype)
    for gtf_line, etype in gene_event.export_events_gtf(edge_len):
        gtf_writer.write(gtf_line, etype)

# to avoid circular dependencies
from lib.var_event import *


def create_event_classes(all_events, b_type):
    """
    Creates list of classes that will be calculated for
    """
    if b_type == 'S':
        options = {'SE': SE, 'RI': RI, 'MX': MXE}
    else:
        options = {'SE': SEv, 'RI': RIv, 'MX': MXEv}
    all_classes = []
    for event in all_events:
        if event == 'SS':
            all_classes.append(ARSS if b_type == 'S' else ARSSv)
            all_classes.append(ALSS if b_type == 'S' else ALSSv)
        elif event == 'FL':
            all_classes.append(ALTL if b_type == 'S' else ALTLv)
            all_classes.append(ALTR if b_type == 'S' else ALTRv)
        else:
            all_classes.append(options[event])
    return all_classes


def make_events(events, my_genome, input_name, output_name, edge_len, logger, b_type, th):
    """
    Controls event creation and writing for each gene
    """
    boundary = 'variable_{}'.format(th) if b_type == 'V' else 'strict'
    Event.file = os.path.basename(input_name)
    my_ioe_writer = EWriter(events, output_name, 'ioe', boundary)
    my_gtf_writer = EWriter(events, output_name, 'gtf', boundary)
    my_events = create_event_classes(events, b_type)
    logger.info("Calculating events")
    for gene, xchr, strand in my_genome:
        my_gene = MGene(gene, xchr, strand)
        logger.debug("Analyzing gene: {}".format(gene.name))
        for event in my_events:
            process_events(my_gene, event, my_ioe_writer, my_gtf_writer, edge_len, th)
    my_ioe_writer.close()
    my_gtf_writer.close()
