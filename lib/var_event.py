import itertools
from lib.event import Event


class VariableEvent(Event):
    """
    PCR2 events - events with thresholds
    """

    def merge_events(self, threshold, overlap_fun):
        """
        Adds new transcripts to existing events.
        function overlap_fun evaluates if two events overlap
        """
        ids = {event: [] for event in self.positive_ids}

        for id1, id2 in itertools.combinations(self.positive_ids, 2):
            if overlap_fun(id1, id2, threshold):
                ids[id1] += [id2]
                ids[id2] += [id1]

        # add new events
        for id1 in ids:
            for id2 in ids[id1]:
                self.positive_ids[id1] += [trans_id for trans_id in self.positive_ids[id2]
                                           if trans_id not in self.positive_ids[id1]]

                self.negative_ids[id1] += [trans_id for trans_id in self.negative_ids[id2]
                                           if trans_id not in self.negative_ids[id1]]

    def export_events_ioe(self):
        """
        Returns a list of possible events
        """
        for event in self.positive_ids:
            pos_trans = ','.join(sorted(self.positive_ids[event]))
            all_trans = ','.join(list(set(sorted(self.positive_ids[event] + self.negative_ids[event]))))
            full_event = '{};{}:{}'.format(self.gene.name, self.etype, event)

            yield ('{}\t{}\t{}\t{}\t{}\n'.format(self.gene.chr, self.gene.name, full_event,
                                                 pos_trans, all_trans),
                   self.etype)


class RIv(VariableEvent):
    """
    Retained intro event. PCR version with threshold
    """
    def __init__(self, gene, threshold):
        Event.__init__(self, gene)
        self.etype = 'RI'

        self.construct_events(gene)

        if len(self.positive_ids) > 1 and threshold:
            self.merge_events(threshold=threshold, overlap_fun=self.overlap)

    def add_putative_events(self, exons, transcript):
        gene = self.gene
        for i in range(len(exons) - 1):
            firstexon = exons[i]
            secondexon = exons[i + 1]
            eventid = '{}:{}-{}:{}'.format(gene.chr, firstexon[1], secondexon[0], gene.strand)

            self.negative_ids[eventid] = self.negative_ids.get(eventid, []) + [transcript]
            search_coord = (firstexon[1], secondexon[0])
            search = self.alt_search.get(search_coord, set())
            search.add(eventid)
            self.alt_search[search_coord] = search

    def add_real_events(self, exons, transcript):
        for i in range(len(exons)):
            search_coord = exons[i]

            for stash_pair in self.alt_search:
                if (search_coord[0] < stash_pair[0]) and (search_coord[1] > stash_pair[1]):
                    eventids = self.alt_search[stash_pair]
                    for eventid in eventids:
                        self.positive_ids[eventid] = self.positive_ids.get(eventid, []) + [transcript]

    def export_events_gtf(self, edge):
        """
        Generator of GTF lines
        """

        strand = self.gene.strand
        for event in self.positive_ids:
            full_event = '{}:{}'.format(self.etype, event)

            e_vals = full_event.replace('-', ':').split(':')

            line1 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[2], strand,
                                           full_event, full_event, 'alternative2')
            yield line1, self.etype

            line2 = self.gtf_string.format(e_vals[3], int(e_vals[3]) + edge, strand,
                                           full_event, full_event, 'alternative2')
            yield line2, self.etype

            line3 = self.gtf_string.format(e_vals[2], e_vals[3], strand, full_event,
                                           full_event, 'alternative1')
            yield line3, self.etype

    @staticmethod
    def overlap(id1, id2, th):
        """
        Returns true if two events match criteria to be merged.
        """

        first = [int(pos) for pos in id1[:-2].replace('-', ':').split(':')[1:]]
        second = [int(pos) for pos in id2[:-2].replace('-', ':').split(':')[1:]]

        if all(map(lambda x: abs(x[0] - x[1]) <= th, zip(first, second))):
            return True
        else:
            return False


class ARSSv(VariableEvent):
    """
    Alternative 'Right' Splice Site. Variable boundary.
    """
    def __init__(self, gene, threshold):
        Event.__init__(self, gene)
        if self.gene.strand == '-':
            self.etype = 'A3'
        else:
            self.etype = 'A5'
        self.positive_coords = {}
        self.negative_coords = {}

        self.construct_events(gene)

        if len(self.positive_ids) > 1 and threshold:
            self.merge_events(threshold=threshold, overlap_fun=self.overlap)

    def add_putative_events(self, exons, transcript):
        for i in range(len(exons) - 1):
            firstexon = exons[i]
            search_coords = firstexon[1]
            self.negative_coords[search_coords] = self.negative_coords.get(search_coords, []) + [transcript]

    def add_real_events(self, exons, transcript):
        for i in range(len(exons) - 1):
            firstexon = exons[i]
            self.create_arss(firstexon[0], firstexon[1], transcript)

    def create_arss(self, coord1, coord2, transcript):
        gene = self.gene
        for neg1 in self.negative_coords:
            if coord1 < neg1 < coord2:
                eventid = '{}:{}:{}:{}'.format(gene.chr, coord2, neg1, gene.strand)
                self.positive_ids[eventid] = self.positive_ids.get(eventid, []) + [transcript]
                self.negative_ids[eventid] = self.negative_coords[neg1]

    def export_events_gtf(self, edge):
        """
        Generator of GTF lines
        """
        strand = self.gene.strand
        for event in self.positive_ids:
            full_event = '{}:{}'.format(self.etype, event)
            e_vals = full_event.replace('-', ':').split(':')

            line1 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[3], strand, full_event, full_event,
                                           'alternative1')
            yield line1, self.etype

            line3 = self.gtf_string.format(int(e_vals[2]) - edge, e_vals[2], strand, full_event, full_event,
                                           'alternative2')
            yield line3, self.etype

    @staticmethod
    def overlap(id1, id2, th):
        """
        Returns true if two events match criteria to be merged.
        """

        first = [int(pos) for pos in id1[:-2].replace('-', ':').split(':')[1:]]
        second = [int(pos) for pos in id2[:-2].replace('-', ':').split(':')[1:]]
        if first[0] == second[0] and abs(first[1] - second[1]) <= th:
            return True
        else:
            return False


class ALSSv(VariableEvent):
    """
    Alternative 'Left' Splice Site. Variable boundary.
    """
    def __init__(self, gene, threshold):
        Event.__init__(self, gene)
        if self.gene.strand == '-':
            self.etype = 'A5'
        else:
            self.etype = 'A3'
        self.positive_coords = {}
        self.negative_coords = {}

        self.construct_events(gene)

        if len(self.positive_ids) > 1 and threshold:
            self.merge_events(threshold=threshold, overlap_fun=self.overlap)

    def add_putative_events(self, exons, transcript):
        for i in range(1, len(exons)):
            firstexon = exons[i]
            search_coords = firstexon[0]
            self.negative_coords[search_coords] = self.negative_coords.get(search_coords, []) + [transcript]

    def add_real_events(self, exons, transcript):
        for i in range(1, len(exons)):
            firstexon = exons[i]
            self.create_alss(firstexon[0], firstexon[1], transcript)

    def create_alss(self, coord1, coord2, transcript):
        gene = self.gene
        for neg1 in self.negative_coords:
            if coord1 < neg1 < coord2:
                eventid = '{}:{}:{}:{}'.format(gene.chr, coord1, neg1, gene.strand)
                self.positive_ids[eventid] = self.positive_ids.get(eventid, []) + [transcript]
                self.negative_ids[eventid] = self.negative_coords[neg1]

    def export_events_gtf(self, edge):
        """
        Generator of GTF lines
        """
        strand = self.gene.strand
        for event in self.positive_ids:
            full_event = '{}:{}'.format(self.etype, event)
            e_vals = full_event.replace('-', ':').split(':')

            line1 = self.gtf_string.format(int(e_vals[2]), int(e_vals[3]) + edge, strand, full_event, full_event,
                                           'alternative1')
            yield line1, self.etype

            line3 = self.gtf_string.format(e_vals[2], int(e_vals[3]) + edge, strand, full_event, full_event,
                                           'alternative2')
            yield line3, self.etype

    @staticmethod
    def overlap(id1, id2, th):
        """
        Returns true if two events match criteria to be merged.
        """
        first = [int(pos) for pos in id1[:-2].replace('-', ':').split(':')[1:]]
        second = [int(pos) for pos in id2[:-2].replace('-', ':').split(':')[1:]]
        if first[0] == second[0] and abs(first[1] - second[1]) <= th:
            return True
        else:
            return False


class SEv(VariableEvent):
    def __init__(self, gene, threshold):
        Event.__init__(self, gene)
        self.etype = 'SE'

        self.construct_events(gene)
        self.clear_redundant()

        if len(self.positive_ids) > 1 and threshold:
            self.merge_events(threshold=threshold, overlap_fun=self.overlap)

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

    @staticmethod
    def overlap(id1, id2, th):
        """
        Returns true if all positions overlap within a certain threshold
        """
        first = [int(pos) for pos in id1[:-2].replace('-', ':').split(':')[1:]]
        second = [int(pos) for pos in id2[:-2].replace('-', ':').split(':')[1:]]
        if all([abs(first[0] - second[0]) <= th,
                first[1] == second[1],
                first[2] == second[2],
                abs(first[3] - second[3]) <= th]):
            return True
        else:
            return False


class MXEv(VariableEvent):
    """
    Mutually exclusive exons event. Variable boundary.
    """
    def __init__(self, gene, threshold):
        Event.__init__(self, gene)
        self.etype = 'MX'
        self.positive_coords = {}
        self.negative_coords = {}

        self.construct_events(gene)

        if len(self.positive_ids) > 1 and threshold:
            self.merge_events(threshold=threshold, overlap_fun=self.overlap)

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

    @staticmethod
    def overlap(id1, id2, th):
        """
        Returns true if two events match criteria to be merged.
        """
        first = [int(pos) for pos in id1[:-2].replace('-', ':').split(':')[1:]]
        second = [int(pos) for pos in id2[:-2].replace('-', ':').split(':')[1:]]
        if all([abs(first[0] - second[0]) <= th,
                abs(first[3] - second[3]) <= th,
                first[1] == second[1],
                first[2] == second[2],
                first[5] == second[5],
                first[6] == second[6]]):
            return True
        else:
            return False


class ALTLv(VariableEvent):
    """
    Alternative on 'left' side of the gene. Variable boundary.
    """
    def __init__(self, gene, threshold):
        Event.__init__(self, gene)
        if self.gene.strand == '-':
            self.etype = 'AL'
        else:
            self.etype = 'AF'
        self.positive_coords = {}
        self.negative_coords = {}

        self.construct_events(gene)

        if len(self.positive_ids) > 1 and threshold:
            self.merge_events(threshold=threshold, overlap_fun=self.overlap)

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

    @staticmethod
    def overlap(id1, id2, th):
        """
        Returns true if two events match criteria to be merged.
        """
        first = [int(pos) for pos in id1[:-2].replace('-', ':').split(':')[1:]]
        second = [int(pos) for pos in id2[:-2].replace('-', ':').split(':')[1:]]
        if all([abs(first[0] - second[0]) <= th,
                abs(first[2] - second[2]) <= th,
                abs(first[3] - second[3]) <= th,
                first[1] == second[1],
                first[4] == second[4]]):
            return True
        else:
            return False


class ALTRv(VariableEvent):
    """
    Alternative on 'right' side of the gene. Variable boundary.
    """
    def __init__(self, gene, threshold):
        Event.__init__(self, gene)
        if self.gene.strand == '-':
            self.etype = 'AF'
        else:
            self.etype = 'AL'
        self.positive_coords = {}
        self.negative_coords = {}

        self.construct_events(gene)

        if len(self.positive_ids) > 1 and threshold:
            self.merge_events(threshold=threshold, overlap_fun=self.overlap)

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

    @staticmethod
    def overlap(id1, id2, th):
        """
        Returns true if two events match criteria to be merged.
        """
        first = [int(pos) for pos in id1[:-2].replace('-', ':').split(':')[1:]]
        second = [int(pos) for pos in id2[:-2].replace('-', ':').split(':')[1:]]
        if all([abs(first[0] - second[0]) <= th,
                abs(first[2] - second[2]) <= th,
                abs(first[5] - second[5]) <= th,
                first[1] == second[1],
                first[4] == second[4]]):
            return True
        else:
            return False
