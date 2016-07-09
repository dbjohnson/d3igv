import random
from xml.etree import ElementTree

import pysam
import requests


BASES = ('A', 'C', 'G', 'T')
UCSC_URL = 'http://genome.ucsc.edu/cgi-bin/das/hg19/dna'

segment_to_reference = {}


class BAM(object):
    def __init__(self, filename):
        self.bam = pysam.AlignmentFile(filename)

    def fetch_segment(self, chrom, start, end):
        reads = self.fetch_reads(chrom, start, end)
        coverage = BAM.pileup(reads)
        try:
            reference_sequence = BAM.fetch_reference_for_segment(chrom, start, end)
        except RuntimeError:
            reference_sequence = BAM.reference_by_consensus(coverage)

        return {'chr': chrom,
                'start': start,
                'end': end,
                'reference': reference_sequence,
                'reads': reads,
                'coverage': coverage}

    def fetch_reads(self, chrom, start, end):
        return [{'sequence': r.query_sequence,
                 'fwd': not r.is_reverse,
                 'quality': r.mapq,
                 'start': r.reference_start}
                for r in self.bam.fetch(chrom, start, end)
                if r.reference_start >= start and r.reference_start + len(r.query_sequence) <= end]

    @staticmethod
    def fetch_reference_for_segment(chrom, start, end):
        segment = '{chrom}:{start},{end}'.format(chrom=chrom, start=start + 1, end=end)
        if segment not in segment_to_reference:
            resp = requests.get(UCSC_URL, params={'segment': segment}, timeout=2)
            if resp.status_code != 200:
                raise RuntimeError("Non-200 response received from UCSC: {}".format(resp))

            e = ElementTree.fromstring(resp.text)
            segment_to_reference[segment] = e.find('SEQUENCE').find('DNA').text.strip().upper()

        return segment_to_reference[segment]

    @staticmethod
    def pileup(reads):
        start = min([r['start'] for r in reads])
        end = max([r['start'] + len(r['sequence']) for r in reads])
        coverage = {pos: {base: 0 for base in BASES}
                    for pos in xrange(start, end)}
        for read in reads:
            for i, base in enumerate(read['sequence']):
                if base in BASES:
                    coverage[read['start'] + i][base] += 1
        return coverage

    @staticmethod
    def reference_by_consensus(pileup):
        reference_sequence = []
        for position in sorted(pileup.keys()):
            coverage = pileup[position]
            reference_sequence.append(max(coverage, key=lambda b: coverage[b]))

        return ''.join(reference_sequence)


#############################################################################


class SimBAM(BAM):
    @staticmethod
    def fetch_segment(read_error_rate=0.03, sequence_length=80, num_reads=100,
                      max_SNPs=4, tumor_content=0.3):
        reference_sequence = [random.choice(BASES)
                              for _ in xrange(sequence_length)]

        SNPs = [SimBAM._make_SNP(reference_sequence, tumor_content, reticle=(i == 0))
                for i in xrange(random.randrange(1, max_SNPs))]

        reads = [SimBAM._make_fake_read(reference_sequence, read_error_rate, SNPs)
                 for _ in xrange(num_reads)]

        return {'chr': 'chr1',
                'start': 0,
                'end': len(reference_sequence),
                'reference': ''.join(reference_sequence),
                'reads': reads,
                'coverage': BAM.pileup(reads)}

    @staticmethod
    def _make_fake_read(reference_sequence, read_error_rate, SNPs):
        read_sequence = reference_sequence[:]
        with_read_errors = SimBAM._induce_read_errors(read_sequence, read_error_rate)
        with_SNPs = SimBAM._induce_SNP_reads(with_read_errors, SNPs)
        return SimBAM._trim_read(with_SNPs, random.choice([True, False]))

    @staticmethod
    def _make_SNP(reference_sequence, tumor_content=0.3,
                  homozygous=random.choice([True, False]),
                  reticle=False):
        pos = len(reference_sequence) / 2 if reticle else random.randrange(len(reference_sequence))
        reference_base = reference_sequence[pos]
        mutant_base = random.choice([b for b in BASES if b != reference_base])
        max_af = tumor_content if homozygous else tumor_content / 2.
        af = (1 + random.random()) / 2 * max_af
        return {'pos': pos,
                'base': mutant_base,
                'af': af}

    @staticmethod
    def _induce_read_errors(read_sequence, read_error_rate):
        return [b if random.random() > read_error_rate
                else random.choice(BASES)
                for b in read_sequence]

    @staticmethod
    def _induce_SNP_reads(read_sequence, SNPs):
        with_SNPs = read_sequence[:]
        for m in SNPs:
            if random.random() < m['af']:
                with_SNPs[m['pos']] = m['base']
        return with_SNPs

    @staticmethod
    def _trim_read(read_sequence, fwd):
        sequence_length = len(read_sequence)
        start = random.randrange(sequence_length * 4 / 5)
        length = (random.randrange(sequence_length) + sequence_length) / 2
        end = min(start + length, sequence_length)
        if not fwd:
            tmp = start
            start = sequence_length - end
            end = sequence_length - tmp
        return {'start': start,
                'sequence': ''.join(read_sequence[start:end + 1]),
                'fwd': fwd}
