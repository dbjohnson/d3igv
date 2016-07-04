import random

BASES = ['A', 'C', 'G', 'T']


def make_SNP(reference_sequence, tumor_content=0.3, homozygous=random.choice([True, False])):
    max_af = tumor_content if homozygous else tumor_content / 2.
    af = (1 + random.random()) / 2 * max_af
    return {'pos': random.randrange(len(reference_sequence)),
            'base': random.choice(BASES),
            'af': af}


def induce_read_errors(read_sequence, read_error_rate):
    return [b if random.random() > read_error_rate
            else random.choice(BASES)
            for b in read_sequence]


def induce_SNP_reads(read_sequence, SNPs):
    with_SNPs = read_sequence[:]
    for m in SNPs:
        if random.random() < m['af']:
            with_SNPs[m['pos']] = m['base']
    return with_SNPs


def trim_read(read_sequence, fwd):
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


def make_fake_read(reference_sequence, read_error_rate, SNPs):
    read_sequence = reference_sequence[:]
    with_read_errors = induce_read_errors(read_sequence, read_error_rate)
    with_SNPs = induce_SNP_reads(with_read_errors, SNPs)
    return trim_read(with_SNPs, random.choice([True, False]))


def generate_fake_reads(read_error_rate=0.03, sequence_length=80, num_reads=100,
                        max_SNPs=4, tumor_content=0.3):
    reference_sequence = [random.choice(BASES) for _ in xrange(sequence_length)]
    n_SNPs = random.randrange(max_SNPs)
    SNPs = [make_SNP(reference_sequence, tumor_content=tumor_content)
            for _ in xrange(n_SNPs)]

    reticle_SNP = make_SNP(reference_sequence, tumor_content)
    reticle_SNP['pos'] = sequence_length / 2
    SNPs.append(reticle_SNP)

    reads = [make_fake_read(reference_sequence, read_error_rate, SNPs)
             for _ in xrange(num_reads)]

    coverage = [{base: 0 for base in BASES}
                for _ in xrange(len(reference_sequence))]
    for read in reads:
        for i, base in enumerate(read['sequence']):
            coverage[i + read['start']][base] += 1

    return {'reference': ''.join(reference_sequence),
            'reads': reads,
            'coverage': coverage}
