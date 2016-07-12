import numpy as np
from pysam import AlignmentFile


bam = AlignmentFile('../bam/demo.bam')
rad = 50
min_peak_coverage = 100

with open('../bam/demo.bed', 'w') as ofh:
    for ref in bam.references:
        pileup = {p.pos: p.n
                  for p in bam.pileup(ref)}
        idx_localmax = [p for p in sorted(pileup)
                        if pileup[p] > min_peak_coverage and
                        pileup[p] == max(map(lambda ii: pileup.get(ii, 0),
                                         xrange(p - rad, p + rad + 1)))]
        # find peaks in localmax
        delta = [idx_localmax[i] - idx_localmax[i - 1]
                 for i in xrange(1, len(idx_localmax))]

        peaks = [idx_localmax[i] for i in xrange(len(delta)) if delta[i] > 10]

        def center_of_mass(peak, rad=50):
            x = range(peak - rad, peak + rad + 1)
            return int(round(np.average(x, weights=[pileup.get(xx, 0) for xx in x])))

        centers = [center_of_mass(p) for p in peaks]
        positions = np.array(sorted(pileup.keys()))
        values = np.array([pileup[p] for p in positions])
        floor = 5
        rad = 100
        positions_over_floor = positions[values > floor]
        for c in centers:
            left = positions_over_floor[positions_over_floor > c - rad][0]
            right = positions_over_floor[positions_over_floor < c + rad][-1]
            ofh.write('{}\n'.format('\t'.join(map(str, [ref, left, right]))))
            print ', '.join(map(str, [ref, left, right]))
