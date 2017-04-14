#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Everything related to extraction of reads in the BAM file - like depth
computation, computing paired-end distances, classify spanning reads and read
length extraction. These signals are then pushed into the Maximum Likelihood
model for the prediction of allele sizes.
"""

import logging
import math
import sys
import string

import numpy as np
import pysam

from collections import defaultdict
from ssw import Aligner
from utils import datafile


SPAN = 1000
FLANKMATCH = 9
DNAPE_ELONGATE = SPAN * 10  # How far do we look beyond the target for paired-end
_complement = string.maketrans('ATCGatcgNnXx', 'TAGCtagcNnXx')


class BamParser:
    '''
    Find TRED repeats from aligned reads bam file
    :inputParams: InputParams object
    '''
    def __init__(self, inputParams):
        self.inputParams = inputParams
        self.logger = logging.getLogger('BamParser')
        self.logger.setLevel(inputParams.getLogLevel())
        self.bam = inputParams.bam
        self.gender = inputParams.gender
        self.depth = inputParams.depth
        self.READLEN = inputParams.READLEN

        # initialize tred-specific things
        self.tred = inputParams.tred
        self.repeatSize = len(self.tred.repeat)
        self.chr = self.tred.chr

        # X-linked TRED
        if self.gender == 'Male' and self.tred.is_xlinked:
            self.ploidy = 1
        else:
            self.ploidy = self.tred.ploidy
        self.logger.debug("Locus: {}; Ploidy={}"\
                        .format(str(self.tred), self.ploidy))

        self.repeat = self.tred.repeat
        self.startRepeat, self.endRepeat = self.tred.repeat_start, self.tred.repeat_end
        self.referenceLen = self.tred.repeat_end - self.tred.repeat_start + 1
        self.fullPrefix, self.fullSuffix = self.tred.prefix, self.tred.suffix

        # Compute REPT cutoff
        self.period = len(self.repeat)
        self.max_units = int(math.ceil(self.READLEN * 1. / self.period))

        # Stores all the read counts for each repeat units
        self.counts = {}
        self.counts["PREF"] = self.counts["POST"] = defaultdict(int)
        for tag in ("FULL", "REPT", "HANG"):
            self.counts[tag] = defaultdict(int)

        self.details = []  # Store read sequences, enabled on logging.INFO

    def _buildDB(self):
        '''
        Build a series of aligners that each uses a reference with varying
        number of repeats - whichever scores the best is the winner.
        '''
        ssws = []
        # Build a list of ssw
        for units in xrange(1, self.max_units + 1):
            target = self.fullPrefix + self.repeat * units + self.fullSuffix
            target_rc = rc(target)
            for seq in (target, target_rc):
                ssw = Aligner(ref_seq=seq,
                              match=1, mismatch=5, gap_open=7, gap_extend=2,  # Strict
                              #match=1, mismatch=4, gap_open=6, gap_extend=1, # BWA-MEM
                              report_secondary=False)
                ssws.append((units, seq, ssw))
        return ssws

    def get_hangs(self, al):
        """
        Determine the type of overlap given query, ref alignment coordinates
        Consider the following alignment between sequence a and b:
        aLhang \              / aRhang
                \------------/
                /------------\
        bLhang /              \ bRhang
        Terminal overlap: a before b, b before a
        Contain overlap: a in b, b in a
        """
        aLhang, aRhang = al.ref_begin, len(al.ref_seq) - al.ref_end - 1
        bLhang, bRhang = al.query_begin, len(al.query_seq) - al.query_end - 1

        s1 = aRhang + bLhang
        s2 = aLhang + bRhang
        s3 = aLhang + aRhang
        s4 = bLhang + bRhang

        return min(s1, s2, s3, s4)

    def _parseReadSW(self, chr, seq, db, verbose=False):
        '''
        Use Smith-Waterman matcher to classify reads and count number of
        repeats. This is the preferred method that allows mismatches (sequencing
        errors or SNPs) inside the read.
        '''
        res = []
        for units, target, ssw in db:
            min_len = min(len(seq), len(target)) / 2
            min_score = max(min_len, 30)
            al = ssw.align(seq, min_score=min_score, min_len=min_len)
            if not al:
                continue

            prefix_read = al.ref_begin < FLANKMATCH
            suffix_read = al.ref_end > len(target) - FLANKMATCH - 1
            hang = self.get_hangs(al)
            hang_read = hang >= FLANKMATCH

            if verbose:
                print >> sys.stderr, units, target
                print >> sys.stderr, str(al).strip()
                print >> sys.stderr, '\n'.join(al.alignment)
                print >> sys.stderr, prefix_read, suffix_read, hang_read, hang
                print >> sys.stderr

            if hang_read:
                tag = "HANG"
            elif prefix_read:
                tag = "PREF"
                if suffix_read:
                    tag = "FULL"
            elif suffix_read:
                tag = "POST"
            elif units >= self.max_units - 1 and units * self.period <= len(seq):
                tag = "REPT"
            else:
                continue
            res.append((al.score, units, tag))

        if not res:
            return

        score, h, tag = max(res, key=lambda x: (x[0], -x[1]))
        self.counts[tag][h] += 1

        s = "{}: h={:>3}, seq={}".format(tag, h, seq)
        self.logger.debug(s)

        if tag == "HANG":
            return
        self.details.append({'tag': tag, 'h': h, 'seq': seq})

    def parse(self, pad=SPAN):
        """
        We examine all reads within interval (WINDOW_START, WINDOW_END). We are
        only interested in two types of reads
        - Unmapped reads (they are really "half-mapped" reads, with the mate
          acting as an anchor)
        - Reads that are aligned close to the repeat region (within distance of
          a read length)
        """
        # TODO: limit by chromosome end
        WINDOW_START = max(0, self.startRepeat - pad)
        WINDOW_END = self.endRepeat + pad
        READ_START = max(0, self.startRepeat - self.READLEN)
        READ_END = self.endRepeat + self.READLEN

        samfile = read_alignment(self.bam)
        db = self._buildDB()

        chr, start, end = self.chr, WINDOW_START, WINDOW_END
        n_unmapped = 0
        if test_fetch(samfile, chr, start, end, self.logger):
            for read in samfile.fetch(chr, start, end):
                if read.is_unmapped:
                    n_unmapped += 1
                else:
                    if read.reference_start < READ_START:
                        continue
                    if read.reference_start > READ_END:
                        continue
                self._parseReadSW(chr=chr, seq=read.query_sequence, db=db)

        self.logger.debug("A total of {} unmapped reads in {}:{}-{}".\
                            format(n_unmapped, chr, start, end))

        for tag in ("FULL", "PREF", "REPT"):
            self.show_counts(tag)

    def show_counts(self, tag):
        countMap = self.counts[tag]
        total = sum(v for v in countMap.values())
        counts = ["{}:{}".format(k, v) for (k, v) in sorted(countMap.items())]
        counts = " ".join(counts)
        self.logger.debug("Counts [{}] (total={}) => {}".format(tag, total, counts))


class BamParserResults:
    '''
    Encapsulates all results: counts from BamParser and calls from different callers
    '''
    def __init__(self, inputParams, tred, counts, details, caller):
        self.inputParams = inputParams
        self.tred = tred
        self.details = details
        self.counts = counts
        self.FDP = sum(counts["FULL"].values())
        self.PDP = sum(counts["PREF"].values())
        self.RDP = max(counts["REPT"].values()) if counts["REPT"] else 0
        self.PEDP = caller.PEDP
        self.PEG = caller.PEG
        self.PET = caller.PET
        self.CI = caller.CI
        self.PP = caller.PP
        self.label = caller.label
        self.alleles = caller.alleles
        self.P_h1 = caller.P_h1
        self.P_h2 = caller.P_h2
        self.P_h1h2 = caller.P_h1h2
        self.P_PEG = caller.P_PEG
        self.P_PET = caller.P_PET


class PEextractor:
    """
    Infer distance paired-end reads spanning a certain region.
    """
    def __init__(self, bp):
        samfile = read_alignment(bp.bam)
        chr = bp.chr
        start = bp.startRepeat
        end = bp.endRepeat
        self.ref = bp.referenceLen

        # Compute the target distribution (defined as paired spanning the CAG repeats)
        pstart = max(start - DNAPE_ELONGATE, 0)
        pend = end + DNAPE_ELONGATE
        cache = {}
        if test_fetch(samfile, chr, pstart, pend, bp.logger):
            cache = defaultdict(list)
            for x in samfile.fetch(chr, pstart, pend):
                if not x.is_paired:
                    continue
                if x.is_unmapped:
                    continue
                if x.is_duplicate:
                    continue
                cache[x.query_name].append(x)

        self.global_lens, self.target_lens = [], []
        tstart = start - FLANKMATCH
        tend = end + FLANKMATCH
        for name, reads in cache.iteritems():
            if len(reads) < 2:
                continue
            a, b = reads[:2]
            if not ((not a.is_reverse) and b.is_reverse):  # Mapped in +, - orientation
                continue

            tlen = self.get_target_length(a, b)
            if tlen >= SPAN:  # Skip pairs that are too distant
                continue
            # Get all pairs where read1 is on left flank and read2 is on right flank (spanning pair)
            if a.reference_start < tstart and b.reference_end > tend:
                self.target_lens.append(tlen)
            else:
                self.global_lens.append(tlen)

        self.MINPE = end - start + 2 * FLANKMATCH + 2

    def get_target_length(self, a, b):
        start, end = a.reference_start, b.reference_end
        if a.query_alignment_start > 0:  # has clips
            start -= a.query_alignment_start
        if b.query_alignment_end < b.query_length:  # has clips
            end += b.query_length - b.query_alignment_end
        return end - start


class BamReadLen:
    """
    Returns the read length in BAM file.
    """
    def __init__(self, bamfile, logger):
        self.bamfile = bamfile
        self.logger = logger

    @property
    def readlen(self, firstN=100):
        sam = read_alignment(self.bamfile)
        rls = []
        for read in sam.fetch():
            rls.append(read.query_length)
            if len(rls) > firstN:
                break
        rmin, rmax = min(rls), max(rls)
        if rmin != rmax:
            self.logger.debug("Read length: min={}bp max={}bp".format(rmin, rmax))
        return rmax


class BamDepth:
    """
    Computes the average depth for a particular region, both for the inference
    of the repeat size and inference of gender.
    """
    def __init__(self, bamfile, ref, logger):
        self.bamfile = bamfile
        self.logger = logger
        self.ref = ref

    def region_depth(self, chr, start, end, verbose=False):
        sam = read_alignment(self.bamfile)
        depths = [c.n for c in sam.pileup(chr, start, end)]
        depth = sum(depths) * 1. / (end - start + 1)
        if verbose:
            self.logger.debug("Depth of region {}:{}-{}: {}"\
                            .format(chr, start, end, depth))
        return depth

    def get_Y_depth(self, N=5):
        UNIQY = datafile("chrY.{}.unique_ccn.gc".format(self.ref))
        fp = open(UNIQY)
        depths = []
        for i, row in enumerate(fp):
            # Some regions still have mapped reads, exclude a few
            if i in (1, 4, 6, 7, 10, 11, 13, 16, 18, 19):
                continue
            if len(depths) >= N:
                break
            c, start, end, gc = row.split()
            start, end = int(start), int(end)
            d = self.region_depth(c, start, end)
            depths.append(d)
        self.logger.debug("Y depths (first {} regions): {}"\
                    .format(N, np.array(depths)))
        return np.median(depths)


def read_alignment(samfile):
    ''' Dispatches BAM/CRAM based on file suffix
    '''
    tag = 'rc' if samfile.endswith(".cram") else 'rb'
    return pysam.AlignmentFile(samfile, tag)


def test_fetch(samfile, chr, start, end, logger):
    try:
        samfile.fetch(chr, start, end)
        return True
    except ValueError:
        logger.error("No reads extracted for region {}:{}-{}".format(chr, start, end))
        return False


def rc(s):
    cs = s.translate(_complement)
    return cs[::-1]
