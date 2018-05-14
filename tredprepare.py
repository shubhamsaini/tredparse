#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Example command to generate commands!! :D
cat large_loci.txt | awk '{print "python tredprepare.py --name "$1"_"$2" --motif "$5" --nummotifs "$4" --startposition "$1":"$2" --fasta ../fasta/human_g1k_v37.fasta"}'

"""

"""
Copyright (c) 2015-2017 Human Longevity Inc.
Author: Haibao Tang <htang@humanlongevity.com>
License: Non-Commercial Use Only. For details, see `LICENSE` file
Prepare details of a STR locus for computation with TREDPARSE.
The script progresses by asking a few key questions:
- Locus name
- Motif
- Expected number of motifs
- Chromosomal start location
"""
import argparse
import sys
import json

from pyfaidx import Fasta
from tredparse.jcvi import mkdir
from tredparse.meta import TREDsRepo


def survey_var(s, cast=str, default=""):
    """ Convenience function to get variable from command line survey.
    """
    res = raw_input("{} [{}]: ".format(s, default))
    if res == "":
        res = default
    return cast(res)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--name", help="Locus Name", required=True, type=str)
    parser.add_argument("--motif", help="Sequence Motif", required=True, type=str)
    parser.add_argument("--nummotifs", help="Number of motifs in reference", required=True, type=str)
    parser.add_argument("--startposition", help="Chromosomal start location (1-based)", required=True, type=str)
    parser.add_argument("--fasta", help="Fasta path", required=True, type=str)
    args = parser.parse_args()

    name = args.name
    repeat = args.motif
    nrepeat = int(args.nummotifs)
    location = args.startposition
    ref = args.fasta

    # Extract sequences
    f = Fasta(ref)
    c, start = location.split(":")
    start = int(start)
    size = len(repeat) * nrepeat
    end = start + size - 1
    tract = f[c][start - 1: end].seq.lower()
    FLANK = 18
    prefix = f[c][start - FLANK - 1: start - 1].seq.upper()
    suffix = f[c][end: end + FLANK].seq.upper()
    print >> sys.stderr, "|".join((prefix, tract, suffix))

    # Populate the tred
    repo = TREDsRepo()
    s = repo.create_tred()
    s["id"] = "{}_{}_{}".format(c, start, repeat)
    s["prefix"] = prefix
    s["suffix"] = suffix
    s["repeat"] = "CAG"
    s["repeat_location"] = "{}:{}-{}".format(c, start, end)
    tred = { name: s }

    # Serialize
    mkdir("sites")
    outfile = "sites/{}.json".format(name)
    fw = open(outfile, "w")
    print >> fw, json.dumps(tred, sort_keys=True, indent=2)
    fw.close()

    print >> sys.stderr
    print >> sys.stderr, "Template json file is written to `{}`".format(outfile)
    print >> sys.stderr, "Please manually fill in the remaining details"
    
    with open("sites/"+name+".json", "r+") as f:
        data = json.load(f)
        data[data.keys()[0]]['inheritance'] = "AD"
        data[data.keys()[0]]['motif'] = repeat
        data[data.keys()[0]]['repeat'] = repeat
        data[data.keys()[0]]['repeat_location.hg19'] = data[data.keys()[0]]['repeat_location']
    with open("sites/"+name+".json", 'w') as f:
        json.dump(data, f, indent=4)


if __name__ == '__main__':
    main()