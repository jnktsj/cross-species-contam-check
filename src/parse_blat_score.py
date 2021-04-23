#!/usr/bin/env python

import sys
import os.path
from argparse import ArgumentParser

def main(args):
    read_dict = dict()
    for line in open(args.query_fasta):
        if line.startswith(">"):
            read_dict.setdefault(line[1:].rstrip("\n"), dict())
    genomes = []
    for blat_output in args.blat_outputs:
        genome = os.path.basename(blat_output).split(".score")[0]
        genomes.append(genome)
        for line in open(blat_output):
            line = line.rstrip("\n").split("\t")
            read_name = ":".join(line[3].split(":")[:-1])
            score = int(line[4])
            read_dict[read_name][genome] = max([score, read_dict[read_name].get(genome, 0)])
    print("\t".join(["query"] + genomes))
    for r in read_dict:
        print("\t".join([r] + [str(read_dict[r].get(g, 0)) for g in genomes]))


if __name__ == "__main__":
    parser = ArgumentParser(description="Parse BLAT score output")
    parser.add_argument("query_fasta", help="query fasta")
    parser.add_argument("blat_outputs", nargs="+", help="BLAT score outputs")
    args = parser.parse_args()
    try:
        main(args)
    except KeyboardInterrupt: pass

