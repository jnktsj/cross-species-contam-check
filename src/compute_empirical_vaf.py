#!/usr/bin/env python3

import sys
import os.path
from argparse import ArgumentParser

import pysam
import pandas as pd


def main(args):
    read_chrom = False
    index = []
    bam = pysam.AlignmentFile(args.bam, "rb")
    columns = ['Chromosome', 'Start_position', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele2']
    for line in open(args.maf):
        line = line.rstrip("\n").split("\t")
        if line[0].startswith("Hugo"):
            index = [line.index(c) for c in columns]
            continue

        if read_chrom and args.chrom != rec[0]:
            break

        rec = [line[i] for i in index]

        if rec[2] != "SNP":
            continue

        if args.chrom == rec[0]:
            read_chrom = True
            t_alt_count = 0
            t_ref_count = 0
            for col in bam.pileup(rec[0], int(rec[1])-1, int(rec[1]), truncate=True, stepper="all", min_base_quality=5, min_mapping_quality=20):
                for r in col.pileups:
                    if r.is_del or r.is_refskip:
                        continue
                    read_name = r.alignment.query_name
                    read_seq = r.alignment.query_sequence
                    # ref base
                    if read_seq[r.query_position] == rec[3]:
                        t_ref_count += 1
                    # alt base
                    if read_seq[r.query_position] == rec[4]:
                        t_alt_count += 1
            rec.pop(2)
            rec.append(str(t_alt_count/float(t_alt_count+t_ref_count)))
            print("\t".join(rec))



if __name__ == "__main__":
    description = "Compute empirical VAF from CGA MAFs"
    parser = ArgumentParser(description=description)
    parser.add_argument("chrom", help="target chromosome")
    parser.add_argument("maf", help="input MAF")
    parser.add_argument("bam", help="BAM")
    args = parser.parse_args()
    try:
        main(args)
    except KeyboardInterrupt: pass
