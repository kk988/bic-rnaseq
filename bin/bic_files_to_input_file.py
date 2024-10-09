#!/usr/bin/env python3

# This code will take in BIC Mapping and make an input file for the NF pipeline
# Requirements: mapping file, access to the fastq (as fastq.gz).

# This will not make a perfect input file. The sample names in nf-rnaseq double as
# the grouping file as well.

# the sample names should only have one underscore in them

#Ex:
# samplename,fastq1,fastq2,strandedness
# ABCD_1,L1_fq1,L1_fq2,auto
# ABCD_1,L2_fq1,L2_fq2,auto
# ABCD_2,fq1,fq2,forward
# ABCD_3,fq1,fq2,forward
# ZYXX_1,fq1,fq2,reverse
# ZYXX_2,fq1,fq2,reverse
# ZYXX_3,fq1,fq2,unstranded

import argparse
import glob
import re

lane_pattern = r"_S.*_(L\d{3})_(R[12])_\d{3}.f"

def parse_args():
    parser = argparse.ArgumentParser(description='Convert BIC mapping and pairing files to input file for NF pipeline')
    parser.add_argument('-m', '--mapping', type=str, required=True, help='BIC mapping file')
    parser.add_argument('-s', '--strandedness', type=str, choices=['auto', 'forward', 'reverse', 'unstranded'], required=True, help='strandedness of the fastq')
    parser.add_argument('-k', '--key', type=str, required=False, help='key file - use if some of the samples are not being compared in DE analysis', default=None)
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file')
    parser.add_argument('--fq_pattern', type=str, required=False, help='pattern to find fq in mapping folder (default: /*_L*_R1_*.fastq.gz)', default='/*_L*_R1_*.fastq.gz')
    return parser.parse_args()

def generate_input_file(args):
    mapping = generate_mapping_dict(args.mapping, args.key)
    write_input_file(mapping, args.strandedness, args.output, args.fq_pattern)


def write_input_file(mapping, strandedness, output, fq_pattern):

    outfile = open(output, 'w')
    print('sample,fastq_1,fastq_2,strandedness', file=outfile)
    for sample in mapping:
        print_input_lines(sample, mapping, strandedness, outfile, fq_pattern)

def print_input_lines(sample, mapping, strandedness, outfile, fq_pattern):
    for fq_path in mapping[sample]:
        fq_r1s = glob.glob(fq_path + fq_pattern)
        if not fq_r1s:
            print("No fastq files found for sample: {} in path {}. Will be skipped.".format(sample, fq_path))
            return
        for fq_r1 in fq_r1s:
            fq_r2 = re.sub(r'_R1_', '_R2_', fq_r1)
            if not glob.glob(fq_r2):
                fq_r2 = ""
            print(sample, fq_r1, fq_r2, strandedness, sep=',', file=outfile)

def generate_mapping_dict(mapping, key=None):
    mapping_dict = {}
    with open(mapping, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[1] not in mapping_dict:
                mapping_dict[line[1]] = [line[3]]
            else:
                mapping_dict[line[1]].append(line[3])

    if key:
        with open(key, 'r') as k:
            for line in k:
                line = line.strip().split('\t')
                if line[1] == "_EXCLUDE_":
                    del mapping_dict[line[0]]

    return mapping_dict

if __name__ == '__main__':
    args = parse_args()
    generate_input_file(args)
