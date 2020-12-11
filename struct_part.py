#!/usr/bin/env python

import os
import sys
import argparse
import itertools
from collections import defaultdict

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Given alignment in FASTA format with secondary structure 
                    annotation as a separate sequence, yield partitioning schemes 
                    in PartitionFinder format.
                    """
                    )
    parser.add_argument('input',
                        help='FASTA alignment with secondary structure as a sequence')
    parser.add_argument('-e', '--encoding', default='(){};-;*', 
                        help='characters used in annotation with different partitions '
                        'split by semicolon')
    parser.add_argument('-l', '--labels', default='stem;loop;excluded',
                        help='partition labels for the output file')
    parser.add_argument('-o', '--out_file', default=None, 
                        help='Output file name. Default: use input as prefix')

    return parser.parse_args()


def find_struct_seq(fname, charset):

    '''
    Given fasta alignment, return first sequence 
    completely matching the charset
    '''

    seqs = defaultdict(str)

    with open(fname) as f:
        for l in f:
            if l.startswith('>'):
                # header
                seq_name = l.strip()[1:]
            else:
                seqs[seq_name] += l.strip()

    for seq in seqs.values():
        if set(seq) - charset == set():
            return seq

def parse_struct(struct_str, partitions):

    charpos = defaultdict(list)

    # initial collection (0-based)
    for i, s in enumerate(struct_str):
        charpos[s].append(i)

    assert set(charpos.keys()) == set(''.join(partitions)), \
            'Unexpected error while defining partitions'

    # combine partitions
    struct_dict = dict()
    for partition in partitions:
        struct_dict[partition] = list()
        for character in partition:
            struct_dict[partition].extend(charpos[character])

    # convert to ranges
    def ranges(i):
        # https://stackoverflow.com/questions/4628333/
        # converting-a-list-of-integers-into-range-in-python
        for a, b in itertools.groupby(enumerate(i), 
                                      lambda pair: pair[1] - pair[0]):
            b = list(b)
            yield b[0][1], b[-1][1]

    struct_ranges = dict()
    for partition in partitions:
        positions = sorted(struct_dict[partition])
        struct_ranges[partition] = list(ranges(positions))

    return struct_ranges

def write_partitionfinder(struct_ranges, out_file, partitions, labels):
    '''
    Write ranges in PartitionFinder2 format,
    converting from 0-based to 1-based coordinates
    '''
    with open(out_file, 'w') as o:
        for partition, label in zip(partitions, labels):
            o.write('{} ='.format(label))
            for i, j in struct_ranges[partition]:
                # 0-based to 1-based coordinates
                o.write(' {}-{}'.format(i + 1,j + 1))
            o.write(';\n')

def main(args):
    
    partitions = args.encoding.split(';')
    labels = args.labels.split(';')
    assert len(partitions) == len(labels), \
        'Number of partitions and labels do not match'
    charset = set(''.join(partitions))

    struct_str = find_struct_seq(args.input, charset)

    struct_ranges = parse_struct(struct_str, partitions)

    if not args.out_file:
        prefix = '.'.join(args.input.split('.')[:-1])
        args.out_file = prefix + '.pf.txt'

    write_partitionfinder(struct_ranges, args.out_file, partitions, labels)
    print('Done!')
            
if __name__ == '__main__':
    main(parse_command_line_arguments())
    
