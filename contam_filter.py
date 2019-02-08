#!/usr/bin/env python

import os
import sys
import pysam
import argparse

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Given BAM alignments of same reads to target and 
                    contamination genome, remove contaminant reads and 
                    perform quality filtering. Contamination is detected by 
                    comparing mean mapping qualities of all aligned 
                    segments for each read.
                    """
                    )
    parser.add_argument('target_bam',
                        help='BAM alignment to target genome')
    parser.add_argument('contam_bam',
                        help='BAM alignment to contamination genome')
    parser.add_argument('-o', '--out_bam', default='test.bam', 
                        help='Output BAM file name')
    parser.add_argument('-s', '--stat_file', default='test.log', 
                        help='Statistics file name')
    parser.add_argument('-q', '--min_quality', 
                        type=int, default=20,
                        help='Minimum quality for filtered fragment. '
                             'Default: 20.')
    parser.add_argument('-l', '--min_length', 
                        type=int, default=20,
                        help='Minimum alignment length for each aligned '
                             'fragment. Default: 20 bp.')
    parser.add_argument('-u', '--keep_unmapped', 
                        action='store_true', default=False,
                        help='Keep unmapped reads, min_quality and min_length '
                             ' be igrored. Note that even with "-q 0 -l 0" '
                             'unmapped reads will be discarded')
    parser.add_argument('-d', '--dry_run', 
                        action='store_true', default=False, 
                        help='Print out commands and exit')

    return parser.parse_args()


def bam_sort(in_file, out_file, name_sort=False, verbose=True, dry_run=False):
    
    if verbose:
        sys.stderr.write('Sorting %s ' % (in_file))
    if name_sort:
        if verbose:
            sys.stderr.write('by read name. Output: %s.\n' % (out_file))
        if not dry_run:
            pysam.sort('-n', '-T', '/tmp/bam_nsort', '-o', out_file, in_file)
    else:
        if verbose:
            sys.stderr.write('by reference coordinate. Output: %s.\n' 
                             % (out_file))
        if not dry_run:
            pysam.sort('-T', '/tmp/bam_sort', '-o', out_file, in_file)


def filter_alns(t_alns, c_alns, f_file,
                q_count, r_count, s_count,
                min_qual=20, min_len=20, 
                keep_unmapped=False):
    """Given pysam alignedSegment lists t_alns and c_alns, from t_alns 
    remove fragments with higher mean mapping quality alignments in c_alns
    with matching read name (all t_alns are assumed to have one read name), 
    apply filters or keep all sequences, 
    write filtered alignments to open f_file, 
    return updated tuple of:
        q_count - number of alignments removed due to quality filter
        r_count - number of alignments removed as contamination

    Note that read might be represented by multiple aligned fragments in case
    of BWA MEM alignment and mapq filter is passed if mean mapping quality 
    among aligned fragments is higher or equal in target compared to 
    contamination alignment.
    """
    
    def mean_mapq(alns):
        """Return mean mapping quality for input alignment list
        """
        sum_mapq = sum([x.mapping_quality for x in alns])
        return sum_mapq / float(len(alns))

    def fragment_pass_filters(aln, min_qual, min_len, keep_unmapped=False):
        """Return True if read passes filters

        min_qual filter is passed if fragment mappinq quality is higher or equal
        than threshold
        
        min_len filter is passed if fragment alignment length on reference is
        higher or equal than threshold

        keep_unmapped is always passed if set to True
        """
        if aln.mapping_quality >= min_qual and aln.reference_length >= min_len:
            return True
        elif keep_unmapped:
            return True
        else:
            return False
    
    # count and remove non-matching contaminant reads     
    read_name = t_alns[0].query_name
    match_c_alns = [c_aln for c_aln in c_alns if c_aln.query_name == read_name]
    # calculate mean mapping qualities
    t_mean_mapq = mean_mapq(t_alns)
    if len(match_c_alns) == 0:
        c_mean_mapq = -1
    else:        
        c_mean_mapq = mean_mapq(c_alns)
    # compare mapqs and apply filters
    if t_mean_mapq >= c_mean_mapq:
        for t_aln in t_alns:
            if fragment_pass_filters(t_aln, min_qual, min_len, 
                                     keep_unmapped):
                f_file.write(t_aln)
                s_count += 1
            else:
                q_count += 1
    else:
        r_count += len(t_alns)

    return (q_count, r_count, s_count)

def higher_read_name(r1, r2):
    '''
    Compare read names: split at ':' and compare each field as int, if possible.
    Returns True is r1 name is higher,
    '''
    r1 = r1.split(':')
    r2 = r2.split(':')
    assert len(r1) == len(r2)
    for i in range(len(r1)):
        try:
            if int(r1[i]) < int(r2[i]):
                return False
            if int(r1[i]) > int(r2[i]):
                return True
        except ValueError:
            if r1[i] < r2[i]:
                return False
            elif r1[i] > r2[i]:
                return True
    else: # equity
        return False

def generate_filenames(t, c, o):

    return {'t': t,
            'c': c,
            'o': o,
            'tns': t + '.namesort.temp',
            'cns': c + '.namesort.temp',
            'ons': o + '.namesort.temp'
            }

def main(args):
    
    for f in (args.target_bam, args.contam_bam, args.out_bam):
        assert f.endswith('.bam')
    fn = generate_filenames(args.target_bam, args.contam_bam, args.out_bam)
   
    bam_sort(fn['t'], fn['tns'], name_sort=True, verbose=False, 
                   dry_run=args.dry_run)
    bam_sort(fn['c'], fn['cns'], name_sort=True, verbose=False, 
                   dry_run=args.dry_run)

    sys.stderr.write('Filtering %s using contaminant alignment %s\n' 
                     % (args.target_bam, args.contam_bam))
    if args.keep_unmapped:
        sys.stderr.write('Unmapped reads retained in the resulting file\n')
    else:
        sys.stderr.write('Filtering options: min_qual = %d, min_len = %d\n'
                         % (args.min_quality, args.min_length))
    sys.stderr.write('Output: %s\n' % args.out_bam)

    if not args.dry_run:
        
        with pysam.AlignmentFile(fn['tns']) as t_file, \
             pysam.AlignmentFile(fn['cns']) as c_file, \
             pysam.AlignmentFile(fn['ons'], 'wb', template=t_file) as f_file:
            
            t_data = t_file.fetch(until_eof=True)
            c_data = c_file.fetch(until_eof=True)
            t_alns = [next(t_data)]
            c_alns = [next(c_data)]
            read_name = t_alns[0].query_name
            (t_count, c_count, q_count, r_count, s_count) = (1, 1, 0, 0, 0)

            for t_aln in t_data:
                t_count += 1
                if t_aln.query_name == read_name: # accumulate targets
                    t_alns.append(t_aln)
                else: # get contaminant reads and compare mapping qualities
                    # previous contam name is higher than analyzed target - is target
                    if higher_read_name(c_alns[0].query_name, read_name): 
                        (q_count, r_count, s_count) = filter_alns(t_alns, c_alns, f_file,
                                q_count, r_count, s_count,
                                min_qual=args.min_quality, 
                                min_len=args.min_length,
                                keep_unmapped=args.keep_unmapped)
                        t_alns = [t_aln]
                        # do not reset previous contam
                        read_name = t_aln.query_name
                        continue
                    # aggregate further contam while name is not higher than in target
                    for c_aln in c_data:
                        c_count += 1
                        if not higher_read_name(c_aln.query_name, read_name): # reads not in target alignment are contaminant
                            c_alns.append(c_aln)
                        else:
                            break
                    # compare reads
                    (q_count, r_count, s_count) = filter_alns(t_alns, c_alns, f_file,
                                            q_count, r_count, s_count,
                                            min_qual=args.min_quality, 
                                            min_len=args.min_length,
                                            keep_unmapped=args.keep_unmapped)
                    # take next pair of target and contam alignments
                    t_alns = [t_aln]
                    c_alns = [c_aln]
                    read_name = t_aln.query_name
            # read the rest of contam after last target
            for c_aln in c_data: 
                c_alns.append(c_aln)
            # compare reads
            (q_count, r_count, s_count) = filter_alns(t_alns, c_alns, f_file,
                                    q_count, r_count, s_count,
                                    min_qual=args.min_quality, 
                                    min_len=args.min_length,
                                    keep_unmapped=args.keep_unmapped)

            # check that read numbers are correct
            assert t_count - q_count - r_count == s_count
            # write log
            with open(args.stat_file, 'w') as log:
                log.write('%d segments in target alignment\n' % (t_count))
                log.write('%d segments in contamination '
                          'alignment\n' % (c_count))
                log.write('%d segments removed as belonging to contaminant '
                              'genome\n' % (r_count))
                if args.keep_unmapped:
                    log.write('Unmapped segments retained in the output\n')
                else:
                    log.write('%d segments removed as not passing filters '
                              '(minimum MAPQ %d, minimum alignment length %d,'
                              ' discard unmapped)\n' % (q_count,
                                                          args.min_quality,
                                                          args.min_length))
                log.write('%d segments retained after contamination removal '
                          'and filtering\n' % (t_count - q_count - r_count))
            
                
    # coordinate-sort output
    bam_sort(fn['ons'], fn['o'], verbose=False, dry_run=args.dry_run)
    
    # clean-up
    for f in (fn['tns'], fn['cns'], fn['ons']):
        if os.path.isfile(f):
            os.unlink(f)
        
if __name__ == '__main__':
    main(parse_command_line_arguments())
    
