import sys
import sam_utils
import argparse
import random
import logging
import numpy as np
from math import floor


## TO DO:
# add support for multiple chromosomed organisms
# add support for strandedness

class ReadSampler(object):
    def __init__(self):
        self.reads = []
        self.total = 0
        self.sampling = False

    def add_read(self, new_read):
        if self.sampling:
            self.convert_to_list()
        self.reads.append(new_read)
        self.total+=1

    def add_reads(self, new_reads):
        if self.sampling:
            self.convert_to_list()
        self.reads.extend(new_reads)
    
    def convert_to_array(self):
        self.reads = np.asarray(self.reads, dtype="int64")
        self.sampling=True

    def convert_to_list(self):
        self.reads = list(self.reads)
        self.sampling = False

    def pull_read(self):
        if not self.sampling:
            self.convert_to_array()
        index = np.random.randint(0, self.total)
        return self.reads[index, :]

    def pull_reads(self, n):
        if not self.sampling:
            self.convert_to_array()
        index = np.random.randint(0, self.total, size=n)
        index = np.sort(index)
        return self.reads[index,:]

    def save_data(self, f):
        if not self.sampling:
            self.convert_to_array()
        np.save(f, self.reads)
    
    def sort_reads(self):
        if not self.sampling:
            self.convert_to_array()
        self.reads = self.reads[self.reads[:,0].argsort()]
    def load_data(self, f):
        self.sampling = True
        self.reads = np.load(f)
        self.total = self.reads.shape[0]



def merge(intervals):
    intervals.sort(key=lambda x: x[0])
    # take the first interval
    merged = [intervals[0]]
    # loop through all the intervals
    for this_interval in intervals:
        if this_interval[0] <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], this_interval[1]))
        else:
            merged.append(this_interval)
    return merged

def get_paired_blocks(r1, r2):
    if r1.TLEN > 0 and r2.TLEN < 0:
        left = r1.get_aligned_blocks()
        right = r2.get_aligned_blocks()
    elif r1.TLEN < 0 and r2.TLEN > 0:
        left = r2.get_aligned_blocks()
        right = r1.get_aligned_blocks()
    elif r1.POS == r2.POS:
        left = r1.get_aligned_blocks()
        right = r2.get_aligned_blocks()
    else:
        raise ValueError("Pair not consistent %s %s"%(r1.QNAME, r2.QNAME))
    total_blocks = []
    if left[-1][1] < right[0][0]:
        total_blocks.append((left[-1][1], right[0][0]))
    total_blocks.extend(left)
    total_blocks.extend(right)
    total_blocks = merge(total_blocks)
    if len(total_blocks) > 1:
        raise RuntimeError("Gapped read found %s %s %s"%(r1.QNAME, r2.QNAME, str(total_blocks)))
    return total_blocks[0]

def create_read_list(samfile):
    read_sampler = ReadSampler()
    for line in samfile:
        line = sam_utils.SamAlignment(line)
        vals = line.get_aligned_blocks()
        if len(vals) > 1:
            logging.info("Skipping gapped read %s %s"%(line.QNAME, str(vals)))     
        read_sampler.add_read(vals[0])
    return read_sampler

def create_read_list_paired(samfile):
    read_sampler = ReadSampler()
    while True: 
        line1 = samfile.readline()
        line2 = samfile.readline()
        if not line2: 
            break
        line1 = sam_utils.SamAlignment(line1)
        line2 = sam_utils.SamAlignment(line2)
        if line1.QNAME != line2.QNAME:
            raise ValueError("Unpaired read or read with more than one pair\
                              encountered. Check your input file. File must\
                              be sorted by read name, every read must have\
                              a single pair and each pair must have one\
                              mapping. %s %s"%(line1.QNAME, line2.QNAME))
        try:
            read_sampler.add_read(get_paired_blocks(line1,line2))
        except ValueError as err:
            logging.error("Skipping pair %s"%err)
        except RuntimeError as err:
            logging.error("Skipping pair %s"%err)
    return read_sampler

def map_read(array, read):
    start, stop = read
    #array[start:stop] += 1
    # the below line implements linear scaling with read length
    array[start:stop] += 100.0/(stop-start)

def sample(read_sampler, n, array):
    for read in read_sampler.pull_reads(n):
        map_read(array, read)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='commands', dest='command')

    # parse
    parse_parser = subparsers.add_parser('parse', help="create a sampling\
            object from a sam file")
    parse_parser.add_argument('samfile', help="Input samfile, this tool does no\
            filtering and will consider every line in the file. Accepts input\
            from stdin if '-' is specified here.")
    parse_parser.add_argument('outpre', help="output prefix to np.save the\
            sampling object data that is created")
    parse_parser.add_argument('--paired', action="store_true", help="Consider\
            the sam file as paired. If this flag is specified then the sam file\
            MUST be pre-filtered to have only ONE alignment per pair. Further,\
            there must be NO unpaired reads in the file and the reads must be\
            sorted by read name.")

    # sample
    sample_parser = subparsers.add_parser('sample', help="sample coverage from a\
            sampling object.")
    sample_parser.add_argument('samplerfile', help="output file from using parse\
            on the samfile of interest")
    sample_parser.add_argument('outpre', help="output file to np.save the\
            numpy array that is created")
    sample_parser.add_argument('array_size',type=int, help="length of genome")
    sample_parser.add_argument('--num_samples', type=int, default=1,
    help="number of full samples to pull from the sampler, default is 1")
    sample_parser.add_argument('--num_reads', type=int, default=None, 
    help="number of reads to pull for each sample. Default is the size of\
            sampling object.")
    sample_parser.add_argument('--identity', action="store_true",
            help="write an array of the actual sample without sampling, ignores\
                    other optional arguments")
    sample_parser.add_argument('--resolution', type=int, default=1,
            help="only keep data for one bp out of this number")
    args = parser.parse_args()

    if args.command == "parse":
        if args.samfile == "-":
            f = sys.stdin
        else:
            f = open(args.samfile, mode="r")
        if args.paired:
            sampler = create_read_list_paired(f)
        else:
            sampler = create_read_list(f)
        f.close()
        sampler.sort_reads()
        sampler.save_data(args.outpre)
    elif args.command == "sample":
        array = np.zeros((args.array_size, args.num_samples))
        sampler = ReadSampler()
        sampler.load_data(args.samplerfile)
        if args.identity:
            for read in sampler.reads:
                map_read(array, read)
            np.save(args.outpre, array[::args.resolution,])
        else:
            for i in xrange(args.num_samples):
                if args.num_reads:
                    num_reads = args.num_reads
                else:
                    num_reads = sampler.total 
                sample(sampler, num_reads, array[:,i])
            np.save(args.outpre, array[::args.resolution,])
