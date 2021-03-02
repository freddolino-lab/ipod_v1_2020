#!/usr/bin/python

# write a gff file from the peaks in a gr file with peak flags

import gfftools
import sys
from peaktools import *
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('gr_file', help="gr file containing the flag locations")
parser.add_argument('output_file', help="name of output .gff file")
parser.add_argument('--chrname', help="name of chromosome to use in output", default="E_coli_genome")
parser.add_argument('--origin', help="name to use in the data origin field", default="E_coli_genome")
parser.add_argument('--sitetype', help="name to use in the site type field", default="Peak")

args=parser.parse_args()


ingr=args.gr_file
outgff=args.output_file

mylib=gfftools.GffData()
peaklist = read_flaggrfile(ingr)

for i in range(peaklist.numpeaks):
  startloc=peaklist.peakstarts[i]
  endloc = peaklist.peakends[i]
  new_gff_entry = gfftools.GffEntry()
  new_gff_entry.genome_name = args.chrname
  new_gff_entry.data_origin = args.origin
  new_gff_entry.site_type = args.sitetype
  new_gff_entry.comments = 'From %s' % ingr
  # the 1s are added since gff files are 1-indexed
  new_gff_entry.start = startloc + 1
  new_gff_entry.end = endloc + 1
  mylib.add_entry(new_gff_entry)

mylib.write_gff_file(outgff)





