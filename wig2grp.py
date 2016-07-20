#!/usr/bin/env python


# wig2grp.py, convert wiggle format to artemis graphic user plot file (grp)
# Copyright (C) 2016  Shengwei Hou
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import os
import sys
import argparse
import numpy as np
from Bio import SeqIO


def get_genome_size(input_file):
    """ parse fasta records in input_file, 
        return the fastaName:genomeSize as a dictionary
    """
    genomeSizeDict = {}
    for fasta in SeqIO.parse(input_file, "fasta"):
        name = fasta.name.strip()
        length = len(fasta.seq)
        genomeSizeDict[name] = length
    return genomeSizeDict


def wig2grp(genomeSizeDict, wig):
    """ convert wig to grp file
    """
    chromArrayDict = {} # {chrom:Array}
    for chrom in genomeSizeDict.keys():
        chromArrayDict[chrom] = np.zeros(genomeSizeDict[chrom])

    with open(wig, "r") as ih:
        for line in ih:
            if line.startswith("variableStep"):
                chrom = line.strip().split("chrom=")[-1].strip()
                if not chrom in genomeSizeDict:
                    print >> sys.stderr, "ERROR: chromosome %s in your wig file was not found in your genome reference! "%chrom
                    sys.exit(1)
            else:
                line = line.strip().split("\t")
                index = int(line[0])-1
                value = float(line[1])
                if index >= len(chromArrayDict[chrom]):
                    print >> sys.stderr, "ATTENTION: Mapped position at %d was out of genome length %d, on chromosome %s"%(index, genomeSizeDict[chrom], chrom)
                    continue
                chromArrayDict[chrom][index] = value

    return chromArrayDict                


def write2grp(genomeSizeDict, fwd=None, rev=None, prefix=None, toInteger=True):
    """ convert wig to grp file
    """

    # only fwd or rev wig
    if (fwd and (not rev)) or (rev and (not fwd)):
        
        if fwd:
            header = "# BASE fwd_grp\n# colour 0:255:0\n"
            wig = fwd
        else:
            header = "# BASE rev_grp\n# colour 0:0:255\n"
            wig = rev

        chromArrayDict = wig2grp(genomeSizeDict, wig)
        for chrom, value_array in chromArrayDict.iteritems():
            pos_array = np.arange(genomeSizeDict[chrom]) + 1
            if prefix:
                outfile = prefix +"_"+ chrom + ".grp"
            else:
                outfile = os.path.basename(wig) +"_"+ chrom +".grp"
            with open(os.path.join(os.getcwd(), outfile), "w") as oh:
                oh.write(header)
                if toInteger:
                    for pos, value in zip(pos_array, value_array):
                        oh.write("%d %d\n"%(pos, int(value)))
                else:
                    for pos, value in zip(pos_array, value_array):
                        oh.write("%d %.2f\n"%(pos, value))
                        
    # both fwd and rev wig, integrate into one grp
    else:
        header = "# BASE fwd_grp rev_grp\n# colour 0:255:0 0:0:255\n"
        fwd_chromArrayDict = wig2grp(genomeSizeDict, fwd)
        rev_chromArrayDict = wig2grp(genomeSizeDict, rev)
        for chrom, fwd_value_array in fwd_chromArrayDict.iteritems():
            pos_array = np.arange(genomeSizeDict[chrom]) + 1
            rev_value_array = rev_chromArrayDict[chrom]
            if prefix:
                outfile = prefix +"_"+ chrom + ".grp"
            else:
                if "Forward" in fwd:
                    outfile = os.path.basename(fwd).replace("Forward", "FwdRev") +"_"+ chrom +".grp"
                else:
                    outfile = os.path.basename(fwd).rsplit(".", 1)[0] +"_FwdRev_"+ chrom +".grp"
            with open(os.path.join(os.getcwd(), outfile), "w") as oh:
                oh.write(header)
                if toInteger:
                    for pos, fwd_value, rev_value in zip(pos_array, fwd_value_array, rev_value_array):
                        oh.write("%d %d %d\n"%(pos, int(fwd_value), int(rev_value)))
                else:
                    for pos, fwd_value, rev_value in zip(pos_array, value_array, rev_value_array):
                        oh.write("%d %.2f %.2f\n"%(pos, fwd_value, rev_value))


def main():
    desc="Convert variableStep forward and reverse wig files into fixedStep grp file, to visualize in Artemis"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-g", "--genomeFasta", required=True, help="genome fasta file that used in the mapping")
    parser.add_argument("-f", "--forward", required=False, help="forward wig file generated by RSeQC::bam2wig.py")
    parser.add_argument("-r", "--reverse", required=False, help="reverse wig file generated by RSeQC::bam2wig.py")
    parser.add_argument("-p", "--prefix", required=False, help="output prefix for grp file")
    parser.add_argument("-t", "--toInteger", required=False, action="store_true", default=True, help="write out integer values")
    args = parser.parse_args()

    genomeSizeDict = get_genome_size(args.genomeFasta)

    if not (args.forward or args.reverse):
        print >> sys.stderr, "\nERROR: At least one wig file need to be given!\n"
        parser.print_help()
        sys.exit(1)
    else:
        write2grp(genomeSizeDict, fwd=args.forward, rev=args.reverse, prefix=args.prefix, toInteger=args.toInteger)

    print "---- DONE ----"


if __name__ == "__main__":
    main()
