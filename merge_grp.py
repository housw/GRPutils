#!/usr/bin/env python


import sys, os
import numpy as np
import argparse


def read_grp(path2grp):
    """ This function used to read in grp file, input file format like follows

        # # BASE    fwd_grp rev_grp
        # # colour  0:204:0 0:0:255
        #
        0.00    0.00   -0.00
        1.00    0.00   -0.00
        2.00    0.00   -0.00
        3.00    0.00   -0.00

        and will return fwd_grp and rev_grp numpy array instances
    """
    fwd_grp = []
    rev_grp = []

    with open(path2grp, "r") as ih:
        for line in ih:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            fwd_grp.append(float(line[1]))
            rev_grp.append(abs(float(line[2])))
    fwd_grp = np.array(fwd_grp)
    rev_grp = np.array(rev_grp)

    return fwd_grp, rev_grp


def average_grpFile(grpFile_list, prefix, toInteger):
    """
    """
    # number of grp files
    n = len(grpFile_list)

    # read in first one
    fwd_grp, rev_grp = read_grp(grpFile_list[0])

    for grp_file in grpFile_list[1:]:
        tmp_fwd_grp, tmp_rev_grp = read_grp(grp_file)
        fwd_grp += tmp_fwd_grp
        rev_grp += tmp_rev_grp

    # initialize array as long as genome length
    pos_array = np.arange(len(fwd_grp))
    pos_array = pos_array+1

    # take average, make rev to minus
    fwd_grp = fwd_grp / float(n)
    rev_grp = -1*rev_grp / float(n)

    if prefix:
        outfile = prefix+".grp"
    else:
        basename = os.path.basename(grpFile_list[0])
        filestem = os.path.splitext(basename)[0]
        outfile = filestem+"_averaged.grp"

    with open(outfile, "w") as oh:
        # to customize the line color, we must use space delimited format
        oh.write("# BASE fwd_grp rev_grp\n")
        oh.write("# colour 0:204:0 0:0:255\n")
        if toInteger:
            for pos, fwd, rev in zip(pos_array, fwd_grp, rev_grp):
                oh.write(str(pos)+" %d %d\n"%(int(fwd), int(rev)))
        else:
            for pos, fwd, rev in zip(pos_array, fwd_grp, rev_grp):
                oh.write(str(pos)+" %.2f %.2f\n"%(fwd, rev))


def main():
    mode_help="Average: sum the value, then take avarage of the sum. \nAppend: append multiple grp columns from each grp file"
    parser = argparse.ArgumentParser(description="Merge individual grp files into one grp file")
    parser.add_argument("-i", "--inputFile", required=True, nargs="+", help="grp files to be merged")
    parser.add_argument("-m", "--mode", required=False, choices=["average", "append"], default="average", help=mode_help)
    parser.add_argument("-p", "--prefix", required=False, help="output prefix for grp file")
    parser.add_argument("-t", "--toInteger", required=False, action="store_true", default=True, help="write out integer values")
    args = parser.parse_args()
    print("[merge_grp] input grp files are: ", ";".join(args.inputFile))

    average_grpFile(grpFile_list=args.inputFile, prefix=args.prefix, toInteger=args.toInteger)



if __name__ =="__main__":
    main()
