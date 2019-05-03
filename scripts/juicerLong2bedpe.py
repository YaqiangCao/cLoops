#!/usr/bin/env pyhthon2.7
#--coding:utf-8--
"""
"""

import argparse, gzip, os, sys
from datetime import datetime


def long2bedpe(fin, fout, ext=75):
    with open(fout, "w") as fo:
        for line in open(fin):
            line = line.split("\n")[0].split()
            nline = [
                line[1],
                max(0,
                    int(line[2]) - ext),
                int(line[2]) + ext,  #pet 1 
                line[5],
                max(0,
                    int(line[6]) - ext),
                int(line[6]) + ext,  #pet 2
                ".",
                ".",
                "+",
                "+"  #other infor
            ]
            if line[0] != "0":
                nline[-2] = "-"
            if line[4] != "0":
                nline[-1] = "-"
            fo.write("\t".join(list(map(str, nline))) + "\n")


def main():
    start_time = datetime.now()
    op = mainHelp()
    if not os.path.isfile(op.fin):
        sys.stderr.write("Error: input file %s not exists!\n" % op.fin)
    if os.path.isfile(op.fout):
        sys.stderr.write("Error: output file %s exists! \n" % op.fout)
    long2bedpe(op.fin, op.fout)
    usedtime = datetime.now() - start_time
    sys.stderr.write("Prcess finished. Used CPU time: %s Bye!\n" % usedtime)


def mainHelp():
    description = "Convert Juicer long format file to bedpe format for cLoops"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        '-i',
        '--input',
        dest="fin",
        type=str,
        help="Input file name, required.")
    parser.add_argument(
        '-o',
        '--out',
        dest="fout",
        required=True,
        type=str,
        help="Output file name, required.")
    op = parser.parse_args()
    return op


if __name__ == '__main__':
    main()
