#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
2017-07-20: loops2juicebox added
2017-08-02: re-design the datastructure
"""

__author__ = "CAO Yaqiang"
__date__ = ""
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#sys
import os, random, gzip

#3rd
import joblib
import numpy as np

#cLoops
from cLoops.utils import callSys
from cLoops.utils import cFlush


class PET(object):
    #cA is the center of left read
    __slots__ = [
        "chromA", "chromB", "startA", "startB", "endA", "endB", "strandA",
        "strandB", "cA", "cB", "distance"
    ]

    def __init__(self, d):
        """
        d is line = line.split( "\n" )[ 0 ].split( "\t" ) from BEDPE file 
        """
        self.chromA = d[0]
        self.startA = int(d[1])
        self.endA = int(d[2])
        self.strandA = d[8]
        self.chromB = d[3]
        self.startB = int(d[4])
        self.endB = int(d[5])
        self.strandB = d[9]
        if self.chromA == self.chromB:
            #adjust the left end and right end to make sure left is alwasy small than right
            if self.startA + self.endA > self.startB + self.endB:
                self.startA, self.startB = self.startB, self.startA
                self.endA, self.endB = self.endB, self.endA
                self.strandA, self.strandB = self.strandB, self.strandA
            self.cA = (self.startA + self.endA) / 2
            self.cB = (self.startB + self.endB) / 2
            self.distance = self.cB - self.cA
        else:
            self.cA, self.cB, self.distance = None, None, None


def parseRawBedpe(fs, fout, cs, cut, logger):
    """
    Get the cis-PETs, organized by chromosomes. Input could be mixed PETs in bedpe.gz. Also change read id to numbers.
    @param fs: bedpe files of replicates, could be .bedpe or .bedpe.gz
    @param fout: output prefix, the name for directory
    @param cs: chroms that wanted, list like ["chr1","chr2"]
    """
    #chroms data
    chroms = {}
    #cis files
    cfs = []
    #distance between PETs mapped to different strands
    ds = []
    i, j, = 0, 0
    for f in fs:
        r = "Parsing PETs from %s, requiring initial distance cutoff > %s" % (
            f, cut)
        logger.info(r)
        if f.endswith(".gz"):
            of = gzip.open(f, "rb")
        else:
            of = open(f)
        for line in of:
            i += 1
            if i % 100000 == 0:
                cFlush("%s PETs processed from %s" % (i, f))
            line = line.split("\n")[0].split("\t")
            if "*" in line and "-1" in line:
                continue
            if len(line) < 6:
                continue
            try:
                pet = PET(line)
            except:
                continue
            #cis reads
            if pet.chromA != pet.chromB:
                continue
            #filtering unwanted PETs in chroms
            if len(cs) > 0 and (not (pet.chromA in cs)):
                continue
            #filtering too close PETs
            if cut > 0 and pet.distance < cut:
                continue
            if pet.chromA not in chroms:
                #f:file, c:count, r:redundant set
                #cf = os.path.join(fout,t[0]+".bedpe.gz")
                #chroms[t[0]] = {"f":gzip.open(cf,"wb"),"c":0,"r":set()}
                cf = os.path.join(fout,
                                  "%s-%s" % (pet.chromA, pet.chromB) + ".txt")
                chroms[pet.chromA] = {"f": open(cf, "w"), "c": 0, "r": set()}
                cfs.append(cf)
            if (pet.cA, pet.cB) in chroms[pet.chromA]["r"]:
                continue
            #in every file, [pointId,cA(x),cB(y)]
            nline = [chroms[pet.chromA]["c"], pet.cA, pet.cB]
            chroms[pet.chromA]["f"].write("\t".join(map(str, nline)) + "\n")
            chroms[pet.chromA]["c"] += 1
            chroms[pet.chromA]["r"].add((pet.cA, pet.cB))
            j += 1
            #collect distances of opsite strand PETs
            if pet.strandA != pet.strandB:
                ds.append(pet.distance)
    print
    del chroms
    r = "Totaly %s PETs from %s, in which %s cis PETs" % (i, ",".join(fs), j)
    logger.info(r)
    return cfs, ds


def txt2jd(f):
    """
    Dump the np.ndarray using joblib.dump for fast access.
    """
    data = []
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        data.append(map(int, line))
    data = np.array(data)
    joblib.dump(data, f.replace(".txt", ".jd"))
    os.remove(f)
    return f.replace(".txt", ".jd")


def parseJd(f, cut=0):
    """
    read data from .jd file 
    """
    key = os.path.split(f)[1].replace(".jd", "")
    key = tuple(key.split("-"))
    mat = joblib.load(f)
    if cut > 0:
        d = mat[:, 2] - mat[:, 1]
        p = np.where(d >= cut)[0]
        mat = mat[p, :]
    return key, mat


def loops2washU(fin, fout, logger, significant=1):
    """
    Convert interaction level loop file to washU long range interactions. 
    Track format according to http://wiki.wubrowse.org/Long-range
    @param fin: interactions in loop file
    @param fout: washU  long-range interaction text file prefix
    @param significant: if set 1, only convert significant loops.
    """
    logger.info("Converting %s to washU long range interaction track." % fin)
    with open(fout, "w") as f:
        for i, line in enumerate(open(fin)):
            if i == 0:
                continue
            line = line.split("\n")[0].split("\t")
            #only using significant results
            if significant and float(line[-1]) < 1:
                continue
            #iva,ivb,ES
            nline = [line[7], line[8], line[1]]
            f.write("\t".join(map(str, nline)) + "\n")
    logger.info(
        "Converting %s to washU long range interaction track finished." % fin)


def parseIv(iv):
    iv = [
        iv.split(":")[0],
        int(iv.split(":")[1].split("-")[0]),
        int(iv.split(":")[1].split("-")[1])
    ]
    return iv


def loops2juice(fin, fout, logger, significant=1):
    """
    Convert interaction level loop file to Juicebox 2D annotation features. 
    The txt file format according to https://github.com/theaidenlab/juicebox/wiki/Loading-Annotations-(Annotations-menu)
    @param fin: interactions in loop file
    @param fout: washU  long-range interaction text file prefix
    @param significant: if set 1, only convert significant loops.
    all p-values are -log10(p) transformed to escape all shown as 0 in juicebox.
    """
    logger.info("Converting %s to Juicebox 2D annotation feature." % fin)
    with open(fout, "w") as f:
        line = [
            "chromosome1", "x1", "x2", "chromosome2", "y1", "y2", "color",
            "observed", "loopId", "FDR", "EnrichmentScore", "distance",
            "-log10(binomal_p-value)", "-log10(poisson_p-value)",
            "-log10(hypergeometric_p-value)"
        ]
        f.write("\t".join(line) + "\n")
        for i, line in enumerate(open(fin)):
            if i == 0:
                continue
            line = line.split("\n")[0].split("\t")
            #only using significant results
            if significant and float(line[-1]) < 1:
                continue
            iva = parseIv(line[7])
            ivb = parseIv(line[8])
            nline = [
                iva[0], iva[1], iva[2], ivb[0], ivb[1], ivb[2], '"0,255,255"',
                line[11], line[0], line[2], line[1], line[4],
                -np.log10(float(line[3])), -np.log10(float(line[9])),
                -np.log10(float(line[6]))
            ]
            f.write("\t".join(map(str, nline)) + "\n")
    logger.info(
        "Converting %s to Juicebox 2D annotation feature finished." % fin)


def jd2washU(fs, fout, cut, ext):
    """
    Convert PETs to washU long range interactions. 
    Track format according to http://wiki.wubrowse.org/Long-range
    @param fs: files of .jd files 
    @param fout: prefix of output files
    """
    print "Converting %s to washU track." % (",".join(fs))
    tmp = str(random.random())
    with open(tmp, "w") as f:
        for fin in fs:
            print "converting %s" % fin
            key, mat = parseJd(fin, cut)
            for t in mat:
                a = (key[0], max([0, t[1] - ext]), t[1] + ext)
                b = (key[1], max([0, t[2] - ext]), t[2] + ext)
                linea = [
                    a[0], a[1], a[2],
                    "%s:%s-%s,1" % (b[0], b[1], b[2]), t[0], "."
                ]
                lineb = [
                    b[0], b[1], b[2],
                    "%s:%s-%s,1" % (a[0], a[1], a[2]), t[0], "."
                ]
                f.write("\t".join(map(str, linea)) + "\n")
                f.write("\t".join(map(str, lineb)) + "\n")
    c1 = "bedtools sort -i %s > %s" % (tmp, fout)
    c2 = "rm %s" % tmp
    c3 = "bgzip %s" % fout
    c4 = "tabix -p bed %s.gz" % fout
    callSys([c1, c2, c3, c4])
    print "Converting %s to washU random accessed track finished." % fout


def jd2hic(fs, fout, cut, org):
    """
    Convert reads level bedpe to HIC.
    Track format according to https://github.com/theaidenlab/juicer/wiki/Pre#file-format
    @param fs: files of .jd files 
    @param fout: prefix of output files
    """
    print "Converting %s to .hic file which could be loaded in juicebox" % (
        ",".join(fs))
    tmp = str(random.random())
    ss = {"+": 0, "-": 1}
    with open(tmp, "w") as f:
        for fin in fs:
            print "converting %s" % fin
            key, mat = parseJd(fin)
            for t in mat:
                line = [0, key[0], t[1], 0, 1, key[1], t[2], 1]
                f.write("\t".join(map(str, line)) + "\n")
    c1 = "juicer_tools pre -d {fin} {fout} {org}".format(
        fin=tmp, fout=fout, org=org)
    c2 = "rm %s" % tmp
    callSys([c1, c2])
    print "Converting %s to juicer's hic file finished." % fout
