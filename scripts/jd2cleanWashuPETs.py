#!/usr/bin/env python
#--coding:utf-8 --
"""
jd2cleanPETs
"""
__date__ = "2017-10-10"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import os,argparse,bisect
from glob import glob
from copy import deepcopy
from datetime import datetime

#3rd library
import numpy as np
import pandas as pd
import joblib
from joblib import Parallel, delayed

#cLoops
from cLoops.utils import getLogger
from cLoops.io import parseJd,parseIv
from cLoops.cModel import getGenomeCoverage

#global settings
global logger
NAME=os.path.basename(__file__)
#epilog for argparse
EPILOG = "Any bug is welcome reported to caoyaqiang0410@gmail.com, chenzhaoxiong@picb.ac.cn"

def jd2cleanWashuPETsHelp():
    """
    Create the command line interface for the script of jd2CleanPETs.
    """
    description = """
        Filtering raw PET to keep only PET located in loop anchors.
        For exmaple:
        %s -d trac -f trac.loop -o trac_clean 
        """%NAME
    parser = argparse.ArgumentParser(description=description, epilog=EPILOG)
    parser.add_argument(
        "-d",
        dest="d",
        required=True,
        type=str,
        help=
        "The directory of cis .jd file."
    )
    parser.add_argument("-f",
                        dest="f",
                        required=True,
                        type=str,
                        help="Loops file called by cLoops.")
    parser.add_argument(
        "-p",
        dest="cpu",
        required=False,
        default=1,
        type=int,
        help=
        "CPU number used to run the job, default is 1,set -1 to use all cpus available. Too many CPU could cause memory error."
    )
    parser.add_argument(
        "-s",
        dest="significant",
        required=False,
        action="store_false",
        help=
        "Whether to only using the PETs located at loop anchors. Default is yes, set this flag to use all potential anchors called in the loop file."
    )
    parser.add_argument(
        "-c",
        dest="chroms",
        required=False,
        default="",
        type=str,
        help=
        "Whether to process limited chroms, specify it as chr1,chr2,chr3, default is processed all chroms."
    )
    parser.add_argument("-o",
                        dest="output",
                        required=True,
                        type=str,
                        help="Output prefix.")
 
    op = parser.parse_args()
    return op


def preDs(f, d, sig=True, chroms=[], ivac=6, ivbc=7):
    """
    Prepare input datasets into pairs, f is the .loop file and d is the directory contains parsed cis-PETs.
    """
    records = {}  #organized by chromosomes.
    if len(chroms) > 0:
        for c in chroms:
            records[c] = {"rs": {}, "f": ""}
    for i, line in enumerate(open(f)):
        if i == 0:
            continue
        line = line.split("\n")[0].split("\t")
        #only using significant loops
        if sig and float(line[-1]) < 1:
            continue
        iva = parseIv(line[ivac])
        ivb = parseIv(line[ivbc])
        if len(chroms) > 0 and iva[0] not in chroms:
            continue
        if iva[0] not in records:
            records[iva[0]] = {"rs": {}, "f": ""}
        records[iva[0]]["rs"][line[0]] = iva + ivb
    for chrom in records.keys():
        if len(records[chrom]["rs"]) == 0:
            del records[chrom]
            continue
        f = os.path.join(d, "%s-%s.jd" % (chrom, chrom))
        if os.path.isfile(f):
            records[chrom]["f"] = f
        else:
            logger.warning(
                "%s not found, however there are loops in that chromosome." %
                f)
            del records[chrom]
    return records


def getCorLink(cs):
    """
    @param cs: [1,2,3,4], a list for the coordinates x or y
    """
    ts = {}
    for i, c in enumerate(cs):
        ts.setdefault(c, []).append(i)
    #ts_keys = sorted(cs)
    ts_keys = np.sort(cs)
    return ts_keys, ts


def checkAnchorOverlap(xa, xb, ya, yb):
    """
    check the overlap of a region for the same chromosome
    """
    if ya <= xa <= yb or ya <= xb <= yb:
        return True
    if xa <= ya <= xb or xa <= yb <= xb:
        return True
    return False

def mergeAnchor(xa, xb,ya,yb):
    """
    merge the overlaped two anchors
    """
    nr = [
        min([xa, ya]),
        max([xb, yb]), 
    ]
    return nr


def mergeAllAnchors(rs):
    """
    Meget all anchors by search overlaps.
    """
    nrs = []
    skips = set()
    for i in range(len(rs)):
        if i in skips:
            continue
        nr = deepcopy(rs[i])
        for j in range(i + 1, len(rs)):
            if j in skips:
                continue
            nrj = rs[j]
            if checkAnchorOverlap(nr[0], nr[1], nrj[0], nrj[1] ):
                skips.add(j)
                nr = mergeAnchor( nr[0], nr[1], nrj[0], nrj[1] )
        nrs.append(nr)
    return nrs


def getAnchors(loops):
    """
    one loop: ('chrX-chrX-126', ['chrX', 55458553, 55459622, 'chrX', 55514054, 55515304])
    """
    anchors = []
    for loopid,r in loops.items():
        anchors.append( [r[1],r[2]] )
        anchors.append( [r[4],r[5]] )
    while True:
        nrs = mergeAllAnchors(anchors)
        if len(nrs) == len(anchors):
            break
        else:
            anchors = nrs
    return anchors


def getAnchorPETs(jdf,loops,pre,cut=0):
    """
    @param jdf,str, file of .jd
    @param loops: dict, 'chr8-chr8-605': ['chr8', 61242502, 61242734, 'chr8', 61244107, 61244150]
    """
    anchors = getAnchors( loops )
    key, mat = parseJd(jdf, cut)
    report ="%s:%s & %s loops,merged %s anchors"%(key,jdf,len(loops),len(anchors))
    logger.info( report )
    xs_keys, xs = getCorLink(mat[:, 1])
    ys_keys, ys = getCorLink(mat[:, 2])
    ps = set()
    for r in anchors:
        #left end
        l_idx = np.searchsorted(xs_keys, r[0],side="left")
        r_idx = np.searchsorted(xs_keys, r[1],side="right")
        for i in range(l_idx, r_idx):
            ps.update( xs[xs_keys[i]] )
        #right end
        l_idx = np.searchsorted(ys_keys, r[0],side="left")
        r_idx = np.searchsorted(ys_keys, r[1],side="right")
        for i in range(l_idx, r_idx):
            ps.update( ys[ys_keys[i]] )
    nmat = mat[list(ps),]
    joblib.dump(nmat, os.path.join(pre,"-".join(key)+".jd"))
    report ="%s:%s raw PETs %s PETs in anchors"%(key,mat.shape[0],nmat.shape[0])
    logger.info( report )
    return len(loops),len(anchors),mat.shape[0],nmat.shape[0]
 

def jd2cleanWashuPETs(f,dir,sig,pre,chroms=[],cpu=1):
    records = preDs(f,dir,sig,chroms=chroms)
    ds = Parallel(n_jobs=cpu)(delayed(getAnchorPETs)(records[key]["f"] , records[key]["rs"], pre) for key in records.keys())
    #stat PETs in anchors
    l,a,n,m, = 0,0,0,0
    for d in ds:
        l += d[0]
        a += d[1]
        n += d[2]
        m += d[3]
    r = m / 1.0 / n
    report = "%s\t%s,loops:%s, anchors:%s,raw PETs: %s, PETs in anchors:%s, ratio:%s"%( f,dir,l,a,n,m,r )
    logger.info( report )

if __name__ == "__main__":
    start = datetime.now()
    fn = os.path.join(os.getcwd(), "%s.log"%NAME)
    logger = getLogger(fn)
    op = jd2cleanWashuPETsHelp()
    if op.chroms == "":
        chroms = []
    else:
        chroms = set(op.chroms.split(","))
    if not os.path.exists( op.output ):
        os.mkdir( op.output )
    jd2cleanWashuPETs( op.f,op.d,op.significant,op.output,chroms=chroms,cpu=op.cpu) 
    usedtime = datetime.now() - start
    r = "%s finished. Bye! Used time:%s" %(NAME, usedtime)
    logger.info(r)
