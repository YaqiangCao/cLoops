#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
Stastical significance is tested for every chromosome using the local permutated background.
2018-02-01: improved data structure for genomecoverage,much faster and less memory than previouse version for significance calling,slightly changed the loops boundary.
2018-03-08: modified ChIA-PET significant loops cutoff
2018-03-16: key change, sliding step changed to the half of the mean anchor size always get more significant loops
2018-03-23: modified merging overlapped loops, significant loops with smaller anchors are first selected
2018-03-26: modified to speed up
2018-03-28: modified the mergeing method, small bugs fixed
"""
__date__ = "2017-03-15"
__modified__ = ""
__email__ = "caoyaqiang@picb.ac.cn chenxingwei@picb.ac.cn"

#general library
import gc

#3rd library
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, binom, poisson, combine_pvalues

#cLoops
from cLoops.io import parseJd, parseIv
from cLoops.utils import cFlush


def getCorLink(cs):
    """
    @param cs: [1,2,3,4], a list for the coordinates x or y
    @rtype: dic, keys is the coordinate, value is the closest next coordinate and points index in this coordinate
    """
    ts = {}
    for i, c in enumerate(cs):
        if c not in ts:
            ts[c] = []
        ts[c].append(i)
    keys = sorted(ts.keys())
    for i in xrange(len(keys) - 1):
        ni = keys[i]
        nj = keys[i + 1]
        ts[ni] = {"next": nj, "points": ts[ni]}
    ts[nj] = {"next": None, "points": ts[nj]}
    return ts


def getGenomeCoverage(f, cut=0):
    """
    Build the genomic model for random access. Could use a lot of memory.
    @param f:.jd file 
    @param cut: distance cutoff for self-ligation PETs.
    """
    key, mat = parseJd(f, cut)
    j = mat.shape[0]
    if j < 2:
        return None, 0
    xs = getCorLink(mat[:, 1])
    ys = getCorLink(mat[:, 2])
    return [xs, ys], j


def getCounts(iv, ts):
    ps = []
    pos = None
    for i in xrange(iv[0], iv[1]):
        if i in ts:
            pos = i
            break
    while pos <= iv[1] and pos != None:
        ps.extend(ts[pos]["points"])
        pos = ts[pos]["next"]
    return set(ps)


def getPETsforRegions(iva, ivb, model):
    raSource = getCounts(iva, model[0])
    raTarget = getCounts(iva, model[1])
    rbSource = getCounts(ivb, model[0])
    rbTarget = getCounts(ivb, model[1])
    ra = len(raSource.union(raTarget))
    rb = len(rbSource.union(rbTarget))
    rab = len(raSource.intersection(rbTarget))
    return ra, rb, rab


def getNearbyPairRegions(iva, ivb, win=5):
    """
    @param iva: [start,end] 
    Get the nearby regions for interacting two locus,win as how many nearby, 6 is enough for interacting more than 100 regions to estimate FDR and others. The mean distance of all the permutated regions is the same to that between iva and ivb.
    """
    ivas, ivbs = [], []
    ca = sum(iva) / 2
    cb = sum(ivb) / 2
    sa = (iva[1] - iva[0]) / 2
    sb = (ivb[1] - ivb[0]) / 2
    step = (sa + sb) / 2
    for i in xrange(0 - win, win + 1):
        if i == 0:
            continue
        niva = [iva[0], iva[1]]
        niva[0] = max([0, ca + i * step - sa])
        niva[1] = max([0, ca + i * step + sa])
        nivb = [ivb[0], ivb[1]]
        nivb[0] = max([0, cb + i * step - sb])
        nivb[1] = max([0, cb + i * step + sb])
        ivas.append(niva)
        ivbs.append(nivb)
    return ivas, ivbs


def getMultiplePsFdr(iva, ivb, model, N, win=5):
    """
    for the interval a and b, searching its nearby windows to estimate FDR and p-values.  
    return ra, rb, rab, es,es_ra,es_rb, fdr, hyp, pop, nbp
    """
    ra, rb, rab = getPETsforRegions(iva, ivb, model)
    hyp = max([1e-300, hypergeom.sf(rab - 1.0, N, ra, rb)])
    ivas, ivbs = getNearbyPairRegions(iva, ivb, win=win)
    #nras is a list for storing points ids for permutated regions
    nras, nrbs = [], []
    for na in ivas:
        nraSource = getCounts(na, model[0])
        nraTarget = getCounts(na, model[1])
        nra = nraSource.union(nraTarget)
        nras.append(nra)
    for nb in ivbs:
        nrbSource = getCounts(nb, model[0])
        nrbTarget = getCounts(nb, model[1])
        nrb = nrbSource.union(nrbTarget)
        nrbs.append(nrb)
    #caculating the permutated background
    rabs, nbps = [], []
    for nra in nras:
        nralen = float(len(nra))
        for nrb in nrbs:
            nrblen = len(nrb)
            nrab = float(len(nra.intersection(nrb)))
            if nrab > 0:
                #collect the value for poisson test
                rabs.append(nrab)
                #collect the possibility for following binomial test
                den = nrab / (nralen * nrblen)
                nbps.append(den)
            else:
                nbps.append(0.0)
                rabs.append(0.0)
    if len(rabs) == 0:
        return ra, rb, rab, np.inf, 0.0, hyp, 0.0, 1e-300, 1e-300,
    rabs = np.array(rabs)
    #local fdr
    fdr = len(rabs[rabs > rab]) / float(len(rabs))
    mrabs = float(np.mean(rabs))
    #enrichment score
    if mrabs > 0:
        es = rab / np.mean(rabs[rabs > 0])
    else:
        es = np.inf
    #simple possion test
    lam = mrabs
    pop = max([1e-300, poisson.sf(rab - 1.0, lam)])
    #simple binomial test
    bp = np.mean(nbps) * ra * rb / N
    nbp = max([1e-300, binom.sf(rab - 1.0, N - rab, bp)])
    return ra, rb, rab, es, fdr, hyp, pop, nbp


def getBonPvalues(ps):
    """
    Return the Bonferroni corrected p-values.
    """
    ps = np.array(ps)
    ps = ps * len(ps)
    ps[ps > 1.0] = 1.0
    return ps


def checkOneEndOverlap(xa, xb, ya, yb):
    """
    check the overlap of a region for the same chromosome
    """
    if (ya <= xa <= yb) or (ya <= xb <= yb) or (ya <= xa <= xb <= yb):
        return True
    if (xa <= ya <= xb) or (xa <= yb <= xb) or (xa <= ya <= yb <= xb):
        return True
    return False


def checkOverlap(ivai,ivbi,ivaj,ivbj):
    """
    check the overlap of two anchors,ra=[chr,left_start,left_end,right_start,right_end]
    """
    if ivai[0] != ivaj[0] or ivbi[0] != ivbj[0]:
        return 
    if checkOneEndOverlap(ivai[1], ivai[2], ivaj[1], ivaj[2]) and checkOneEndOverlap(
            ivbi[1], ivbi[2], ivbj[1], ivbj[2]):
        return True
    return False


def removeDup(ds, bpcut=1e-5):
    """
    Remove overlapped called loops, keep the more significant one for multiple eps result. 
    @param:ds, from getIntSig
    """
    uniqueds = {}
    reds = {}
    rekeys = set()
    keys = ds.keys()
    for i in xrange(len(keys) - 1):
        keyi = keys[i]
        if keyi in rekeys:
            continue
        ivai = parseIv(ds[keyi]["iva"])
        ivbi = parseIv(ds[keyi]["ivb"])
        #1 means unique loops
        flag = 1
        #collect overlapped loops
        for j in xrange(i + 1, len(keys)):
            keyj = keys[j]
            if keyj in rekeys:
                continue
            ivaj = parseIv(ds[keyj]["iva"])
            ivbj = parseIv(ds[keyj]["ivb"])
            flagj = checkOverlap(ivai,ivbi,ivaj,ivbj)
            #there is overlapped loops,collect them
            if flagj:
                if keyi not in reds:
                    reds[keyi] = [keyi]
                    rekeys.add(keyi)
                reds[keyi].append(keyj)
                rekeys.add(keyj)
                flag = 0
        #collect unique loops
        if flag:
            uniqueds[keyi] = ds[keyi]
    #for overlapped loops, choose the more significant ones
    for key in reds.keys():
        ts = {}
        for t in reds[key]:
            if ds[t]["binomial_p-value"] > bpcut:
                continue
            #ts[t] = ds[t]["binomial_p-value"]
            #first select the significant loops, then select the loops with smaller anchors and higher density
            ts[t] = float(ds[t]["rab"]) / ds[t]["ra"] / ds[t]["rb"]
        """
        Used for debugging
            iva = parseIv(ds[t]["iva"])
            ivb = parseIv(ds[t]["ivb"])
            a = iva[2]-iva[1]+ivb[2]-ivb[1]
            b = float(ds[t]["rab"])/ds[t]["ra"]/ds[t]["rb"]
            c = float(ds[t]["rab"])/a
            print t
            print a,b,c,ds[t]["rab"],ds[t]["ra"],ds[t]["rb"],ds[t]["ES"],ds[t]["binomial_p-value"],ds[t]["poisson_p-value"]
        print 
        """
        if len(ts) == 0:
            continue
        ts = pd.Series(ts)
        ts.sort_values(inplace=True, ascending=False)
        uniqueds[ts.index[0]] = ds[ts.index[0]]
    return uniqueds


def getIntSig(f, records, minPts, discut):
    """
    @param:discut, distance cutoff determined for self-ligation pets.
    """
    print "Starting estimate significance for %s candidate interactions in %s" % (
        len(records), f)
    model, N = getGenomeCoverage(f, discut)
    print "Genomic coverage model built from %s" % f
    if N == 0:
        print "No cis-PETs parsed as requiring distance cutoff >%s from %s" % (
            discut, f)
        return None
    #print "records:",len(records) #used for debuging
    ds = {}
    i = 0
    for r in records:
        chrom = r[0]
        key = "%s-%s-%s" % (r[0], r[3], i)
        iva = [max(0, r[1]), r[2]]
        ivb = [max(0, r[4]), r[5]]
        #filter loops
        distance = abs(sum(ivb) / 2.0 - sum(iva) / 2.0)
        if distance < discut:
            continue
        ra, rb, rab = getPETsforRegions(iva, ivb, model)
        #filter clusters contain many self-ligation PETs within distance cutoff
        #if rab < min(minPts):
        if rab < max(minPts):
            continue
        i += 1
        if i % 100 == 0:
            cFlush("%s interaction p-values estimated for %s" % (i, f))
        ra, rb, rab, es, fdr, hyp, pop, nbp = getMultiplePsFdr(
            iva, ivb, model, N)
        #this part should be furthur modified, as for most ideable data, there are no noise, so the es should be inf, however, not possible
        ds[key] = {
            "distance": distance,
            "ra": ra,
            "rb": rb,
            "rab": rab,
            "ES": es,
            "FDR": fdr,
            "hypergeometric_p-value": hyp,
            "poisson_p-value": pop,
            "binomial_p-value": nbp,
            "iva": "%s:%s-%s" % (chrom, iva[0], iva[1]),
            "ivb": "%s:%s-%s" % (chrom, ivb[0], ivb[1])
        }
    #memory usage
    del model
    gc.collect()
    print
    #print "records before remove duplicates:",len(ds) #used for debuging
    if len(ds.keys()) == 0:
        return None
    ds = removeDup(ds)
    if len(ds.keys()) == 0:
        return None
    #print "records after remove duplicates:",len(ds) #used for debuging
    ds = removeDup(ds)
    if len(ds.keys()) == 0:
        return None
    #print "records after remove duplicates again:",len(ds) #used for debuging
    ds = pd.DataFrame(ds).T
    ds["poisson_p-value_corrected"] = getBonPvalues(ds["poisson_p-value"])
    ds["binomial_p-value_corrected"] = getBonPvalues(ds["binomial_p-value"])
    ds["hypergeometric_p-value_corrected"] = getBonPvalues(
        ds["hypergeometric_p-value"])
    return ds


def markIntSig(ds,
               escut=2.0,
               fdrcut=1e-2,
               bpcut=1e-5,
               ppcut=1e-5,
               hypcut=1e-10):
    """
    gpcut is general p-value cutoff for binomial test, poisson test and hypergeometric test.
    """
    #filter data according to cutoffs
    #larger enrichment score
    a = ds["ES"]
    a = a[a >= escut]
    #smaller FDR
    b = ds.loc[a.index, "FDR"]
    b = b[b <= fdrcut]
    #smaller hypergeometric result
    c = ds.loc[b.index, "hypergeometric_p-value"]
    c = c[c <= hypcut]
    #smaller poisson and binomial
    d = ds.loc[c.index, "poisson_p-value"]
    d = d[d <= ppcut]
    e = ds.loc[d.index, "binomial_p-value"]
    e = e[e <= bpcut]
    rs = e.index
    ns = pd.Series(data=np.zeros(ds.shape[0]), index=ds.index)
    ns[rs] = 1.0
    ds["significant"] = ns
    return ds


def markIntSigHic(ds, escut=2.0, fdrcut=0.01, bpcut=1e-5, ppcut=1e-5):
    """
    For HiChIP/HiC data, hypergeometric test is not working, poisson and binomial works well. For mouse data, pcut=1e-3 maybe better
    """
    #filter data according to cutoffs
    #larger enrichment score
    a = ds["ES"]
    a = a[a >= escut]
    #smaller FDR
    b = ds.loc[a.index, "FDR"]
    b = b[b < fdrcut]
    #smaller poisson and binomial result
    c = ds.loc[b.index, "poisson_p-value"]
    c = c[c <= ppcut]
    d = ds.loc[b.index, "binomial_p-value"]
    d = d[d <= bpcut]
    e = c.index.intersection(d.index)
    ns = pd.Series(data=np.zeros(ds.shape[0]), index=ds.index)
    ns[e] = 1.0
    ds["significant"] = ns
    return ds
