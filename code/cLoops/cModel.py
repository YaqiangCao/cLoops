#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
Stastical significance is tested for every chromosome using the local permutated background.
2017-03-24: try except added for estimating reads in a region just in case.
2017-06-07: fix small bug as all self-ligation clusters due to distance cutoff ; also change the PETs of a chromosome to the filtered numbers according to distance cutoff.
2017-06-13: modified ES caculation in case the mean ES is zero, especially for larger window size shifting
2017-06-15: modified combined p-values as geometric mean
2017-06-20: modified the hypergeometric test rab to rab-1, all as following to avoid p-values too small; if no noise nearby, the ES is set to rab;combined-p-value is cancled; pending change N to N' share the same distance distribution,this will lead to too much burden on computation.
2017-07-19: modified HiChIP significant cutoff, for all datasets, eps=5000,minPts=20 without twice option is a good parameter.
2017-08-02: modified genomecoverage model, not using HTSeq now to save memories, also, much faster than previouse method
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
from cLoops.io import parseJd
from cLoops.utils import cFlush


def getGenomeCoverage(f, cut=0):
    """
    Build the genomic model for random access. Could use a lot of memory.
    @param f:.jd file 
    @param cut: distance cutoff for self-ligation PETs.
    """
    key, mat = parseJd(f, cut)
    j = mat.shape[0]
    if j == 0:
        return None, 0
    m = max([np.max(mat[:, 1]), np.max(mat[:, 2])])
    model = [False] * (m + 1000000)  #+10 just in case boundary escape
    for t in mat:
        if model[t[1]] == False:
            model[t[1]] = []
        if model[t[2]] == False:
            model[t[2]] = []
        model[t[1]].append(t[0])
        model[t[2]].append(0 - t[0])
    return model, j * 2


def getCounts(iv, model):
    """
    Get reads ids for a region.
    """
    rs = []
    for i in xrange(iv[0], iv[1]):
        if model[i] != False:
            rs.extend(model[i])
    rs = list(np.abs(list(rs)))
    return rs


def getPETsforRegions(iva, ivb, model):
    ra = getCounts(iva, model)
    rb = getCounts(ivb, model)
    rab = set(ra).intersection(set(rb))
    return len(ra), len(rb), len(rab)


def getNearbyPairRegions(iva, ivb, win=6):
    """
    @param iva: [start,end] 
    Get the nearby regions for interacting two locus,win as how many nearby, 6 is enough for interacting more than 100 regions to estimate FDR and others. The mean distance of all the permutated regions is the same to that between iva and ivb.
    """
    ivas, ivbs = [], []
    ca = sum(iva) / 2
    cb = sum(ivb) / 2
    sa = (iva[1] - iva[0]) / 2
    sb = (ivb[1] - ivb[0]) / 2
    step = min([sa, sb])
    #step = max([sa, sb])
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


def getMultiplePsFdr(iva, ivb, model, N, win=6):
    """
    for the interval a and b, searching its nearby windows to estimate FDR and p-values. THe idea that using  matched nearby windows, which could have similar distance with a & b, needs too many windows. 
    return ra, rb, rab, es, fdr, hyp, chyp, pop, nbp
    """
    ra, rb, rab = getPETsforRegions(iva, ivb, model)
    #simple hypergeometric test, the idea using cis_a + cis_b + trans_a+trans_b as M and cis_a+cis_b as N fails with all p-value as 1
    hyp = hypergeom.sf(rab - 1.0, N, ra, rb)
    ivas, ivbs = getNearbyPairRegions(iva, ivb, win=win)
    hyps, rabs, nbps = [], [], []
    for na in ivas:
        try:
            nra = getCounts(na, model)
        except:
            continue
        nralen = float(len(nra))
        if nralen == 0:
            continue
        for nb in ivbs:
            try:
                nrb = getCounts(nb, model)
            except:
                continue
            if len(nrb) == 0:
                continue
            nrab = len(set(nra).intersection(set(nrb)))
            #collect the value for poisson test
            rabs.append(nrab)
            #collect the nearby hypergeometric test result
            nhyp = hypergeom.sf(nrab - 1.0, N, nralen, len(nrb))
            hyps.append(nhyp)
            #collect the possibility for following binomal test
            den = nrab / (nralen * len(nrb))
            nbps.append(den)
    if len(rabs) == 0:
        return ra, rb, rab, np.inf, 0.0, hyp, 0.0, 0.0, 0.0,
    hyps, rabs = np.array(hyps), np.array(rabs)
    #local fdr
    fdr = len(rabs[rabs > rab]) / float(len(rabs))
    mrabs = float(np.mean(rabs))
    #enrichment score
    if mrabs > 0:
        es = rab / mrabs
    else:
        es = np.inf
    #es = rab / max([np.mean(rabs),float(np.percentile(rabs,90))])
    #es = rab / float(np.percentile(rabs,90))
    #corrected hypergeometric fdr
    chyp = len(hyps[hyps < hyp]) / float(len(hyps))
    #simple possion test, the idea benefits from MACS as using dynamic lambda
    lam = mrabs
    pop = poisson.sf(rab - 1.0, lam)
    #simple binomal test
    bp = np.mean(nbps) * ra * rb / N
    #nbp = binom.sf(rab, N, bp)
    nbp = binom.sf(rab - 1.0, N - rab, bp)
    return ra, rb, rab, es, fdr, hyp, chyp, pop, nbp


def getBonPvalues(ps):
    """
    Return the Bonferroni corrected p-values.
    """
    ps = np.array(ps)
    ps = ps * len(ps)
    ps[ps > 1.0] = 1.0
    return ps


def getIntSig(f, records, minPts, discut):
    """
    @param:discut, distance cutoff determined for self-ligation pets.
    """
    print "Starting estimate significance for interactions in %s" % f
    model, N = getGenomeCoverage(f, discut)
    print "Genomic coverage model built from %s" % f
    if N == 0:
        print "No cis-PETs parsed as requiring distance cutoff >%s from %s" % (
            discut, f)
        return None
    ds = {}
    i = 0
    for r in records:
        chrom = r[0]
        key = "%s-%s-%s" % (r[0], r[3], i)
        iva = [r[1] - 1, r[2] + 1]
        ivb = [r[4] - 1, r[5] + 1]
        #filter loops
        distance = abs(sum(ivb) / 2.0 - sum(iva) / 2.0)
        if distance < discut:
            continue
        ra, rb, rab = getPETsforRegions(iva, ivb, model)
        #filter clusters contain many self-ligation PETs within distance cutoff
        if rab < minPts:
            continue
        i += 1
        if i % 100 == 0:
            cFlush("%s interaction p-values estimated for %s" % (i, f))
        ra, rb, rab, es, fdr, hyp, chyp, pop, nbp = getMultiplePsFdr(
            iva, ivb, model, N)
        #this part should be furthur modified, as for most ideable data, there are no noise, so the es should be inf, however, not possible
        if es == "None":
            continue
        ds[key] = {
            "distance": distance,
            "ra": ra,
            "rb": rb,
            "rab": rab,
            "ES": es,
            "FDR": fdr,
            "hypergeometric_p-value": hyp,
            "hypergeometric_local_FDR": chyp,
            "poisson_p-value": pop,
            "binomal_p-value": nbp,
            "iva": "%s:%s-%s" % (chrom, iva[0], iva[1]),
            "ivb": "%s:%s-%s" % (chrom, ivb[0], ivb[1])
        }
    #memory usage
    del model
    gc.collect()
    print
    if len(ds.keys()) == 0:
        return None
    ds = pd.DataFrame(ds).T
    ds["poisson_p-value_corrected"] = getBonPvalues(ds["poisson_p-value"])
    ds["binomal_p-value_corrected"] = getBonPvalues(ds["binomal_p-value"])
    ds["hypergeometric_p-value_corrected"] = getBonPvalues(
        ds["hypergeometric_p-value"])
    return ds


def markIntSig(ds, escut=1.0, fdrcut=0.05, hpfdrcut=0.05, gpcut=1e-5):
    """
    gpcut is general p-value cutoff for binomal test, poisson test and hypergeometric test.
    """
    #filter data according to cutoffs
    #larger enrichment score
    a = ds["ES"]
    a = a[a >= escut]
    #smaller FDR
    b = ds.loc[a.index, "FDR"]
    b = b[b <= fdrcut]
    #smaller local hypergeometric test FDR
    c = ds.loc[b.index, "hypergeometric_local_FDR"]
    c = c[c <= hpfdrcut]
    #smaller hypergeometric result
    d = ds.loc[c.index, "hypergeometric_p-value"]
    d = d[d <= gpcut]
    #smaller  poisson or binomal,poisson maybe better
    e = ds.loc[d.index, "poisson_p-value"]
    f = ds.loc[d.index, "binomal_p-value"]
    e = e[e <= gpcut]
    f = f[f <= gpcut]
    rs = e.index.union(f.index)
    ns = pd.Series(data=np.zeros(ds.shape[0]), index=ds.index)
    ns[rs] = 1.0
    ds["significant"] = ns
    return ds


def markIntSigHic(ds, escut=2.0, fdrcut=0.05, pcut=1e-5):
    """
    For HiChIP/HiC data, hypergeometric test is not working, poisson and binomal works well. For mouse data, pcut=1e-3 maybe better
    """
    #filter data according to cutoffs
    #larger enrichment score
    a = ds["ES"]
    a = a[a >= escut]
    #smaller FDR
    b = ds.loc[a.index, "FDR"]
    b = b[b <= fdrcut]
    #smaller poisson and binomal result
    c = ds.loc[b.index, "poisson_p-value"]
    c = c[c <= pcut]
    d = ds.loc[b.index, "binomal_p-value"]
    d = d[d <= pcut]
    e = c.index.intersection(d.index)
    ns = pd.Series(data=np.zeros(ds.shape[0]), index=ds.index)
    ns[e] = 1.0
    ds["significant"] = ns
    return ds
