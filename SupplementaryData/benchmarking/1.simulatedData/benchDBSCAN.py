#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
benchDBSCAN.py
2015-07-18: Modified as add replicates to get std.
"""

__author__ = "CAO Yaqiang"
__date__ = "2015-04-15"
__modified__ = ""
__email__ = "caoyaqiang0410@gmail.com"

#general library
import time

#3rd library
import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4, 2.75)
mpl.rcParams["figure.dpi"] = 100
mpl.rcParams["savefig.transparent"] = True
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["font.size"] = 10.0
mpl.rcParams["font.sans-serif"] = "Arial"
mpl.rcParams["savefig.format"] = "pdf"
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
import pylab
import brewer2mpl
#import seaborn as sns
colors =  brewer2mpl.get_map( 'Set3', 'qualitative',8 ).mpl_colors

#own library
from cDBSCAN import cDBSCAN
from kDBSCAN import kDBSCAN


def getDataD(r=1, samples=1000, std1=0.2, centercounts=10,rmax=15):
    centers = []
    for i in xrange(centercounts):
        x,y = np.random.randint(0-rmax,rmax),np.random.randint(0-rmax,rmax)
        centers.append([x,y])
    #data as true clusters
    X, xlabels = make_blobs(n_samples=samples,
                            centers=centers,
                            cluster_std=std1,
                            random_state=0)

    fig,ax = pylab.subplots()
    ax.scatter(X[:,0],X[:,1],color="red")
    Y = []
    if r > 0:
        for i in xrange(int(samples*r)):
            x,y = np.random.random()*rmax,np.random.random()*rmax
            rx = np.random.random()
            ry = np.random.random()
            if rx < 0.5:
                x = 0-x
            if ry < 0.5:
                y = 0-y
            Y.append([x,y])
        ylabels = np.zeros(len(Y)) -1.0
    Y = np.array(Y)
    ax.scatter(Y[:,0],Y[:,1],color="gray")
    pylab.savefig("data.pdf")



def getData(r=1, samples=1000, std1=0.2, centercounts=100,rmax=5000):
    centers = []
    for i in xrange(centercounts):
        x,y = np.random.randint(0-rmax,rmax),np.random.randint(0-rmax,rmax)
        centers.append([x,y])
    #data as true clusters
    X, xlabels = make_blobs(n_samples=samples,
                            centers=centers,
                            cluster_std=std1,
                            random_state=0)

    Y = []
    if r > 0:
        for i in xrange(int(samples*r)):
            x,y = np.random.random()*rmax,np.random.random()*rmax
            rx = np.random.random()
            ry = np.random.random()
            if rx < 0.5:
                x = 0-x
            if ry < 0.5:
                y = 0-y
            Y.append([x,y])
        ylabels = np.zeros(len(Y)) -1.0
    Y = np.array(Y)
    if r > 0:
        mat = np.concatenate((X, Y), axis=0)
    else:
        mat = X
    x = [  ]
    for i in xrange(len( mat )):
        x.append( [ i,mat[ i,0 ],mat[ i,1 ] ] )
    x = np.array( x )
    mat = pd.DataFrame({"X": mat[:, 0], "Y": mat[:, 1]}, )
    xlabels = list(xlabels)
    #if r > 0:
    #    xlabels.extend(list(ylabels))
    return mat, xlabels,x



def runCompare(samples, r, eps=0.2, minSamples=5, repeat=3):
    mat, labels,x = getData(r=r, samples=samples, std1=eps)
    rs = [  ]
    for i in xrange( 0,repeat ):
        s = time.clock(  )
        db = cDBSCAN(x, eps, minSamples)
        e = time.clock()
        t1 = e - s
        s = time.clock(  )
        db2 = kDBSCAN(x, eps, minSamples)
        e = time.clock()
        t2 = e - s
        tr = t1 / t2
        rs += [ tr ]
    rs = ",".join( map( str,rs ) )
    labela = pd.Series(db.labels)
    labelb = pd.Series(db2.labels)
    na = pd.Series(np.zeros(len(labels))-1.0)
    nb = pd.Series(np.zeros(len(labels))-1.0)
    for i in na.index:
        if i in labela.index:
            na[i] = labela[i]
    for i in nb.index:
        if i in labelb.index:
            nb[i] = labelb[i]
    m1 = metrics.adjusted_rand_score(labels, na)
    m2 = metrics.adjusted_rand_score(labels, nb)
    c1 = len(set(labela.values))
    c2 = len(set(labelb.values)) 
    print "Ratio as %s finished! Time ratio as %s,m1 :%s,m2:%s" % (r, rs,m1,m2)
    return rs, m1, m2, c1, c2


def bench( fn = "bench_cDBSCAN.txt" ):
    data = {}
    rs = np.arange(0.0, 1, 0.1)
    rs = list(rs)
    rs.extend(list(np.arange(1, 20, 2)))
    rs.extend(list(np.arange(20, 56, 5)))
    samples = 10000
    ds = Parallel(n_jobs=-1)(delayed(runCompare)(samples, r, 0.2, 5, 10)
                             for r in rs)
    for i, d in enumerate(ds):
        #adi is short for adjusted rand score
        data[rs[i]] = {
            "TimeRatio": d[0],
            "ads_cDBSCAN": d[1],
            "ads_kDBSCAN": d[2],
            "clusters_cDBSCAN": d[3],
            "clusters_kBSCAN": d[4]
        }
    mat = pd.DataFrame(data).T
    mat.to_csv(fn,
               sep="\t",
               index_label="Noise/Signal Ratio")



def plotBench( fn = "bench_cDBSCAN.txt" ):
    data = pd.read_table( fn,index_col=0 )
    #fig, ax = pylab.subplots( figsize=( 10,6 ) )
    fig, ax = pylab.subplots( )
    width = 0.6
    samples = map( str,data.index )
    ind = np.arange( len( data.index ) )
    rs = data.loc[ :,"TimeRatio" ].values
    nrs = [  ]
    for r in rs:
        nrs.append( map( float,r.split( "," ) ) )
    nrs = np.array(nrs)
    rs = 1.0/nrs
    #rs = np.log2( nrs )
    rsm = rs.mean( axis=1 )
    rst = rs.std( axis=1 )
    ax.bar( ind,rsm,yerr=rst,width=width,color=colors[ 0 ],error_kw={ "ecolor":"k" } )
    #ax.set_ylabel( "log2(Time of kDBSCAN / cDBSCAN)" )
    ax.set_ylabel( "Time of kDBSCAN / cDBSCAN" )
    ax.set_xticks( ind + width/2 )
    ax.set_xticklabels( samples,rotation=45,ha="right" )
    ax.set_xlabel( "Noise/Signal Ratio" )
    #ax.hlines( 1,0,len( ind ) )
    ax2 = ax.twinx(  )
    umr = data.loc[ :,"ads_cDBSCAN" ].values
    mmr = data.loc[ :,"ads_kDBSCAN" ].values
    ax2.plot( ind+width/2,umr,color="b",marker="s",label="ARS for cDBSCAN" )
    ax2.plot( ind+width/2,mmr,color="r",marker="h",label="ARS for kDBSCAN" )
    for t in ax2.get_yticklabels(  ):
        t.set_color( "r" )
    #ax2.set_ylim( [ 0.986,0.988 ] )
    ax2.set_ylabel( "Adjusted Rand Score (ARS)" )
    leg = ax2.legend( loc="best",fancybox=True )
    #leg.get_frame(  ).set_alpha( 0.5 )
    pylab.savefig( "%s.pdf"%fn.split( "." )[ 0 ],dpi=1000,bbox_inches="tight" )

#getDataD(r=2,samples=200)
#bench(fn = "bench_cDBSCAN.txt")
plotBench( fn = "bench_cDBSCAN.txt" )
