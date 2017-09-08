#!/usr/bin/env python2.7
#--coding:utf-8 --

"""
chiaIntersCmpAll.py
2015-07-18: Modified as caculating std.
2015-08-25: Modified as cDBSCAN.
2017-05-22: Modified for cDBSCAN vs. kDBSCAN.
"""


__author__="CAO Yaqiang"
__date__="2017-05-22"
__modified__=""
__email__="caoyaqiang0410@gmail.com"




#general library
import glob,os,time,sys

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
import pylab
import numpy as np
import pandas as pd
from joblib import Parallel,delayed
from sklearn.cluster import DBSCAN
from sklearn import metrics
import brewer2mpl
colors =  brewer2mpl.get_map( 'Set3', 'qualitative',8 ).mpl_colors


#my own library
from Biolibs.rel.General.logger import getlogger
from cDBSCAN import cDBSCAN
from kDBSCAN import kDBSCAN


#global settings
#logger
date=time.strftime(' %Y-%m-%d',time.localtime( time.time())) 
logger = getlogger( fn=os.getcwd()+"/"+date.strip()+"_"+os.path.basename(__file__)+".log" )


def flush( i ):
    if i%1000 == 0:
        report = "\r%dk pairs parsed"%( i/1000 )
        sys.stdout.write( report )
        sys.stdout.flush(  )


def readMat(dir="GM12878_CTCF/"):
    data = {}
    fs = glob.glob(dir+"*.bedpe")
    for f in fs:
        for i, line in enumerate(open(f)):
            flush(i)
            line = line.split("\n")[0].split("\t")
            t = [line[0], line[3]]
            t.sort()
            t = tuple(t)
            if t not in data:
                data[t] = []
            s1 = int( line[ 1 ] )
            e1 = int( line[ 2 ] )
            s2 = int( line[ 4 ] )
            e2 = int( line[ 5 ] )
            s = ( s1+e1 ) /2
            e = ( s2+e2 )/2
            rid = int( line[ 6 ] )
            d = [ rid,s,e ]
            data[t].append(d)
        print
    return data



def getADI( a,b,minPts ):
    for label in set( a.values ):
        c = a[ a==label ]
        if c.shape[ 0 ] < minPts:
            a[ a==label ] = -1
    for label in set( b.values ):
        c = b[ b==label ]
        if c.shape[ 0 ] < minPts:
            b[ b==label ] = -1
    a = a[ a!=-1 ]
    b = b[ b!=-1 ]
    ind = list( set( a.index ).union( b.index ) )
    na = pd.Series( np.zeros( len( ind ) )-1,index=ind )
    na[ a.index ] = a.values
    nb = pd.Series( np.zeros( len( ind ) )-1,index=ind )
    nb[ b.index ] = b.values
    m = metrics.adjusted_rand_score(nb,na)
    c1 = len( set( a.values ) )
    c2 = len( set( b.values ) )
    return m,c1,c2



def runCompare( key,mat,eps,minPts,repeat=5 ):
    report = "Clustering %s and %s using eps as %s minPts as %s" % (key[0], key[1], eps,minPts)
    logger.info(report)
    mat = np.array( mat )
    nmat = pd.DataFrame(mat[ :,1: ], index=mat[ :,0 ], columns=["X", "Y"])
    t1,t2,trs,ms,c1s,c2s = [  ],[  ],[  ],[  ],[  ],[  ]
    for i in xrange( 0,repeat ):
        s = time.clock(  )
        db = cDBSCAN(mat, eps, minPts)
        e = time.clock()
        t = e - s
        t1 += [ t ]
        s = time.clock(  )
        db2 = kDBSCAN(mat, eps, minPts)
        e = time.clock()
        t = e - s
        t2 += [ t ]
        trs += [ t1[ -1 ] / t2[ -1 ] ]
        #ADI
        a = pd.Series(db.labels)
        b = pd.Series(db2.labels)
        m,c1,c2 = getADI( a,b,minPts )
        ms += [ m ]
        c1s += [ c1 ]
        c2s += [ c2 ] 
    t1 = ",".join( map( str,t1 ) )
    t2 = ",".join( map( str,t2 ) )
    trs = ",".join( map( str,trs ) )
    ms = ",".join( map( str,ms ) )
    c1s = ",".join( map( str,c1s ) )
    c2s = ",".join( map( str,c2s ) )
    logger.info("Clustering finished for %s and %s! Time ratio as %s, ADI as %s" % (key[0], key[1],trs,ms) )
    return t1, t2, trs, ms, c1s, c2s


def bench( cs=None,fn = "bench_iDBSCAN.txt",eps=750,minPts=5,repeat=5):
    if cs == None: return
    data = readMat("GM12878_CTCF/")
    keys = [  ]
    for c in cs:
        key = ( c,c )
        if key in data.keys(  ):
            keys.append( ( c,c ) )
    ds = Parallel(n_jobs=-1)(delayed(runCompare)(key,data[ key ],eps,minPts,repeat)
                             for key in keys)
    data = {  }
    for i,d in enumerate(ds):
        key = "%s-%s"%( keys[ i ][ 0 ],keys[ i ][ 1 ] )
        #adi is short for adjusted rand score
        data[key] = {
            "cDBSCAN_Time":d[ 0 ],
            "kDBSCAN_Time":d[ 1 ],
            "TimeRatio": d[ 2 ],
            "ads": d[3],
            "cDBSCAN_cluster":d[ 4 ],
            "kDBSCAN_cluster":d[ 5 ],
        }
    mat = pd.DataFrame(data).T
    mat.to_csv(fn,
            sep="\t",
            index_label="chr")


def plotBench(cs, fn="bench_cDBSCAN.txt"  ):
    data = pd.read_table( fn,index_col=0 )
    keys = [  ]
    for c in cs:
        key = "%s-%s"%( c,c )
        if key in data.index:
            keys.append( key )
    data = data.loc[ keys,: ]
    data.index = [ i.split( "-" )[ 0 ] for i in data.index ]
    fig, ax = pylab.subplots( )
    width = 0.6
    samples = map( str,data.index )
    ind = np.arange( len( data.index ) )
    rs = data.loc[ :,"TimeRatio" ].values
    nrs = [  ]
    for r in rs:
        nrs.append( map( float,r.split( "," ) ) )
    rs = np.array( nrs )
    rs = 1.0/rs
    rs = np.log2(rs)
    #rs = 0.0-np.log10( rs )
    rsm = rs.mean( axis=1 )
    print rsm
    rst = rs.std( axis=1 )
    ax.bar( ind,rsm,width=width,yerr=rst,color=colors[ 0 ],error_kw={ "ecolor":"k" } )
    ax.set_ylabel( "Time of kDBSCAN/cDBSCAN,log2" )
    ax.set_xticks( ind + width/2 )
    ax.set_xticklabels( samples,rotation=45,ha="right" )
    ax.set_xlabel( "Chromosomes sorted by descending PETs counts" )
    #ax.set_ylim( [ 0.0,0.02 ] )
    ax2 = ax.twinx(  )
    umr = data.loc[ :,"ads" ].values
    nrs = [  ]
    for r in umr:
        nrs.append( map( float,r.split( "," ) ) )
    nrs = np.array( nrs )
    umrm = nrs.mean( axis=1 )
    umrt = nrs.std( axis=1 )
    ax2.errorbar( ind+width/2,umrm,yerr=umrt,ecolor="k",color="b",marker="s",label="ARS between cDBSCAN and kDBSCAN" )
    for t in ax2.get_yticklabels(  ):
        t.set_color( "r" )
    ax2.set_ylim( [ 0.9,1.1 ] )
    ax2.set_ylabel( "Adjusted Rand Score (ARS)" )
    leg = ax2.legend( loc="upper right",fancybox=True )
    #leg.get_frame(  ).set_alpha( 0.5 )
    pylab.savefig( "%s.pdf"%fn.split( "." )[ 0 ],dpi=1000,bbox_inches="tight" )




def main(  ):
    #sorted by chromsize,from big to small
    cs = "chrM,chr21,chr22,chr20,chr19,chr16,chr18,chr15,chr17,chr14,chr13,chr9,chrX,chr10,chr8,chr11,chr12,chr7,chr6,chr5,chr4,chr3,chr2,chr1"
    #cs = "chr2,chr1"
    #cs = "chrM,chr21,chr22,chr20,chr19,chr16,chr18,chr15,chr17,chr14,chr13,chr9,chrX,chr10,chr8,chr11,chr12,chr7,chr6,chr5,chr4,chr3"
    cs = cs.split(",")
    cs.reverse()
    eps = 750
    minPts = 5
    #bench(cs=cs,fn="bench_cDBSCAN_chr1chr2.txt",eps=eps,minPts=minPts,repeat=2  )
    #bench(cs=cs,fn="bench_cDBSCAN.txt",eps=eps,minPts=minPts  )
    plotBench(cs, fn="bench_cDBSCAN.txt" )
    


if __name__ == "__main__":
    main(  )
