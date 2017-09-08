#!/usr/bin/env python2.7
#--coding:utf-8 --
"""
"""
from collections import Counter
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree as KDTree

#from scipy.spatial import KDTree as KDTree


class kDBSCAN:
    """
    KDTree based DBSCAN.
    """

    def __init__(self, mat, eps, minPts):
        """
        @param mat: the raw or normalized [pointId,X,Y] data matrix
        @type mat : np.array

        @param eps: The clustering distance threshold, key parameter in DBSCAN.
        @type eps: float

        @param minPts: The min point in neighbor to define a core point, key 
                parameter in DBSCAN.
        @type minPts: int

        """
        #: build the data in the class for global use
        self.eps = eps
        self.minPts = minPts
        self.mat, self.tree, self.cs = self.buildTree(mat)
        self.callDBScan()

    def buildTree(self, mat):
        """
        Build the query tree using KDTree.
        """
        nmat = pd.DataFrame({"x": mat[:, 1], "y": mat[:, 2]}, index=mat[:, 0])
        tree = KDTree(nmat.values)
        cs = pd.Series(-np.ones(mat.shape[0]), index=nmat.index)
        return nmat, tree, cs

    def callDBScan(self):
        """
        Do DBSCAN clustering by go through all points in the sets.
        """
        #: clustering id, noise is -2 and unclassified point is -1.
        clusterId = 0
        for key in self.cs.index:
            if self.cs[key] == -1:
                if self.expandCluster(key, clusterId):
                    clusterId += 1
        #remove the noise and unclassified points
        cs = self.cs
        cs = cs[cs != -1]
        for key in set(cs.values):
            t = cs[cs == key]
            if len(t.index) < self.minPts:
                cs = cs[cs != key]
        self.labels = cs

    def expandCluster(self, pointKey, clusterId):
        """
        Search connection for given point to others.
        
        @param pointKey: the key in self.dataPoints
        @type pointKey: 
       
        @param clusterId: the cluster id for the current
        @type clusterId: int

        @return: bool
        """
        seeds = self.regionQuery(pointKey)
        if len(seeds) < self.minPts:
            self.cs[pointKey] = -2
            return False
        else:
            for key in seeds:
                self.cs[key] = clusterId
            while len(seeds) > 0:
                currentP = seeds[0]
                result = self.regionQuery(currentP)
                if len(result) >= self.minPts:
                    for key in result:
                        if self.cs[key] in [-1, -2]:
                            if self.cs[key] == -1:
                                seeds.append(key)
                            self.cs[key] = clusterId
                del (seeds[0])
            return True

    def regionQuery(self, pointKey):
        """
        Find the related points to the queried point, city block distance is used,based on kdtree
        
        @param pointKey: the key in self.dataPoints
        @type pointKey:
        
        @return: list
        """
        p = self.mat.loc[pointKey, ].values
        rs = self.tree.query_ball_point(p, self.eps, p=1)
        result = list(self.mat.index[rs])
        return result
