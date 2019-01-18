#--coding:utf-8 --
"""
"""
from copy import deepcopy


def checkOneEndOverlap(xa, xb, ya, yb):
    """
    check the overlap of a region for the same chromosome
    """
    if (ya <= xa <= yb) or (ya <= xb <= yb):
        return True
    if (xa <= ya <= xb) or (xa <= yb <= xb):
        return True
    return False


def checkOverlap(ivai, ivbi, ivaj, ivbj,margin=None):
    """
    check the overlap of two anchors,ra=[chr,left_start,left_end,right_start,right_end], if there is small margin
    """
    if checkOneEndOverlap(ivai[0], ivai[1], ivaj[0], ivaj[1]) and checkOneEndOverlap( ivbi[0], ivbi[1], ivbj[0], ivbj[1]): 
        return True
    if margin != None:
        if checkOneEndOverlap(ivai[0], ivai[1], ivaj[0]-margin, ivaj[1]-margin) and checkOneEndOverlap( ivbi[0], ivbi[1], ivbj[0], ivbj[1]): return True
        if checkOneEndOverlap(ivai[0], ivai[1], ivaj[0]+margin, ivaj[1]+margin) and checkOneEndOverlap( ivbi[0], ivbi[1], ivbj[0], ivbj[1]): return True
        if checkOneEndOverlap(ivai[0], ivai[1], ivaj[0]+margin, ivaj[1]+margin) and checkOneEndOverlap( ivbi[0], ivbi[1], ivbj[0]-margin, ivbj[1]-margin): return True
        if checkOneEndOverlap(ivai[0], ivai[1], ivaj[0]+margin, ivaj[1]+margin) and checkOneEndOverlap( ivbi[0], ivbi[1], ivbj[0]+margin, ivbj[1]+margin): return True
    return False


def merge(rs,margin=None):
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
            if checkOverlap([nr[0],nr[1]], [nr[2],nr[3]], [nrj[0],nrj[1]], [nrj[2],nrj[3]],margin):
                skips.add(j)
                nr[0] = min([nr[0],nrj[0]])
                nr[1] = max([nr[1],nrj[1]])
                nr[2] = min([nr[2],nrj[2]])
                nr[3] = max([nr[3],nrj[3]])
                nr[4].extend(nrj[4])
        nrs.append(nr)
    return nrs


def mergeClusters( rs,margin=None ):
    """
    Merge overlapped clusters, 
    rs = [  xmin,xmax,ymin,ymax,[pid1,pid2,pid3] ] 
    """
    i = 0
    while True:
        nrs = merge(rs,margin)
        if len(nrs) == len(rs):
            break
        rs = nrs
    return nrs


class blockDBSCAN:
    """
    The major class of the blockDBSCAN algorithm, belong to CAO Yaqiang.

    """

    def __init__(self, mat, eps, minPts,mergeMargin=False):
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
        self.mergeMargin = mergeMargin
        #: cell width, city block distance
        self.cw = self.eps
        #: build the square index for quick neighbor search
        self.buildGrids(mat)
        #: get the points for all neighbors
        self.buildGridNeighbors()
        #: remove noise grids
        self.removeNoiseGrids()
        #: get the points for all neighbors
        self.buildGridNeighbors()
        #: convert grids into points
        self.centerGrids()
        #: get the clusters
        self.callClusters()
        self.getLabels()
        del self.Gs, self.Gs2, self.Gs3, self.ps

    def getDist(self, x, y):
        """
        Basic function 1, city block distance funciton.
        """
        d = abs(x[0] - y[0]) + abs(x[1] - y[1])
        return d

    def getNearbyGrids(self, cell):
        """
        Basic funciton 2, 9 grid as searching neghbors, grid width is eps.
        """
        x, y = cell[0], cell[1]
        #keys = [(x, y),
        keys = [(x, y - 1), (x, y + 1), (x - 1, y), (x + 1, y), (x - 1, y - 1),
                (x - 1, y + 1), (x + 1, y - 1), (x + 1, y + 1)]
        #keys = [(x, y), (x, y - 1), (x, y + 1), (x - 1, y), (x - 1, y - 1),
        #        (x - 1, y + 1), (x + 1, y), (x + 1, y - 1), (x + 1, y + 1),
        #        (x, y + 2), (x, y - 2), (x + 1, y + 2), (x + 1, y - 2),
        #        (x - 1, y + 2), (x - 1, y - 2), (x + 2, y), (x + 2, y + 1),
        #        (x + 2, y - 1), (x - 2, y), (x - 2, y + 1), (x - 2, y - 1)]
        ncells = []
        for key in keys:
            if key in self.Gs:
                ncells.append(key)
        return ncells

    def buildGrids(self, mat):
        """
        Algorithm 1: Construct the grids.
        @param mat: the raw or normalized [pointId,X,Y] data matrix
        """
        minX, minY = mat[0][1], mat[0][2]
        for t in mat:
            minX = min([minX, t[1]])
            minY = min([minY, t[2]])
        Gs = {}
        ps = {}
        for d in mat:
            nx = int((d[1] - minX) / self.cw) + 1
            ny = int((d[2] - minY) / self.cw) + 1
            Gs.setdefault((nx, ny), []).append(d[0])
            #last elements marks the class, initially -1 as noise
            ps[d[0]] = [d[1], d[2], nx, ny, -1]
        self.Gs, self.ps = Gs, ps

    def buildGridNeighbors(self):
        """
        Algorithm 2 : Grid index with all neighbor points.
        """
        Gs2 = {}
        for cell in self.Gs.keys():
            nps = []
            nps.extend(self.Gs[cell])
            for cellj in self.getNearbyGrids(cell):
                nps.extend(self.Gs[cellj])
            Gs2[cell] = nps
        self.Gs2 = Gs2

    def removeNoiseGrids(self):
        """
        Algorithm 3: Remove noise grid according to KNN and get the obvious core points and core grids.
        """
        #: noise cells without neighbors
        tode = set()
        #: noise cells with neighbors
        tode2 = set()
        for cell in self.Gs.keys():
            if len(self.Gs2[cell]) < self.minPts:
                tode2.add(cell)
        #KNN to noise cells with neighbors
        for cell in tode2:
            cells = self.getNearbyGrids(cell)
            ncells = set(cells) & tode2
            #all neighbor cells are noise
            if len(cells) == len(ncells):
                tode.add(cell)
        for cell in tode:
            for p in self.Gs[cell]:
                del self.ps[p]
            del self.Gs[cell]
    
    def centerGrids(self):
        """
        Algorithm 4: convert each grid into a point, and then cluster grids.
        """
        Gs3 = {}
        for cell in self.Gs.keys():
            pids = self.Gs[cell]
            ps = []
            x, y = 0, 0
            for pid in pids:
                x += self.ps[pid][0]
                y += self.ps[pid][1]
            x = x / len(pids)
            y = y / len(pids)
            Gs3[cell] = [x, y, len(pids), -1]
            #print(cell,x,y,len(pids))
        self.Gs3 = Gs3

    def callClusters(self):
        """
        Algorithm 4: Do DBSCAN clustering by go through all points in the sets.
        """
        #: clustering id, noise is -2 and unclassified point is -1.
        clusterId = 0
        for key in self.Gs3:
            #print("begining",key,self.Gs3[key])
            if self.Gs3[key][-1] == -1:
                if self.expandCluster(key, clusterId):
                    clusterId += 1

    def getLabels(self, ):
        #remove the noise and unclassified points
        cs = {}
        #label each point
        for c in self.Gs3.keys():  #visit each cell
            if self.Gs3[c][-1] == -2:
                continue
            cid = self.Gs3[c][-1]
            for p in self.Gs[c]:
                cs.setdefault(cid, []).append(p)
        #collected clusters and get their boundary
        rs = []
        for cid, ps in cs.items():
            if len(ps) < self.minPts:
                continue
            for i, p in enumerate(ps):
                x = self.ps[p][0]
                y = self.ps[p][1]
                if i == 0:
                    xmin = x
                    xmax = x
                    ymin = y
                    ymax = y
                else:
                    if x < xmin:
                        xmin = x
                    if x > xmax:
                        xmax = x
                    if y < ymin:
                        ymin = y
                    if y > ymax:
                        ymax = y
            rs.append([xmin, xmax, ymin, ymax, ps]) #boundary and points id
        #merge clusters
        if self.mergeMargin:
            rs = mergeClusters( rs,self.eps /2 )
        labels = {}
        for i, r in enumerate(rs):
            for p in r[-1]:
                labels[p] = i
        self.labels = labels

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
        #print("init seeds:",seeds)
        if len(seeds) == 1 and self.Gs3[pointKey][2] < self.minPts:
            self.Gs3[pointKey][-1] = -2
            return False
        else:
            for key in seeds:
                self.Gs3[key][-1] = clusterId
                #print("assign keyid",key,self.Gs3[key])
            while len(seeds) > 0:
                currentP = seeds[0]
                result = self.regionQuery(currentP)
                #print("current and find",currentP,result)
                if len(result) >= 2:
                    for key in result:
                        if self.Gs3[key][-1] in [-1, -2]:
                            if self.Gs3[key][2] >= self.minPts and self.Gs3[
                                    key][-1] == -1:
                                seeds.append(key)
                            self.Gs3[key][-1] = clusterId
                            #print("assigned clusterid",key,self.Gs3[key])
                del (seeds[0])
                #print(seeds)
                #print("\n")
            return True

    def getGridDist(self, keya, keyb):
        """
        """
        for p in self.Gs[keya]:
            x = (self.ps[p][0], self.ps[p][1])
            for q in self.Gs[keyb]:
                y = (self.ps[q][0], self.ps[q][1])
                if self.getDist(x, y) <= self.eps:
                    return True
        return False

    def regionQuery(self, pointKey):
        """
        Find the related points to the queried point, city block distance is used.
        @param pointKey: the key in self.dataPoints
        @type pointKey:
        @return: list
        """
        p = self.Gs3[pointKey]
        x = (p[0], p[1])
        #scan square and get nearby points.
        result = [pointKey]
        for q in self.getNearbyGrids(pointKey):
            if q == pointKey:
                continue
            qq = self.Gs3[q]
            y = (qq[0], qq[1])
            if self.getDist(x, y) <= self.eps:
                result.append(q)
            else:
                if self.getGridDist(pointKey, q):
                    result.append(q)
        return result
