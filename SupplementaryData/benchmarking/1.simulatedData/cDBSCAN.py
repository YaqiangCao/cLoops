#--coding:utf-8 --
"""
"""


class cDBSCAN:
    """
    The major class of the cDBSCAN algorithm, belong to CAO Yaqiang, CHEN Xingwei & AI Daosheng.
    Algorithm complexity analysis is as following:
    buildGrids: 2N,N for find the minal X and Y, N for coordinates convertion.
    removeNoiseGrids:

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
        #: get the clusters
        self.callClusters()
        del self.Gs, self.Gs2, self.ps

    def getDist(self, p, q):
        """
        Basic function 1, city block distance funciton.
        """
        x = self.ps[p]
        y = self.ps[q]
        d = abs(x[0] - y[0]) + abs(x[1] - y[1])
        #euclidean distance ,just in case.
        #d = np.sqrt((x[0] - y[0])**2 + (x[1] - y[1])**2)
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
            Gs.setdefault((nx, ny), [])
            Gs[(nx, ny)].append(d[0])
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

    def callClusters(self):
        """
        Algorithm 4: Do DBSCAN clustering by go through all points in the sets.
        """
        #: clustering id, noise is -2 and unclassified point is -1.
        clusterId = 0
        for key in self.ps:
            if self.ps[key][-1] == -1:
                if self.expandCluster(key, clusterId):
                    clusterId += 1
        #remove the noise and unclassified points
        labels = {}
        cs = {}
        for p in self.ps.keys():
            c = self.ps[p][-1]
            if c == -2:
                continue
            labels[p] = c
            if c not in cs:
                cs[c] = []
            cs[c].append(p)
        for key in cs.keys():
            if len(cs[key]) < self.minPts:
                for p in cs[key]:
                    del labels[p]
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
        if len(seeds) < self.minPts:
            self.ps[pointKey][-1] = -2
            return False
        else:
            for key in seeds:
                self.ps[key][-1] = clusterId
            while len(seeds) > 0:
                currentP = seeds[0]
                result = self.regionQuery(currentP)
                if len(result) >= self.minPts:
                    for key in result:
                        if self.ps[key][-1] in [-1, -2]:
                            if self.ps[key][-1] == -1:
                                seeds.append(key)
                            self.ps[key][-1] = clusterId
                del (seeds[0])
            return True

    def regionQuery(self, pointKey):
        """
        Find the related points to the queried point, city block distance is used.
        
        @param pointKey: the key in self.dataPoints
        @type pointKey:
        
        @return: list
        """
        p = self.ps[pointKey]
        x = p[2]
        y = p[3]
        #scan square and get nearby points.
        result = [pointKey]
        for q in self.Gs2[(x, y)]:
            if q == pointKey:
                continue
            if self.getDist(pointKey, q) <= self.eps:
                result.append(q)
        return result
