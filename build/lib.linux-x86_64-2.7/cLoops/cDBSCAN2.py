#--coding:utf-8 --
"""
"""
import bisect


class cDBSCAN:
    """
    The major class of the cDBSCAN algorithm, belong to CAO Yaqiang, CHEN Xingwei, AI Daosheng and CHEN Zhaoxiong.
    2018-02-02: all CHEN Zhaoxiong's improvement, very little difference from previouse.
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
        self.buildGrid(mat)
        #: get the clusters
        self.queryGrid()
        del self.Grid

    def getNearbyCells(self, index):
        """
        Basic funciton 2, 9 grid as searching neghbors, grid width is eps.
        """
        x, y = index[0], index[1]
        keys = [(x, y - 1), (x, y + 1), (x - 1, y), (x + 1, y), (x - 1, y - 1),
                (x - 1, y + 1), (x + 1, y - 1), (x + 1, y + 1)]
        #keys = [(x, y), (x, y - 1), (x, y + 1), (x - 1, y), (x - 1, y - 1),
        #        (x - 1, y + 1), (x + 1, y), (x + 1, y - 1), (x + 1, y + 1),
        #        (x, y + 2), (x, y - 2), (x + 1, y + 2), (x + 1, y - 2),
        #        (x - 1, y + 2), (x - 1, y - 2), (x + 2, y), (x + 2, y + 1),
        #        (x + 2, y - 1), (x - 2, y), (x - 2, y + 1), (x - 2, y - 1)]
        nindex = []
        for key in keys:
            if key in self.Grid:
                nindex.append(key)
        return nindex

    def buildGrid(self, mat):
        """
        Algorithm 1: Construct the grids.
        @param mat: the raw or normalized [pointId,X,Y] data matrix
        """
        Grid = {}
        Gorder = {'x': {}, 'y': {}}
        Gtype = {}
        self.axisindex = {'x': 0, 'y': 1}
        for d in mat:
            #modification
            #rotate coordinate system by 45 degree
            x = d[1] - d[2]
            y = d[1] + d[2]
            nx = int(x / self.cw) + 1
            ny = int(y / self.cw) + 1
            Grid.setdefault((nx, ny), [])
            #gird types {0: sparse cell, 1: crowded cell,
            #2: core cell (crowded cell assigned to a cluster)
            #-1: noise cell or edge cell}
            Grid[(nx, ny)].append([x, y, d[0], -1])
        self.Grid = Grid
        for index, cell in Grid.iteritems():
            Gorder['x'][index] = Grid[index]
            nearpnum = len(cell)
            if nearpnum >= self.minPts:
                #crowded cell
                Gtype[index] = 1
                continue
            for near_index in self.getNearbyCells(index):
                nearpnum += len(Grid[near_index])
            if nearpnum < self.minPts:
                #noise cell or edge cell
                Gtype[index] = -1
            else:
                #sparse cell
                Gtype[index] = 0
        #get pt orders
        noisecell = []
        for index in Grid:
            noiseflag = all([
                Gtype[near_index] == -1
                for near_index in self.getNearbyCells(index)
            ])
            if Gtype[index] == -1 and noiseflag:
                #check noise cell
                #noise cell should not waste time sorting
                noisecell.append(index)
                continue
            Gorder['x'][index].sort(key=lambda x: x[0])
            Gorder['y'][index] = sorted(Grid[index], key=lambda x: x[1])
        #remove noise cell
        for index in noisecell:
            del Grid[index]
            del Gtype[index]
        self.Grid = Grid
        self.Gtype = Gtype
        self.Gorder = Gorder

    def queryGrid(self):
        clusterId = 0
        clusters = {}
        for index, cell in self.Grid.iteritems():
            #omit edge cells and core cells
            if self.Gtype[index] in [-1, 2]:
                continue
            border_pts = {}
            #start cell check
            clusters[clusterId] = []
            if self.Gtype[index] == 1:
                #crowded cell must be core cell in one cluster
                border_pts[index] = self.Grid[index]
            else:
                #sparse cell first round checking
                #omit pt which already assigned to other cluster
                pts = [x for x in cell if x[-1] == -1]
                adjacent_pts, flag = self.getSparseCellNeighbor(pts, index)
                if flag:
                    #at leaset one pt in cell meet requirement
                    for p in pts:
                        p[-1] = clusterId
                    clusters[clusterId].extend(pts)
                    border_pts = adjacent_pts
                else:
                    #skip this cell if no core pt found
                    continue
            #Breadth First Search
            while len(border_pts) > 0:
                nindex = sorted(border_pts.keys())[0]
                #print nindex
                #print sorted(border_pts.keys())
                ncell = self.Grid[nindex]
                if self.Gtype[nindex] == 1:
                    #crowded cell process
                    self.Gtype[nindex] = 2
                    for p in ncell:
                        p[-1] = clusterId
                    clusters[clusterId].extend(ncell)
                    self.updatePtDict(border_pts,
                                      self.getCrowdedCellNeighbor(nindex))
                elif self.Gtype[nindex] == 0:
                    #sparse cell process
                    adjacent_pts, flag = self.getSparseCellNeighbor(
                        border_pts[nindex], nindex)
                    if flag:
                        #at leaset one pt in cell meet requirement
                        for p in ncell:
                            if p[-1] == -1:
                                p[-1] = clusterId
                                clusters[clusterId].append(p)
                        self.updatePtDict(border_pts, adjacent_pts)
                    else:
                        #all queried pt are border pt in cluster
                        for p in border_pts[nindex]:
                            p[-1] = clusterId
                        clusters[clusterId].extend(border_pts[nindex])
                else:
                    #edge cell process
                    #can't expand to nearby cells
                    for p in border_pts[nindex]:
                        p[-1] = clusterId
                    clusters[clusterId].extend(border_pts[nindex])
                del border_pts[nindex]

            #release pts if cluster size is too small
            if len(clusters[clusterId]) < self.minPts:
                for p in clusters[clusterId]:
                    p[-1] = -1
                del clusters[clusterId]
            else:
                clusterId += 1
        self.labels = {}
        #f = open("cD_v4.txt", 'w')
        for cid, cluster_pts in clusters.iteritems():
            #f.write("\t".join(map(lambda x: str(x[-2]), sorted(clusters[cid], key=lambda x: x[-2]))) + "\n")
            for p in cluster_pts:
                self.labels[p[-2]] = cid
        #f.close()

    def getCrowdedCellNeighbor(self, index):
        adj_pts = {}
        for axis in ['x', 'y']:
            for delta in [-1, 1]:
                if axis == 'x':
                    newindex = (index[0] + delta, index[1])
                else:
                    newindex = (index[0], index[1] + delta)
                if newindex not in self.Grid or self.Gtype[newindex] == 2:
                    #crowded cell can't expand to core cell
                    continue
                if delta == -1:
                    edgept = self.Gorder[axis][index][0]
                else:
                    edgept = self.Gorder[axis][index][-1]
                newresult = [
                    x
                    for x in self.binSearchAdjPt(newindex, edgept, axis, delta)
                    if x[-1] == -1
                ]
                if newresult: adj_pts[newindex] = newresult

        edge_pts = self.findEdgePts(index)
        for delta in [(-1, -1), (-1, 1), (1, -1), (1, 1)]:
            newindex = (index[0] + delta[0], index[1] + delta[1])
            if newindex not in self.Grid or self.Gtype[newindex] == 2:
                #crowded cell can't expand to core cell
                continue
            for p in edge_pts[delta]:
                newresult = self.overlapPtList(
                    self.binSearchAdjPt(newindex, p, 'x', delta[0]),
                    self.binSearchAdjPt(newindex, p, 'y', delta[1]))
                if self.Gtype[newindex] == 1 and len(newresult) > 0:
                    #In crowded cell one hit equals to all hit
                    adj_pts[newindex] = self.Grid[newindex]
                    break
                if newindex in adj_pts:
                    pre_ids = set([x[-2] for x in adj_pts[newindex]])
                    #for x in newresult:
                    #    if x[-2] not in pre_ids and x[-1] == -1:
                    #        adj_pts[newindex].append(x)
                    adj_pts[newindex].extend([
                        x for x in newresult
                        if x[-2] not in pre_ids and x[-1] == -1
                    ])
                else:
                    newresult = [x for x in newresult if x[-1] == -1]
                    if newresult: adj_pts[newindex] = newresult
        return adj_pts

    def findEdgePts(self, index):
        #find edge pts in crowded cell
        order = {'x': self.Gorder['x'][index], 'y': self.Gorder['y'][index]}
        upleft = [order['x'][0]]
        downleft = [order['x'][0]]
        upflag = True
        downflag = True
        for i in order['x'][1:]:
            if upflag:
                j = upleft[-1]
                if i[1] > j[1]:
                    if i[0] == j[0]:
                        upleft[-1] = i
                    else:
                        upleft.append(i)
                if i[1] == order['y'][-1][1]:
                    upflag = False
            if downflag:
                j = downleft[-1]
                if i[1] < j[1]:
                    if i[0] == j[0]:
                        downleft[-1] = i
                    else:
                        downleft.append(i)
                if i[1] == order['y'][0][1]:
                    downflag = False
            if not (upflag or downflag):
                break
        upright = [order['x'][-1]]
        downright = [order['x'][-1]]
        upflag = True
        downflag = True
        for i in order['x'][-1::-1]:
            if upflag:
                j = upright[-1]
                if i[1] > j[1]:
                    if i[0] == j[0]:
                        upright[-1] = i
                    else:
                        upright.append(i)
                if i[1] == order['y'][-1][1]:
                    upflag = False
            if downflag:
                j = downright[-1]
                if i[1] < j[1]:
                    if i[0] == j[0]:
                        downright[-1] = i
                    else:
                        downright.append(i)
                if i[1] == order['y'][0][1]:
                    downflag = False
            if not (upflag or downflag):
                break
        return {
            (-1, -1): downleft,
            (-1, 1): upleft,
            (1, -1): downright,
            (1, 1): upright
        }

    def getSparseCellNeighbor(self, seedpts, index):
        cell_pt_num = len(self.Grid[index])
        adj_pts_order = {}
        totalresult = {}
        ####
        #pts = copy.copy(seedpts)
        pts = seedpts[:]
        flag = False
        while pts:
            p = pts.pop()
            p_adjacent = {}
            for delta in [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1),
                          (1, 0), (1, 1)]:
                newindex = (index[0] + delta[0], index[1] + delta[1])
                if newindex not in self.Grid:
                    continue
                if delta[1] == 0:
                    #only need to compare x axis
                    p_adjacent[newindex] = self.binSearchAdjPt(
                        newindex, p, 'x', delta[0])
                elif delta[0] == 0:
                    #only need to compare y axis
                    p_adjacent[newindex] = self.binSearchAdjPt(
                        newindex, p, 'y', delta[1])
                else:
                    #find overlap
                    p_adjacent[newindex] = self.overlapPtList(
                        self.binSearchAdjPt(newindex, p, 'x', delta[0]),
                        self.binSearchAdjPt(newindex, p, 'y', delta[1]))
            n = sum([len(cellpts) for cellpts in p_adjacent.values()])
            if n + cell_pt_num >= self.minPts:
                #meet minPts requirement
                self.updatePtDict(totalresult, p_adjacent, checkPt=True)
                if not flag:
                    #if any pt meet the requirement, other unassigned pts not
                    #already in seedpts should be added to check list
                    seedPtIds = set([x[-2] for x in seedpts])
                    pts.extend([
                        x for x in self.Grid[index]
                        if x[-1] == -1 and x[-2] not in seedPtIds
                    ])
                    flag = True
        return totalresult, flag

    def updatePtDict(self, dictA, dictB, checkPt=False):
        #update dictA by adding items from dictB to dictA
        for index, pts in dictB.iteritems():
            if checkPt:
                pts = [x for x in pts if x[-1] == -1]
            if pts:
                if index in dictA:
                    pre_ids = [x[-2] for x in dictA[index]]
                    #for p in pts:
                    #    if p[-2] not in pre_ids:
                    #        dictA[index].append(p)
                    dictA[index].extend(
                        [x for x in pts if x[-2] not in pre_ids])
                else:
                    dictA[index] = pts

    def binSearchAdjPt(self, index, q_pt, axis, delta):
        #binary search
        #p
        pts = self.Gorder[axis][index]
        if delta == 0:
            return pts
        axispos = self.axisindex[axis]
        posarray = [p[axispos] for p in pts]
        xpos = q_pt[axispos] + self.eps * delta
        if delta == 1:
            index = bisect.bisect_right(posarray, xpos)
            return pts[0:index]
        elif delta == -1:
            index = bisect.bisect_left(posarray, xpos)
            return pts[index:]

    def overlapPtList(self, ptlistA, ptlistB):
        ptdictA = {x[-2]: x for x in ptlistA}
        newkeys = set(ptdictA.keys()) & set([x[-2] for x in ptlistB])
        return [ptdictA[id] for id in newkeys]
