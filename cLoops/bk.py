def checkOneEndOverlap(xa, xb, ya, yb):
    """
    check the overlap of a region for the same chromosome
    """
    if ya <= xa <= yb or ya <= xb <= yb:
        return True
    if xa <= ya <= xb or xa <= yb <= xb:
        return True
    return False


def checkOverlap(ra, rb):
    """
    check the overlap of two anchors,ra=[chr,left_start,left_end,chr,right_start,right_end]
    """
    if checkOneEndOverlap(ra[1], ra[2], rb[1], rb[2]) and checkOneEndOverlap(
            ra[4], ra[5], rb[4], rb[5]):
        return True
    return False


def mergeOverlap(ra, rb):
    """
    merge the overlaped two anchors,ra=[chr,left_start,left_end,right_start,right_end]
    """
    nr = [
        ra[0],
        min([ra[1], rb[1]]),
        max([ra[2], rb[2]]), ra[3],
        min([ra[4], rb[4]]),
        max([ra[5], rb[5]])
    ]
    if nr[2] <= nr[1] or nr[5] <= nr[4]:
        return None
    return nr


def filterClusterByDis(data, cut):
    """
    Filter inter-ligation clusters by distances
    """
    for key in data:
        nr = []
        for r in data[key]["records"]:
            d = (r[4] + r[5]) / 2 - (r[1] + r[2]) / 2
            if d >= cut:
                nr.append(r)
        data[key]["records"] = nr
    return data

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

