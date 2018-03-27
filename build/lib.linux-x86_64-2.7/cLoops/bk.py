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
