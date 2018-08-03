"""
   Convex Hull Assignment: COSC262 (2017)
   Student Name: Matthew Belworthy Lewthwaite
   Usercode: mbe84
"""
from time import time,clock

def readDataPts(filename, N):
    """Reads the first N lines of data from the input file
          and returns a list of N tuples
          [(x0,y0), (x1, y1), ...]
    """
    file = open(filename, 'r')
    lines = file.readlines()
    listPts = []
    for i in range(N):
        current_line = lines[i]
        point = current_line.split()
        listPts.append((float(point[0]),float(point[1])))
    return listPts

def giftwrap(listPts):
    """Returns the convex hull vertices computed using the
          giftwrap algorithm as a list of 'h' tuples
          [(u0,v0), (u1,v1), ...]    
    """
    PtList = listPts[:]
    k = listPts.index(min(listPts, key=lambda t: (t[1], -t[0]))) #this is start point
    i, v, = 0, 0
    PtList.append(PtList[k]) 
    n = len(PtList) 
    chull = [PtList[k]]
    while k != len(listPts):
        PtList[i], PtList[k] = PtList[k], PtList[i]
        minAngle = 361
        for j in range(i + 1, len(PtList)):
            angle = theta(PtList[i], PtList[j])
            if angle < minAngle and angle > v and i != j:
                minAngle = angle
                k = j
        
        chull.append(PtList[k])
        i += 1
        v = minAngle   
    return chull[:-1]



def theta(pointA, pointB):
    """Returns the angle using the theta approximation. If 0 returns 360.
       For purposes of the giftwrap algorithm
       """
    dx = pointB[0] - pointA[0]
    dy = pointB[1] - pointA[1]
    if abs(dx) < 1.e-6 and abs(dy) < 1.6e-6:
        t = 0
    else:
        t = dy/(abs(dx) + abs(dy))
    
    if dx < 0:
        t = 2 - t
    elif dy < 0:
        t = 4 + t
    if t == 0:
        t = 4
    return t * 90

def grahamscan(listPts):
    """Returns the convex hull vertices computed using the
         Graham-scan algorithm as a list of 'h' tuples
         [(u0,v0), (u1,v1), ...]  
    """
    k = listPts.index(min(listPts, key=lambda t: (t[1], -t[0])))
    ptSort = sorted(listPts, key=lambda x: theta(listPts[k], x))
    ptSort.remove(listPts[k])
    ptSort = [listPts[k]] + ptSort
    
    chull = []
    chull.append(listPts[k])
    chull.append(ptSort[1])
    chull.append(ptSort[2])
    for i in range(3, len(ptSort)):
        while not isCCW(chull[-2], chull[-1], ptSort[i]):
            chull.pop()
        if theta(listPts[k], ptSort[i]) == theta(listPts[k], chull[-1]):
            var1 = ptSort[i][0] + ptSort[i][1] - listPts[k][0] - listPts[k][0]
            var2 = chull[-1][0] + chull[-1][1] - listPts[k][0] - listPts[k][0]
            if var1 < var2:
                chull.pop()
                chull.append(ptSort[i])
        else:
            chull.append(ptSort[i])
    return chull

def amethod(listPts):
    """Returns the convex hull vertices computed using the
         Graham-scan algorithm as a list of 'h' tuples
         [(u0,v0), (u1,v1), ...] this 
    """
    points = sorted(listPts, key=lambda x: x[0])
    if len(points) <= 1:
        return points

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    return lower[1:-1] + upper[:-1] + lower[:1]



def lineFn(ptA, ptB, ptC):
    return (ptB[0]-ptA[0]) * (ptC[1] - ptA[1]) - (ptB[1] - ptA[1]) * (ptC[0] - ptA[0])

def isCCW(ptA, ptB, ptC):
    return lineFn(ptA, ptB, ptC) > 0

def main():
    listPts = readDataPts('Set_B.dat', 30000) 
    timerGiftWrap = clock()
    print(giftwrap(listPts))      #You may replace these three print statements
    print("Giftwrap time: ", clock()-timerGiftWrap)
    timerGrahamScan = clock()
    print (grahamscan(listPts))   #with any code for validating your outputs
    print("Graham Scan time: ", clock()-timerGrahamScan)
    timerAMethod = clock()
    print (amethod(listPts)) 
    print("A method time: ", clock()-timerAMethod)

if __name__ == "__main__":
    main()