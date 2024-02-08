import matplotlib.pyplot as plt
import numpy as np

class Vector:
  def __init__(self, inicial , final):
    self.inicial = inicial
    self.final = final

class Punto:
    def __init__(self, x , y):
        self.x = x
        self.y = y
        
    def __repr__(self):
        return f"x:{self.x} y:{self.y}"
        
    def __str__(self):
        return f"x:{self.x} y:{self.y}"
    
def plot_vecs(cH): 
    for i in range (len(cH)):
        plt.plot(cH[i].x, cH[i].y, 'bo', linestyle="--")

    
def leftmost(points):
    minim = 0
    for i in range(1,len(points)):
        if points[i][0] < points[minim][0]:
            minim = i
        elif points[i][0] == points[minim][0]:
            if points[i][1] > points[minim][1]:
                minim = i
    return minim

def det(p1, p2, p3):
    """ 
    > 0: CCW turn
    < 0 CW turn
    = 0: colinear
    """
    return (p2[0] - p1[0]) * (p3[1] - p1[1]) \
        -(p2[1] - p1[1]) * (p3[0] - p1[0])

def rotation(punto1, punto2):
    ans = punto1.x*punto2.y - punto1.y*punto2.x
    if ans > 0:
        return -1
    if ans < 0:
        return 1
    return 0

def turn(punto1, punto2, punto3):
    v1 = Punto(punto2.x - punto1.x, punto2.y - punto1.y)
    v2 = Punto(punto3.x - punto1.x, punto3.y - punto1.y)
    rot = rotation(v1, v2)
    
    print(rot)
    return rot

points = [(2,1),(1,2),(2,2),(5,4),(3,4),(1,2)]

hull = []
l = leftmost(points)
leftMost = points[l]
currentVertex = leftMost
hull.append(currentVertex)
nextVertex = points[1]
index = 2
nextIndex = -1
while True:
    c0 = currentVertex
    c1 = nextVertex

    checking = points[index]
    c2 = checking

    crossProduct = det(currentVertex, nextVertex, checking)
    if crossProduct < 0:
        nextVertex = checking
        nextIndex = index
    index += 1
    if index == len(points):
        if nextVertex == leftMost:
            break
        index = 0
        hull.append(nextVertex)
        currentVertex = nextVertex
        nextVertex = leftMost
print(hull)

def create_points(ct, min = 0, max = 50):
    return [[randint(min, max), randint(min, max)]\
        for _ in range(ct)] 

def scatter_plot(coords, convex_hull = None):
    xs, ys = zip(*coords) #unzip into x and y coordinates
    plt.scatter(xs, ys)

    if convex_hull:
        for i in range(1, len(convex_hull) + 1):
            if i == len(convex_hull): i = 0 #wrap
            c0 = convex_hull[i-1]
            c1 = convex_hull[i]
            plt.plot((c0[0], c1[0]), (c0[1], c1[1]), 'r')
    plt.show()
    
scatter_plot(points, hull)