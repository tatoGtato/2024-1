import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True


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
    
def leftmost(points):
    minim = 0
    for i in range(1,len(points)):
        if points[i].x < points[minim].x:
            minim = i
        elif points[i].x == points[minim].x:
            if points[i].y > points[minim].y:
                minim = i
    return minim

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



def convexHull(p):
    #Sort the points by x-coordinate
    p.sort(key=lambda punto: punto.x)
    #Create lsuper with the fist two values of p
    lSuper = [p[0], p[1]]
    lLower = []
    L = []
    for i in range (3, len(p)):
        lSuper.append(p[i])
        #While lSuper contains more than 2 points and 
        #the las three points dont make a right turn
        while (len(lSuper) > 2 and turn(lSuper[-3:][0],lSuper[-3:][1],lSuper[-3:][2]) != 1):
            #Remove the middle point of the last three
            lSuper.remove(lSuper[-3:][1])
    #Put in the last two point into Llower
    lLower.append(p[-2:][0])
    lLower.append(p[-2:][1])

    for i in reversed(range(len(p)-2, 1)):
        lLower.append(p[i])
        
        while (len(lLower) > 2 and turn(lLower[-3:][0],lLower[-3:][1],lLower[-3:][2]) != 1):
            lLower.remove(lLower[-3:][1])
               
    lLower.remove(lLower[0])  
    lLower.remove(lLower[-1]) 
    
    L = lSuper + lLower
    return L

data = np.array([
    [2, 3],
    [4, 6],
    [1, 2],
    [2,1]
])
x, y = data.T
plt.scatter(x,y)
plt.show()

print(convexHull(P))
P = [Punto(2,3), Punto(4,6), Punto(1,2), Punto(2,1)]


# plt.xlim(0, 10)
# plt.ylim(0, 10)
# plt.show()