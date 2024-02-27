import matplotlib.pyplot as plt
import numpy as np
import random
import time


""" 
CLASSES
""" 


class Punto:
    def __init__(self, x , y):
        self.x = x
        self.y = y
        
    def __repr__(self):
        return f"(x:{self.x},y:{self.y})"
        
    def __str__(self):
        return f"(x:{self.x},y:{self.y})"
    
    def __lt__(self, other):
        if (self.x < other.x):
            return self.x < other.x
        elif (self.x == other.x):
            return self.y < other.y
    
    def __gt__(self, other):
        if (self.x > other.x):
            return self.x > other.x
        elif(self.x == other.x):
            return self.y > other.y

""" 
CLASSES
""" 
 
""" 
FUNCTIONS
""" 
    
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
    return rot

def convexHull(points):
    hull = []  
    points.sort()
    print(points)
    lUpper = []
    lLower = []
    n = len(points)
    
    #PUT THE POINTS P1 & P2 IN LUPPER
    lUpper.append(points[0])
    lUpper.append(points[1])
    
    for i in range (2, n):
        lUpper.append(points[i])
        #LUPPER CONTAINS MORE THAN 2 AND LAST 3 DONT MAKE A RIGHT TURN
        while(len(lUpper) > 2 and turn(lUpper[-3:][0],lUpper[-3:][1],lUpper[-3:][2]) != 1):
            lUpper.remove(lUpper[-3:][1])
    
    #PUT THE POINTS Pn & Pn-1 IN Llower
    lLower.append(points[-2:][1])
    lLower.append(points[-2:][0])
    
    for i in range (n-1, 0-1, -1):
        lLower.append(points[i])
        
        #LLOWER CONTAINS MORE THAN 2 AND LAST 3 DONT MAKE A RIGHT TURN
        while(len(lLower) > 2 and turn(lLower[-3:][0],lLower[-3:][1],lLower[-3:][2]) != 1):
            lLower.remove(lLower[-3:][1])
    
    #REMOVE THE FIRST AND LAST OF LLOWER TO AVOID DUPLICATION
    lLower.remove(lLower[0])
    lLower.remove(lLower[-1])
    
    print("UPPER: " + str(lUpper))
    print("LOWER: " + str(lLower))
    
    hull = lUpper + lLower
    
    return hull

def plotPoints(points, hull):
    dataL = []
    fig, (ax1, ax2) = plt.subplots(2)
    
    #PLOT 1
    for i in points:
        dataL.append([i.x, i.y])
    
    data = np.array(dataL)  
    x, y = data.T
    ax1.scatter(x,y)
    
    ax1.set_xlim(left=-1, right=102)
    ax1.set_ylim(bottom=-1, top=102)
    
    #PLOT 2
    ax2.scatter(x,y)
    x_values = []
    y_values = []
    
    for i in hull:
        x_values.append(i.x)
        y_values.append(i.y)
        
    x_values.append((hull[0]).x)
    y_values.append((hull[0]).y)    
    
    ax2.plot(x_values, y_values, 'bo', linestyle="--")
    
    ax2.set_xlim(left=-1, right=102)
    ax2.set_ylim(bottom=-1, top=102)

def generateRandomPoints(n):
    points = []
    
    for i in range (0, n):
        x = random.randint(0,10)
        y = random.randint(0,10)
        points.append(Punto(x,y))
        
    return points

""" 
FUNCTIONS
"""

"""
TESTS
"""

def generate_test_case():
    """Genera un caso de prueba con un conjunto de 9 puntos aleatorios."""
    points = [Punto(random.randint(0, 100), random.randint(0, 100)) for _ in range(9)]
    return points

def run_test_case(points):
    """Ejecuta el algoritmo convexHull en un conjunto de puntos y mide el tiempo."""
    start_time = time.time()
    convex_hull = convexHull(points)
    end_time = time.time()
    execution_time = end_time - start_time
    return convex_hull, execution_time

def print_test_case_result(test_case_number, points, convex_hull, execution_time):
    """Imprime el resultado de un caso de prueba."""
    print(f"Test Case {test_case_number}:")
    print("Input:")
    for point in points:
        print(point)
    print("Output (Convex Hull):")
    for point in convex_hull:
        print(point)
    print(f"Execution Time: {execution_time} seconds")
    print()

# Ejecutar las pruebas
for i in range(6):
    points = generate_test_case()
    convex_hull, execution_time = run_test_case(points)
    print_test_case_result(i+1, points, convex_hull, execution_time)
    plotPoints(points, convex_hull)


# P = generateRandomPoints(9)
# print(P)
# hull = convexHull(P)

# print("HULL: " + str(hull))



plt.show()