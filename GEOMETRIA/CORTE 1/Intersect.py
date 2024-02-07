import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True

class Punto:
  def __init__(self, x , y):
    self.x = x
    self.y = y

class Vector:
  def __init__(self, inicial , final):
    self.inicial = inicial
    self.final = final

def plot_vecs(vec1, vec2): 
    point1i = vec1.inicial
    point1f = vec1.final
    
    point2i = vec2.inicial
    point2f = vec2.final

    x1_values = [point1i.x, point1f.x]
    y1_values = [point1i.y, point1f.y]
    
    x2_values = [point2i.x, point2f.x]
    y2_values = [point2i.y, point2f.y]

    plt.plot(x1_values, y1_values, 'bo', linestyle="--")
    plt.plot(x2_values, y2_values, 'bo', linestyle="--")

def plot_rot(p1, p2):
    x1_values = [p1.x, 0]
    y1_values = [p1.y, 0]
    
    x2_values = [p2.x, 0]
    y2_values = [p2.y, 0]

    plt.plot(x1_values, y1_values, 'bo', linestyle="--")
    plt.plot(x2_values, y2_values, 'bo', linestyle="--")

    
def plot_turn(p1, p2, p3):
    x1_values = [p1.x, p2.x]
    y1_values = [p1.y, p2.y]
    
    x2_values = [p2.x, p3.x]
    y2_values = [p2.y, p3.y]

    plt.plot(x1_values, y1_values, 'bo', linestyle="--")
    plt.plot(x2_values, y2_values, 'bo', linestyle="--",color='red')

    

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

def intersect(v1, v2):
    p1 = v1.inicial
    p2 = v1.final
    
    p3 = v2.inicial
    p4 = v2.final
    
    if turn(p1, p2, p3) == turn(p1, p2, p4) != 0:
        return False
    
    if turn(p1, p2, p3) == - turn(p1, p2, p4) != 0:
        if turn(p3, p4, p1) != turn(p3, p4, p2):
            return True
        return False
    
    if turn(p1, p2, p3) == 0:
        x1 = p1.x
        x2 = p2.x
        x3 = p3.x
        if min(x1, x2) <= x3 <= max(x1, x2):
            return True
        if turn(p1, p2, p4) == 0:
            x4 = p4[0]
            if min(x1, x2) <= x4 <= max(x1, x2):
                return True
            return False
        return False
    
    if turn(p1, p2, p4) == 0:
        x1 = p1.x
        x2 = p2.x
        x4 = p4.x
        if min(x1, x2) <= x4 <= max(x1, x2):
            return True
        if turn(p1, p2, p3) == 0:
            x3 = p3[0]
            if min(x1, x2) <= x3 <= max(x1, x2):
                return True
            return False
        return False

p1 = Punto(2,2)
p2 = Punto(4,6)
p3 = Punto(1,2)
p4 = Punto(5,6)

v1 = Vector(p1, p2)
v2 = Vector(p3, p4)



print(intersect(v1, v2))
plot_vecs(v1, v2)


plt.xlim(0, 10)
plt.ylim(0, 10)
plt.show()