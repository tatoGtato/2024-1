import math as m
import numpy as np
from numpy.linalg import inv

A = (
    [1,0,0,0],
    [1,4,1,0],
    [0,1,4,1],
    [0,0,0,1],
)

b= (
    [0],
    [3*(m.exp(2)-m.exp(1))-3*(m.exp(1)-1)],
    [3*(m.exp(3)-m.exp(2))-3*(m.exp(2)-m.exp(1))],
    [0],
)

 
Ainv =  inv(A)

# Las matriz X es la matriz de c's 
X = np.dot(Ainv,b)
 
 
# print resulted matrix
print(res)