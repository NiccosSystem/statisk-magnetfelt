import math
import numpy as np
from scipy.misc import derivative
from matplotlib import pyplot as plt

def P_D(func, index, points):
    points_copy = points[:] #Copy points as not to alter the points
    def as_func_of(x):      #Need to input the function as a one variable function
        points_copy[index] = x #The variable is the index'th element
        return func(*points_copy) #Now it is only a function of points[index]. Need star to input a list as arguments
    return derivative(as_func_of, points[index], dx = 1e-6) #See manual to scipy.misc.derivative

def get_values(file):
    values = []
    data = open(file, 'r')
    for line in data:
        if line != "":
            values.append(float(line))
    return np.array(values)

def gauss(x_val, b2_val, variables, uncertains, func):
    unc = []
    assigned = False
    reci = np.reciprocal(b2_val)
    i = 0
    while i < len(variables):
        partial = np.array(P_D(stromsloyfe, i, variables))
        step1 = partial*uncertains[i]
        c_unc = (10**4 * np.multiply(step1, reci))**2
        i += 1
        if assigned:
            unc += c_unc
            continue
        unc = c_unc
        assigned = True
    return np.sqrt(unc)


def stromsloyfe(x, R, I, mu, N):
    return ((N*mu*I)/(2*R))*(1+(x/R)**2)**(-3/2)

x_values = np.array(get_values("1-x.txt"))
b2_values = np.array(get_values("1-b2.txt"))

k_mu = 1.26 * 10**(-6)
k_N = 330
k_I = 1
k_R = 0.07

d_mu = 0.01 * 10**(-6)
d_N = 0
d_I = 0.01
d_R = 0.0001
d_x = 0.0005

test = gauss(x_values, b2_values, [x_values, k_R, k_I, k_mu, k_N], [d_x, d_R, d_I, d_mu, d_N], stromsloyfe)
test = np.multiply(test, np.divide(b2_values, b2_values[9]))
print(test)