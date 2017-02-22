import numpy as np
from scipy.misc import derivative
from matplotlib import pyplot as plt
import math

#test
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
        partial = np.array(P_D(func, i, variables))
        step1 = partial*uncertains[i]
        c_unc = (10**4 * np.multiply(step1, reci))**2
        i += 1
        if assigned:
            unc += c_unc
            continue
        unc = c_unc
        assigned = True
    return np.sqrt(unc)

def solf(z, R, I, mu, N, l):
    return ((N*mu*I)/(2*l))*((z/((z**2+R**2)**(1/2)))+((l-z)/(((l-z)**2 + R**2)**(1/2))))

x_values = get_values("3-x.txt")
b2_values = get_values("3-b2.txt")

k_mu = 1.26 * 10**(-6)
k_N = 368
k_I = 1
k_R = 0.05
k_l = 0.397
constants = [x_values, k_R, k_I, k_mu, k_N, k_l]

d_mu = 0.01 * 10**(-6)
d_N = 0
d_I = 0.01
d_R = 0.0001
d_x = 0.0005
d_l = 0.0001
uncertainties = [d_x, d_R, d_I, d_mu, d_N, d_l]

#delx_value = P_D(stromsloyfe, 0, [x_values, k_R, k_I, k_mu, k_N])
#delI_value = P_D(stromsloyfe, 2, [x_values, k_R, k_I, k_mu, k_N])
#delR_value = P_D(stromsloyfe, 1, [x_values, k_R, k_I, k_mu, k_N])
#print(delx_value)

gauss_s = gauss(x_values, b2_values, constants, uncertainties, solf)

print(gauss_s)



        #  uncx = delx_value[i] * d_x * b2_val[i] ** (-1)
        #  uncI = delI_value[i] * d_I * b2_val[i] ** (-1)
        #  uncR = delR_value[i] * d_R * b2_val[i] ** (-1)
        #  squared = uncx ** 2 + uncI ** 2 + uncR ** 2
        #  unc.append(10 ** 4 * math.sqrt(squared))
        # i += 1

b2_unitless = np.array(b2_values)/b2_values[7]

x_unitless = np.divide(np.array(x_values), k_R)

c_unc = np.multiply(b2_unitless, np.array(gauss_s))


z = np.linspace(-2, 12, 500)

sf = ((z/((1+z**2)**(1/2))) + (k_l-z*k_R)/(((k_l-z*k_R)**2 + k_R**2)**(1/2))) / (k_l/((k_l**2 + k_R**2)**(1/2)))

print(c_unc)
#print(current_unc)
#plt.plot(x_space, (1+x_space**2)**(-3/2), '-', color='black', label=r"$B$")
plt.plot(z, sf, '-')
plt.errorbar(x_unitless, b2_unitless, yerr=np.array(c_unc), linestyle="None", marker='.', color='red', label='Målepunkter')
plt.xlabel(r"$x/R$")
plt.ylabel(r"$B(x/R)/B(0)$")
plt.legend()
plt.show()





