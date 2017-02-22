import numpy as np
from scipy.misc import derivative
from matplotlib import pyplot as plt
import math

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

def ar(x, R, I, mu, N, a):
    return ((N*mu*I)/(2*R))*((1+((x-a/2)/R)**2)**(-3/2)+(1+((x+a/2)/R)**2)**(-3/2))

x_values = get_values("2-x.txt")
b2r_values = get_values("2-r-b2.txt")
b22r_values = get_values("2-2r-b2.txt")
b205r_values = get_values("2-05r-b2.txt")

k_mu = 1.26 * 10**(-6)
k_N = 330
k_I = 1
k_R = 0.07
constants = [x_values, k_R, k_I, k_mu, k_N]
print([]+constants+[5])

a_r = k_R
a_2r = 2*k_R
a_05r = 0.5*k_R

d_mu = 0.01 * 10**(-6)
d_N = 0
d_I = 0.01
d_R = 0.0001
d_x = 0.0005
d_a = 0.001
uncertainties = [d_x, d_R, d_I, d_mu, d_N, d_a]

#delx_value = P_D(stromsloyfe, 0, [x_values, k_R, k_I, k_mu, k_N])
#delI_value = P_D(stromsloyfe, 2, [x_values, k_R, k_I, k_mu, k_N])
#delR_value = P_D(stromsloyfe, 1, [x_values, k_R, k_I, k_mu, k_N])
#print(delx_value)

gauss_r = gauss(x_values, b2r_values, []+constants+[a_r], uncertainties, ar)
gauss_2r = gauss(x_values, b22r_values, []+constants+[a_2r], uncertainties, ar)
gauss_05r = gauss(x_values, b205r_values, []+constants+[a_05r], uncertainties, ar)

print(gauss_r)
print(gauss_2r)
print(gauss_05r)



        #  uncx = delx_value[i] * d_x * b2_val[i] ** (-1)
        #  uncI = delI_value[i] * d_I * b2_val[i] ** (-1)
        #  uncR = delR_value[i] * d_R * b2_val[i] ** (-1)
        #  squared = uncx ** 2 + uncI ** 2 + uncR ** 2
        #  unc.append(10 ** 4 * math.sqrt(squared))
        # i += 1

b2r_unitless = np.array(b2r_values)/b2r_values[11]
b22r_unitless = np.array(b22r_values)/b22r_values[11]
b205r_unitless = np.array(b205r_values)/b205r_values[11]

x_unitless = np.divide(np.array(x_values), k_R)

c_uncr = np.multiply(b2r_unitless, np.array(gauss_r))
c_unc2r = np.multiply(b22r_unitless, np.array(gauss_2r))
c_unc05r = np.multiply(b205r_unitless, np.array(gauss_05r))


x_space = np.linspace(-3, 3, 200)

fr = (1/2)*((1+(x_space-(a_r/(2*k_R)))**2)**(-3/2)+(1+(x_space+(a_r/(2*k_R)))**2)**(-3/2)) * ((1+(a_r/(2*k_R))**2)**(3/2))
f2r = (1/2)*((1+(x_space-(a_2r/(2*k_R)))**2)**(-3/2)+(1+(x_space+(a_2r/(2*k_R)))**2)**(-3/2)) * ((1+(a_2r/(2*k_R))**2)**(3/2))
f05r = (1/2)*((1+(x_space-(a_05r/(2*k_R)))**2)**(-3/2)+(1+(x_space+(a_05r/(2*k_R)))**2)**(-3/2)) * ((1+(a_05r/(2*k_R))**2)**(3/2))

print(c_uncr)
#print(current_unc)
#plt.plot(x_space, (1+x_space**2)**(-3/2), '-', color='black', label=r"$B$")
plt.plot(x_space, fr, '-')
plt.plot(x_space, f2r, '-')
plt.plot(x_space, f05r, '-')
plt.errorbar(x_unitless, b2r_unitless, yerr=np.array(c_uncr), linestyle="None", marker='.', color='black', label='Målepunkter')
plt.errorbar(x_unitless, b22r_unitless, yerr=np.array(c_unc2r), linestyle="None", marker='.', color='blue', label='Målepunkter')
plt.errorbar(x_unitless, b205r_unitless, yerr=np.array(c_unc05r), linestyle="None", marker='.', color='red', label='Målepunkter')
plt.xlabel(r"$x/R$")
plt.ylabel(r"$B(x/R)/B(0)$")
plt.legend()
plt.show()



# Take the partial derivative of the function with respect to x1, and input the values in the list [A,B,C]. Since Python is 0-indexed (lists/arrays start at 0),
#  the x1 is the 0th variable in the function. The first delx1_value gives the partial derivative in one point.
# The second, delx1_list, gives the partial derivative as an array of 100 elements where x1 varies from -2 to 2.
#delx1_value = P_D(function_f, 0, [x1_value, x2_constant, x3_constant])
#delx1_list = P_D(function_f, 0, [x1_list, x2_constant, x3_constant])
#delx2 = P_D(function_f, 1, [x1_list, x2_constant, x3_constant])
#if delx1_value == delx1_list[74]:
#    print('This works') #the list should be the same as the value since we chose number 74 as an example.

#plot the result as a function of x1_list
#plt.plot(x1_list, (x3_constant*(1+x1_list**2/x2_constant**2)**(-3.0/2)))
#plt.plot(x1_list, delx1_list)
#plt.plot(x1_list, delx2)
#plt.show()



