import numpy as np
from scipy.misc import derivative
from matplotlib import pyplot as plt
import math

#Tar partiellderivert av funksjon
def P_D(func, index, points):
    points_copy = points[:] #Copy points as not to alter the points
    def as_func_of(x):      #Need to input the function as a one variable function
        points_copy[index] = x #The variable is the index'th element
        return func(*points_copy) #Now it is only a function of points[index]. Need star to input a list as arguments
    return derivative(as_func_of, points[index], dx = 1e-6) #See manual to scipy.misc.derivative

#Henter inn data fra fil og legger det i en numpy-array
def get_values(file):
    values = []
    data = open(file, 'r')
    for line in data:
        if line != "":
            values.append(float(line))
    return np.array(values)

#Funksjonen for en strømsløyfe
def stromsloyfe(x, R, I, mu, N):
    return ((N*mu*I)/(2*R))*(1+(x/R)**2)**(-3/2)

#x-verdier i meter og b2-verdier
x_values = get_values("1-x.txt")
b2_values = get_values("1-b2.txt")

#Konstante verdier for mu, N, I og R
k_mu = 1.26 * 10**(-6)
k_N = 330
k_I = 1
k_R = 0.07

#Usikkerheten i de forskjellige variablene
d_mu = 0.01 * 10**(-6)
d_N = 0
d_I = 0.05
d_R = 0.0001
d_x = 0.0005

#De partiellderiverte
delx_value = P_D(stromsloyfe, 0, [x_values, k_R, k_I, k_mu, k_N])
delI_value = P_D(stromsloyfe, 2, [x_values, k_R, k_I, k_mu, k_N])
delR_value = P_D(stromsloyfe, 1, [x_values, k_R, k_I, k_mu, k_N])
print(delx_value)

#Her lages gauss-verdier for alle variablene med henhold til måleverdier NB!! SKAL EGENTLIG VÆRE TEORETISKE VERDIER, IKKE MÅLEVERDIER!!
unc = []
i = 0
while i < len(x_values):
    squared = (delx_value[i] * d_x * (10**(-4)*b2_values[i]**(-1)))**2 + (delI_value[i] * d_I * (10**(-4)*b2_values[i])**(-1))**2 + (delR_value[i] * d_R * (10**(-4)*b2_values[i])**(-1))**2
    unc.append(math.sqrt(squared))
    i += 1

#Multipliser elemtnvis alle gauss-utregninger med måleverdi for å få absolutt usikkerhet NB! SE OVENFOR
current_unc = np.multiply(b2_values, unc)

#Lage linspace for teoretiske verdier
x_space = np.linspace(-0.2, 0.2, 500)

#Plotte shit
print(unc)
print(current_unc)
plt.plot(x_space, 10**4 * ((k_N*k_mu*k_I)/(2*k_R))*(1+(x_space/k_R)**2)**(-3/2), '-', color='black', label=r"$B$")
plt.errorbar(x_values, b2_values, yerr=np.array(current_unc), linestyle="None", marker='.', color='black', label='Målepunkter')
plt.xlabel(r"$x$ (m)", size=18)
plt.ylabel(r"$B(x)$ (Gauss)", size=18)
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



