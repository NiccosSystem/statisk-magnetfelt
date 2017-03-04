import numpy as np
from scipy.misc import derivative
from matplotlib import pyplot as plt
import math

#partiellderivert
def P_D(func, index, points):
    points_copy = points[:] #Copy points as not to alter the points
    def as_func_of(x):      #Need to input the function as a one variable function
        points_copy[index] = x #The variable is the index'th element
        return func(*points_copy) #Now it is only a function of points[index]. Need star to input a list as arguments
    return derivative(as_func_of, points[index], dx = 1e-6) #See manual to scipy.misc.derivative

#leser inn filer i numpy array
def get_values(file):
    values = []
    data = open(file, 'r')
    for line in data:
        if line != "":
            values.append(float(line))
    return np.array(values)

#Regner ut gauss for alle verdier og alle variabler kombinert
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

#helmholtsspolefunksjon
def ar(x, R, I, mu, N, a):
    return ((N*mu*I)/(2*R))*((1+((x-a/2)/R)**2)**(-3/2)+(1+((x+a/2)/R)**2)**(-3/2))

#verdiene for de forskjellige a-ene
x_values = get_values("2-x.txt")
b2r_values = get_values("2-r-b2.txt")
b22r_values = get_values("2-2r-b2.txt")
b205r_values = get_values("2-05r-b2.txt")

#konstante verdier
k_mu = 1.26 * 10**(-6)
k_N = 330
k_I = 1
k_R = 0.07

#x-verdier pluss konstante verdier
constants = [x_values, k_R, k_I, k_mu, k_N]
print([]+constants+[5])

#forskjellige a-verdiene
a_r = k_R
a_2r = 2*k_R
a_05r = 0.5*k_R

#usikkerhetene
d_mu = 0
d_N = 0
d_I = 0.05
d_R = 0.0001
d_x = 0.0005
d_a = 0.001
uncertainties = [d_x, d_R, d_I, d_mu, d_N, d_a]

#delx_value = P_D(stromsloyfe, 0, [x_values, k_R, k_I, k_mu, k_N])
#delI_value = P_D(stromsloyfe, 2, [x_values, k_R, k_I, k_mu, k_N])
#delR_value = P_D(stromsloyfe, 1, [x_values, k_R, k_I, k_mu, k_N])
#print(delx_value)

#gauss for de forskjellige a-ene
gauss_r = gauss(x_values, b2r_values, []+constants+[a_r], uncertainties, ar)
gauss_2r = gauss(x_values, b22r_values, []+constants+[a_2r], uncertainties, ar)
gauss_05r = gauss(x_values, b205r_values, []+constants+[a_05r], uncertainties, ar)

print(gauss_r)
print(gauss_2r)
print(gauss_05r)

#gaussene ganget med måleverdier, NB! skal egentlig ganges med teoretiske verdier etc
c_uncr = np.multiply(b2r_values, np.array(gauss_r))
c_unc2r = np.multiply(b22r_values, np.array(gauss_2r))
c_unc05r = np.multiply(b205r_values, np.array(gauss_05r))


x_space = np.linspace(-0.2, 0.2, 200)

#funksjonene for de forskjellige a-verdiene
fr = ((k_N*k_mu*k_I)/(2*k_R))*((1+((x_space-k_R/2)/k_R)**2)**(-3/2)+(1+((x_space+k_R/2)/k_R)**2)**(-3/2)) * 10**4
f2r = ((k_N*k_mu*k_I)/(2*k_R))*((1+((x_space-k_R)/k_R)**2)**(-3/2)+(1+((x_space+k_R)/k_R)**2)**(-3/2)) * 10**4
f05r = ((k_N*k_mu*k_I)/(2*k_R))*((1+((x_space-k_R/4)/k_R)**2)**(-3/2)+(1+((x_space+k_R/4)/k_R)**2)**(-3/2)) * 10**4

print(c_uncr)
#plotte shit
plt.plot(x_space, fr, '-', color='0', label=r'$B(a=R)$')
plt.plot(x_space, f2r, '-', color='0.33', label=r'$B(a=2R)$')
plt.plot(x_space, f05r, '-', color='0.66', label=r'$B(a=R/2)$')
plt.errorbar(x_values, b2r_values, yerr=np.array(c_uncr), linestyle="None", marker='.', color='0', label="Målepunkter")
plt.errorbar(x_values, b22r_values, yerr=np.array(c_unc2r), linestyle="None", marker='.', color='0.33')
plt.errorbar(x_values, b205r_values, yerr=np.array(c_unc05r), linestyle="None", marker='.', color='0.66')
plt.xlabel(r"$x$ (m)", size=18)
plt.ylabel(r"$B(x)$ (Gauss)", size=18)
plt.legend()
plt.show()


