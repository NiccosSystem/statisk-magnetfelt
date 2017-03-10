import numpy as np
from matplotlib import gridspec
from scipy.misc import derivative
from matplotlib import pyplot as plt
import math

#partiellderivert-tager
def P_D(func, index, points):
    points_copy = points[:] #Copy points as not to alter the points
    def as_func_of(x):      #Need to input the function as a one variable function
        points_copy[index] = x #The variable is the index'th element
        return func(*points_copy) #Now it is only a function of points[index]. Need star to input a list as arguments
    return derivative(as_func_of, points[index], dx = 1e-6) #See manual to scipy.misc.derivative

#tekstfil-verdier inn i numpy array
def get_values(file):
    values = []
    data = open(file, 'r')
    for line in data:
        if line != "":
            values.append(float(line))
    return np.array(values)

#gauss for gitte verdier
def gauss(x_val, b2_val, variables, uncertains, func):
    unc = []
    assigned = False
    reci = np.reciprocal(b2_val)
    i = 0
    while i < len(variables):
        partial = np.array(P_D(func, i, variables))
        print("partial {0}".format(partial))
        step1 = partial*uncertains[i]
        c_unc = (10**4 * np.multiply(step1, reci))**2
        i += 1
        if assigned:
            unc += c_unc
            continue
        unc = c_unc
        assigned = True
    return np.sqrt(unc)

#funksjon for solenoide
def solf(z, R, I, mu, N, l):
    return ((N*mu*I)/(2*l))*((z/((z**2+R**2)**(1/2)))+((l-z)/(((l-z)**2 + R**2)**(1/2))))

#x- og b-verdier
x_values = get_values("3-x.txt")
b2_values = get_values("3-b2.txt")

#konstanter
k_mu = 1.26 * 10**(-6)
k_N = 368
k_I = 1
k_R = 0.05
k_l = 0.397

#x-verdier pluss konstanter inn i array
#constants = [x_values, k_R, k_I, k_mu, k_N, k_l]

#usikkerheter inni array
d_mu = 0
d_N = 0
d_I = 0.05
d_R = 0.0001
d_x = 0.0005
d_l = 0.0001
uncertainties = [d_x, d_R, d_I, d_mu, d_N, d_l]

#gauss for solenoide
#gauss_s = gauss(x_values, b2_values, constants, uncertainties, solf)

#print(gauss_s)

#gange gauss med funksjonsverdier NB! skal egentlig ganges med teoretiske verdier etc
#c_unc = np.multiply(b2_values, np.array(gauss_s))

#linspace for teoretiske verdier
z = np.linspace(-0.10, 0.60, 500)

#funksjon for solenoide

sf = 10**4 * ((k_N*k_mu*k_I)/(2*k_l))*((z/((z**2+k_R**2)**(1/2)))+((k_l-z)/(((k_l-z)**2 + k_R**2)**(1/2))))
#print(len(z), len(sf))
constants = [z, k_R, k_I, k_mu, k_N, k_l]
gauss_teoretisk = gauss(z, sf, constants, uncertainties, solf)
c_unc = np.multiply(sf, np.array(gauss_teoretisk))

avvik = b2_values - (10**4 * ((k_N*k_mu*k_I)/(2*k_l))*((x_values/((x_values**2+k_R**2)**(1/2)))+((k_l-x_values)/(((k_l-x_values)**2 + k_R**2)**(1/2)))))
#plotte shit
#print(c_unc)

#f, axs = plt.subplots(2, sharex=True)

#axs[0].plot(z, sf, '-', color='black', label=r"$B$")
#axs[0].scatter(x_values, b2_values, marker='o', color='black', label='Målepunkter')
#plt.xlabel(r"$x$ (m)", size=18)
#plt.ylabel(r"$B(x)$ (Gauss)", size=18)
#f.add_subplot(111, frameon=False)
#plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
#plt.legend()

xlim = [min(x_values)*1.05, max(x_values)*1.05]

fur = plt.figure(1)
gsr = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
ax1r = plt.subplot(gsr[0])
ax2r = plt.subplot(gsr[1], sharex=ax1r)
ax1r.plot(z, sf, '-', color='0', label=r'$B(a=R)$')
ax1r.scatter(x_values, b2_values, marker='.', color='0', label="Målepunkter")
ax2r.plot(z, np.abs(c_unc), '-', color="0")
ax2r.plot(z, -np.abs(c_unc), '-', color="0")
ax2r.scatter(x_values, avvik, marker='.', color="0")
fur.text(0.04, 0.5, r"$B(x)$ (Gauss)", va='center', rotation='vertical', size=18)
ax1r.legend()
plt.xlim(xlim)
plt.xlabel(r"$x$ (m)", size=18)
#plt.tight_layout()

#plt.subplot()
plt.show()





