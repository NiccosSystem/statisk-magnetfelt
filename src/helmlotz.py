import numpy as np
from matplotlib import gridspec
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

#print([]+constants+[5])

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



#print(gauss_r)
#print(gauss_2r)
#print(gauss_05r)

#gaussene ganget med måleverdier, NB! skal egentlig ganges med teoretiske verdier etc

#x_space = np.linspace(-0.2, 0.2, 200)
x_space = np.linspace(min(x_values), max(x_values), 200)

#funksjonene for de forskjellige a-verdiene
fr = ((k_N*k_mu*k_I)/(2*k_R))*((1+((x_space-k_R/2)/k_R)**2)**(-3/2)+(1+((x_space+k_R/2)/k_R)**2)**(-3/2)) * 10**4
f2r = ((k_N*k_mu*k_I)/(2*k_R))*((1+((x_space-k_R)/k_R)**2)**(-3/2)+(1+((x_space+k_R)/k_R)**2)**(-3/2)) * 10**4
f05r = ((k_N*k_mu*k_I)/(2*k_R))*((1+((x_space-k_R/4)/k_R)**2)**(-3/2)+(1+((x_space+k_R/4)/k_R)**2)**(-3/2)) * 10**4

constants = [x_space, k_R, k_I, k_mu, k_N]

#gauss for de forskjellige a-ene
gauss_r = gauss(x_values, fr, []+constants+[a_r], uncertainties, ar)
gauss_2r = gauss(x_values, f2r, []+constants+[a_2r], uncertainties, ar)
gauss_05r = gauss(x_values, f05r, []+constants+[a_05r], uncertainties, ar)

c_uncr = np.multiply(fr, np.array(gauss_r))
c_unc2r = np.multiply(f2r, np.array(gauss_2r))
c_unc05r = np.multiply(f05r, np.array(gauss_05r))

avvikr = b2r_values - (((k_N*k_mu*k_I)/(2*k_R))*((1+((x_values-k_R/2)/k_R)**2)**(-3/2)+(1+((x_values+k_R/2)/k_R)**2)**(-3/2)) * 10**4)
avvik2r = b22r_values - (((k_N*k_mu*k_I)/(2*k_R))*((1+((x_values-k_R)/k_R)**2)**(-3/2)+(1+((x_values+k_R)/k_R)**2)**(-3/2)) * 10**4)
avvik05r = b205r_values - (((k_N*k_mu*k_I)/(2*k_R))*((1+((x_values-k_R/4)/k_R)**2)**(-3/2)+(1+((x_values+k_R/4)/k_R)**2)**(-3/2)) * 10**4)

#print(c_uncr)
#plotte shit

xlim = [min(x_values)*1.05, max(x_values)*1.05]

fur = plt.figure(1)
gsr = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
ax1r = plt.subplot(gsr[0])
ax2r = plt.subplot(gsr[1], sharex=ax1r)
ax1r.plot(x_space, fr, '-', color='0', label=r'$B(a=R)$')
ax1r.scatter(x_values, b2r_values, marker='.', color='0', label="Målepunkter")
ax2r.plot(x_space, np.abs(c_uncr), '-', color="0")
ax2r.plot(x_space, -np.abs(c_uncr), '-', color="0")
ax2r.scatter(x_values, avvikr, marker='.', color="0")
fur.text(0.04, 0.5, r"$B(x)$ (Gauss)", va='center', rotation='vertical', size=18)
ax1r.legend()
plt.xlim(xlim)
plt.xlabel(r"$x$ (m)", size=18)
#plt.tight_layout()

fu2r = plt.figure(2)
gs2r = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
ax12r = plt.subplot(gs2r[0])
ax22r = plt.subplot(gs2r[1], sharex=ax12r)
ax12r.plot(x_space, f2r, '-', color='0', label=r'$B(a=2R)$')
ax12r.scatter(x_values, b22r_values, marker='.', color='0', label="Målepunkter")
ax22r.plot(x_space, np.abs(c_unc2r), '-', color="0")
ax22r.plot(x_space, -np.abs(c_unc2r), '-', color="0")
ax22r.scatter(x_values, avvik2r, marker='.', color="0")
fu2r.text(0.04, 0.5, r"$B(x)$ (Gauss)", va='center', rotation='vertical', size=18)
ax12r.legend()
plt.xlim(xlim)
plt.xlabel(r"$x$ (m)", size=18)
#plt.tight_layout()

fu05r = plt.figure(3)
gs05r = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
ax105r = plt.subplot(gs05r[0])
ax205r = plt.subplot(gs05r[1], sharex=ax105r)
tline05r = ax105r.plot(x_space, f05r, '-', color='0', label=r'$B(a=R/2)$')
sctr05r = ax105r.scatter(x_values, b205r_values, marker='.', color='0', label="Målepunkter")
ax205r.plot(x_space, np.abs(c_unc05r), '-', color="0")
ax205r.plot(x_space, -np.abs(c_unc05r), '-', color="0")
ax205r.scatter(x_values, avvik05r, marker='.', color="0")
fu05r.text(0.04, 0.5, r"$B(x)$ (Gauss)", va='center', rotation='vertical', size=18)
ax105r.legend()
plt.xlabel(r"$x$ (m)", size=18)
plt.xlim(xlim)
#plt.tight_layout()

#plt.plot(x_space, fr, '-', color='0', label=r'$B(a=R)$')
#plt.plot(x_space, f2r, '-', color='0.33', label=r'$B(a=2R)$')
#plt.plot(x_space, f05r, '-', color='0.66', label=r'$B(a=R/2)$')
#plt.errorbar(x_values, b2r_values, yerr=np.array(c_uncr), linestyle="None", marker='.', color='0', label="Målepunkter")
#plt.errorbar(x_values, b22r_values, yerr=np.array(c_unc2r), linestyle="None", marker='.', color='0.33')
#plt.errorbar(x_values, b205r_values, yerr=np.array(c_unc05r), linestyle="None", marker='.', color='0.66')


#plt.ylabel(r"$B(x)$ (Gauss)", size=18)

plt.show()


