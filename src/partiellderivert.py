import numpy as np
from scipy.misc import derivative
from matplotlib import pyplot as plt

def P_D(func, index, points):
    points_copy = points[:] #Copy points as not to alter the points
    def as_func_of(x):      #Need to input the function as a one variable function
        points_copy[index] = x #The variable is the index'th element
        return func(*points_copy) #Now it is only a function of points[index]. Need star to input a list as arguments
    return derivative(as_func_of, points[index], dx = 1e-6) #See manual to scipy.misc.derivative

def function_f(x1,x2,x3): #Defines a function depending on three variables
    return x3*(1+x1**2/x2**2)**(-3.0/2)

x1_list = np.linspace(-2,2,100) #Creates an array of 100 elements which is evenly distributed between -2 and 2
x1_value = x1_list[74] #one of the values in x1_list
x2_constant = 0.05 #x2 does not change, but the function is dependent on this value
x3_constant = 0.75 #x3 does not change, but the function is dependent on this value

# Take the partial derivative of the function with respect to x1, and input the values in the list [A,B,C]. Since Python is 0-indexed (lists/arrays start at 0),
#  the x1 is the 0th variable in the function. The first delx1_value gives the partial derivative in one point.
# The second, delx1_list, gives the partial derivative as an array of 100 elements where x1 varies from -2 to 2.
delx1_value = P_D(function_f, 0, [x1_value, x2_constant, x3_constant])
delx1_list = P_D(function_f, 0, [x1_list, x2_constant, x3_constant])
delx2 = P_D(function_f, 1, [x1_list, x2_constant, x3_constant])
if delx1_value == delx1_list[74]:
    print('This works') #the list should be the same as the value since we chose number 74 as an example.

#plot the result as a function of x1_list
plt.plot(x1_list, delx1_list)
plt.plot(x1_list, delx2)
plt.show()