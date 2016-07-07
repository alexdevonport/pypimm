__author__ = 'alex'
import numpy as np
from math import pi
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')
print(plt.style.available)

def dsinplus_sp(x, p0, p1, p2, p3, p4, p5):
    """
    Exactly the same as dsin, but with separate args because that's what Scipy's curve_fit takes.
    :param x:
    :param p:
    :return:
    """
    u = np.piecewise(x, [x < p5*1E-9, x >= p5*1E-9], [0, 1])
    x = np.array(x)
    t1 = p0
    t2 = np.cos((x-p5*1E-9)*2*pi*p1*1E9)
    t3 = np.exp(-(x-p5*1E-9)/p2*1E9)
    t4 = p3
    t5 = p3 * np.exp(-(x-p5*1E-9)/p4*1E9)
    y = t1*t2*t3 + t4*t5
    y = np.multiply(u, y)
    return y

x = np.linspace(0, 3E-9, 512)

test1 = dsinplus_sp(x , 1,  1,  1,  0.3,  1,  1)
plt.plot(x, test1)
plt.xlabel('time (s)')
plt.ylabel('signal (arb. units)')
plt.title('Simulated Data')
plt.show()