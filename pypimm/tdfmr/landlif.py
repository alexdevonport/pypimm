import scipy.integrate
import numpy as np
import matplotlib.pyplot as plt
import math
import time
from math import pi
__author__ = 'alex'


"""
def landau_lifshitz(m, t, hx, hy, hz, l, g):
    mx, my, mz = m
    if callable(hx):
        hx = hx(t)
    if callable(hy):
        hy = hy(t)
    if callable(hz):
        hz = hz(t)
    dmx = g*(hy*mz - hz*my) - l*(my*(hy*mx - hx*my) - mz*(hx*mz - hz*mx))
    dmy = g*(hz*mx - hx*mz) - l*(mz*(hz*my - hy*mz) - mx*(hy*mx - hx*my))
    dmz = g*(hx*my - hy*mx) - l*(mx*(hx*mz - hz*mx) - my*(hz*my - hy*mz))
    return [dmx, dmy, dmz]
"""

def landau_lifshitz(m, t, hx, hy, hz, l, g):
    if callable(hx):
        hx = hx(t)
    if callable(hy):
        hy = hy(t)
    if callable(hz):
        hz = hz(t)
    h = (hx, hy, hz)
    landau1 = tuple_scale(tuple_cross(m, h), -abs(g))
    landau2 = tuple_scale(tuple_cross(m,tuple_cross(m, h)), -l)
    return list(tuple_add(landau1, landau2))

def tuple_add(a, b):
    a1, a2, a3 = a
    b1, b2, b3 = b
    return (a1 + b1, a2 + b2, a3 + b3)

def tuple_scale(a, s):
    a1, a2, a3 = a
    return (s*a1, s*a2, s*a3)

def tuple_cross(a, b):
    a1, a2, a3 = a
    b1, b2, b3 = b
    c1 = a2*b3 - a3*b2
    c2 = a3*b1 - a1*b3
    c3 = a1*b2 - a2*b1
    return (c1, c2, c3)

def delayed_exp_approach(t, t0, tau):
    if t < t0:
        return 0
    else:
        return 1 - np.exp(-(t - t0) / tau)

def slow_approach(t, t0, tau):
    t = np.array(t)
    return np.divide(1, 1 + np.exp(-(t - t0)/tau))

hstim = 1.75
ms = 800.0
hb = 70

g0 = 2 * pi * 28E-3
m0 = [ms, 0, 0]
hx0, hy0, hz0 = hb, lambda t1: hstim*slow_approach(t1, 1, 0.035), 0
l0 = 30E-3 * g0 / (4 * pi * ms * ms)
#l0 = 0
#damplambda = 130E-3
#l0 = damplambda / (4*pi*1E-7 * ms*ms)
print(l0)

# set up ODE solver
ll = lambda t1, m: landau_lifshitz(m, t1, hx0, hy0, hz0, l0, g0)
ode15s = scipy.integrate.ode(ll)
ode15s.set_integrator('vode', method='bdf', order=3, nsteps=8000, max_step=1E-3)
ode15s.set_initial_value(m0, 0)
npts = 512
endt = 10
dt = endt/npts
sol = []
t = []
startTime = time.time()
while ode15s.successful() and ode15s.t < endt:
    a = ode15s.integrate(ode15s.t+dt)
    t.append(ode15s.t)
    sol.append(a)
sol = np.asarray(sol)
elapsedTime = (time.time() - startTime) * 1E3
print('Solution of LL took {:3.2f} ms.'.format(elapsedTime))
mxsol = sol[:,0]
mysol = sol[:,1]
mzsol = sol[:,2]

hyvect =  hy0(t)


dmy = np.diff(mysol)
phi = np.arctan2(mysol, mxsol)

#plt.plot(t, mxsol, label='mx')
#plt.plot(t, mysol, label='$M_y$')
#plt.plot(t, mzsol, label='mz')
#plt.plot(t, hxvect, label='Hx')
#plt.plot(t, phi, label=r'$\phi$')
plt.plot(t, hyvect, label=r'$H_y$')
plt.plot(t[:-1], dmy, label=r'$\frac{d}{dt}M_y$')
plt.legend(loc='best')
plt.xlabel('time (ns)')
plt.grid()
plt.show()