import scipy.integrate
import numpy as np
import matplotlib.pyplot as plt
import time
from math import pi
__author__ = 'alex'


def lldt(m, t, hx, hy, hz, l, g):
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
    return [dmx*1E-5, dmy*1E-5, dmz*1E-5]

"""
def landau_lifshitz(m, t, hx, hy, hz, l, g):
    if callable(hx):
        hx = hx(t)
    if callable(hy):
        hy = hy(t)
    if callable(hz):
        hz = hz(t)
    mz = m[2]
    h = (hx, hy, hz)
    landau1 = tuple_scale(tuple_cross(m, h), -abs(g))
    landau2 = tuple_scale(tuple_cross(m,tuple_cross(m, h)), -l)
    return list(tuple_scale(tuple_add(landau1, landau2), 1E-5))
"""

def tuple_add(a, b):
    a1, a2, a3 = a
    b1, b2, b3 = b
    return a1 + b1, a2 + b2, a3 + b3


def tuple_scale(a, s):
    a1, a2, a3 = a
    return s*a1, s*a2, s*a3


def tuple_cross(a, b):
    a1, a2, a3 = a
    b1, b2, b3 = b
    c1 = a2*b3 - a3*b2
    c2 = a3*b1 - a1*b3
    c3 = a1*b2 - a2*b1
    return c1, c2, c3


def delayed_exp_approach(t, t0, tau):
    if t < t0:
        return 0
    else:
        return 1 - np.exp(-(t - t0) / tau)


def sigmoid_approach(t, t0, tau):
    t = np.array(t)
    return np.divide(1, 1 + np.exp(-(t - t0)/tau))

oe2am = 1000.0 / (4*pi)
emucc2am = 1000.0

hstim = 4 * oe2am
ms = 800.0 * emucc2am
hb = 20 * oe2am
hk0 = 10 * oe2am
g0 = 2 * pi * 28
m0 = [ms, 0, 0]
hx0, hy0, hz0 = hb + hk0, lambda t1: hstim*sigmoid_approach(t1, 1, 0.035), 0
#landau = 8
gilbert = 0.5
l0 = g0 * gilbert / ms
print(l0)
#l0 = 0.002 # alpha
# set up ODE solver
ll = lambda t1, m: lldt(m, t1, hx0, hy0, hz0, l0, g0)
ode15s = scipy.integrate.ode(ll)
ode15s.set_integrator('lsoda')
ode15s.set_initial_value(m0, 0)
npts = 100000
endt = 10
dt = endt/npts
sol = []
t = []
startTime = time.time()
nint = 0
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
dmx = np.diff(mxsol)
dmz = np.diff(mzsol)
phi = np.arctan2(mysol, mxsol)

#plt.plot(t, mxsol, label='mx')
#plt.plot(t, mysol, label='$M_y$')
#plt.plot(t, mzsol, label='mz')
#plt.plot(t, hxvect, label='Hx')
#plt.plot(t, phi, label=r'$\phi$')
#plt.plot(t, hyvect, label=r'$H_y$')
plt.plot(t[:-1], np.log(np.abs(dmy)), label=r'$\frac{d}{dt}M_y$')
plt.legend(loc='best')
plt.xlabel('time (ns)')
plt.grid()
plt.show()