from math import sin, cos, pi, sqrt
from scipy.special import itairy, airy as sp_airy
from scipy.integrate import dblquad, quad
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines

plt.ion()

phi_ce = pi/2
dphi = 3.27
#beam width/power, quantum effiency parameter  
chi_e = 0.5
#xi = 10
wavelength = 0.8
outputEnergy = 0.1
time = 40
beamDiameter = 15
watts = outputEnergy * (1/time) * (1/(((beamDiameter/2)**2)*pi)) * (10 ** 15) * 100
xi = 8.55 * sqrt(wavelength**2 * (watts/(10**20)))
#quantum energy parameter
b = chi_e / xi


def airy(z):
    if z > 1000: return 0 
    return sp_airy(z)[0]


def airyp(z):
    if z > 1000: return 0
    return sp_airy(z)[1]


def hdot (phi, phi_ce, dphi):
    if abs(phi) > dphi: return 0
    arg1 = phi + phi_ce
    sin1 = sin(arg1)
    cos1 = cos(arg1)
    arg2 = pi*phi/2/dphi
    sin2 = sin(arg2)
    cos2 = cos(arg2)
    return -sin1*cos2*cos2 - pi/dphi * cos1*cos2*sin2


def zfunc2 (t, tau, hdot_f, chi_e):
    hdot = hdot_f(tau)
    return (abs(hdot) * chi_e)**(-2/3), hdot
def zfunc1 (t, tau, hdot_f, chi_e):
    hdot = hdot_f(tau)
    return (1/(1-t) / abs(hdot) / chi_e)**(2/3), hdot
def zfunc (t, tau, hdot_f, chi_e):
    return (t/(1-t) / abs(hdot_f(tau)) / chi_e)**(2/3)

def ai1 (z):
    return 1/3 - itairy (z)[0]

# Use transform from Numerical Recipes eq 4.4.3 to handle the integrable
# singularities at t=0.
def integrand_1 (u, tau, hdot_f, chi_e):
    u2 = u*u
    t = u*u2
    zp,_ = zfunc1 (t, tau, hdot_f, chi_e)
    z = zp*u2
    return u2*ai1(z) + 2*airyp(z)/zp

def int_1 (chi_e = chi_e, umax=1):
    hdot_f = lambda tau, phi_ce=phi_ce, dphi=dphi: hdot(tau, phi_ce, dphi)
    ans, err = dblquad (lambda u, tau, hdot_f=hdot_f, chi_e=chi_e: integrand_1(u, tau, hdot_f, chi_e),
                        -dphi, dphi, # tau
                        0, umax) # u
    return (ans*3, err*3)


def integrand_1b (t, tau, hdot_f, chi_e):
    z = zfunc (t, tau, hdot_f, chi_e)
    return ai1(z) + 2*airyp(z)/z


def xx (t, tau, hdot_f, chi_e):
    z = zfunc (t, tau, hdot_f, chi_e)
    return ai1(z)


def int_1b (chi_e = chi_e, tmax=1):
    hdot_f = lambda tau, phi_ce=phi_ce, dphi=dphi: hdot(tau, phi_ce, dphi)
    ans, err = dblquad (lambda t, tau, hdot_f=hdot_f, chi_e=chi_e: integrand_1b(t, tau, hdot_f, chi_e),
                        -dphi, dphi, # tau
                        0.5, tmax) # t
    return (ans, err)


def integrand_2 (t, tau, hdot_f, chi_e):
    zp,_ = zfunc1 (t, tau, hdot_f, chi_e)
    t23 = t**(2/3)
    z = zp * t23
    return t23*t23/(1-t)*airyp(z)/zp


def int_2 (chi_e = chi_e):
    hdot_f = lambda tau, phi_ce=phi_ce, dphi=dphi: hdot(tau, phi_ce, dphi)
    ans, err = dblquad (lambda t, tau, hdot_f=hdot_f, chi_e=chi_e: integrand_2(t, tau, hdot_f, chi_e),
                        -dphi, dphi, # tau
                        0, 1) # t
    return (ans, err)


def integrand_3 (t, tau, hdot_f, chi_e):
    zp, hdot = zfunc1 (t, tau, hdot_f, chi_e)
    sgn = 1 if hdot >= 0 else -1
    t23 = t**(2/3)
    z = zp * t23
    return t23 / sqrt(zp) * sgn * airy(z)


def int_3 (chi_e = chi_e):
    hdot_f = lambda tau, phi_ce=phi_ce, dphi=dphi: hdot(tau, phi_ce, dphi)
    ans, err = dblquad (lambda t, tau, hdot_f=hdot_f, chi_e=chi_e: integrand_3(t, tau, hdot_f, chi_e),
                        -dphi, dphi, # tau
                        0, 1) # t
    return (ans, err)


def P (xi_zeta, chi_e=chi_e):
    return int_1(chi_e)[0] + int_2(chi_e)[0] + xi_zeta*int_3(chi_e)[0]


def integrand_4 (t, tau, hdot_f, chi_e):
    zp, hdot = zfunc1 (t, tau, hdot_f, chi_e)
    sgn = 1 if hdot >= 0 else -1
    t23 = t**(2/3)
    z = zp * t23
    return t23 / sqrt(zp) * sgn * airy(z) / (1-t)


def int_4 (chi_e = chi_e):
    hdot_f = lambda tau, phi_ce=phi_ce, dphi=dphi: hdot(tau, phi_ce, dphi)
    ans, err = dblquad (lambda t, tau, hdot_f=hdot_f, chi_e=chi_e: integrand_4(t, tau, hdot_f, chi_e),
                        -dphi, dphi, # tau
                        0, 1) # t
    return (ans, err)


#Important Function
def xi_zeta2 (xi_zeta, chi_e=chi_e):
    return (xi_zeta * int_1(chi_e)[0] + int_4(chi_e)[0]) / P(xi_zeta, chi_e)

def plot (f, lo, hi, nx=100):
    x = np.linspace(lo, hi, nx)
    y = list(map(f,x))
    plt.plot (x, y)
    return


def sgn(x):
    if x > 0: return 1
    return -1
def hdotter (phi_ce=phi_ce, dphi=dphi):
    return quad (lambda tau: sgn(hdot(tau, phi_ce, dphi)),
                 -dphi, dphi)

#first number is initial polarization
#second is Chi_e
print(xi_zeta2(0, 0.03))