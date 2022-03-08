#!/usr/bin/env python3
""" Copyright © 2022 Borys Olifirov
Experiments with calcium sensors kinetic models.

Conda env 'model'

Global parameters:
t - time (s)
r - radius of compartment (1 um)
A - surface area of compartment (20* pi um^2)
V - volume of compartment (10*pi um^3)
d - thinkness of shells (0.04 um)
I - Ca current (130 pA)
F - Faraday's constant (96485 e/mol)
C_f - free Ca concentration (nM)
C_r - resting Ca level (50 nM)
C_t - total Ca concentration (nM)
D_ca - cytoplasm diffusion constant for Ca (220 um^2/s)
v_max - maximum pump velocity (75 (pmol cm^2)/s)
K_m - pump Michaelis-Menten constant (1 uM)

Sensor parameters:
C_rp - total reaction partner concentration (uM)
C_rp_f - free reaction partner concentration (uM)
C_rp_b - Ca-bound form of reaction partner concentration (uM)
K_d - reaction partner dissociation constant (uM)
k_on - reaction partner association rate constant (1/(M*s))
k_off - reaction partner dissociation rate constant (1/s)


Useful links:
https://www.cfm.brown.edu/people/dobrush/am33/SymPy/index.html
https://www.sympy.org/scipy-2017-codegen-tutorial/notebooks/20-ordinary-differential-equations.html

"""

import numpy as np
import sympy as sym
from sympy.plotting import plot
# import matplotlib
# import matplotlib.pyplot as plt


# # simple exponential equation 
# sym.init_printing()
# t, l = sym.symbols('t lambda')
# y = sym.Function('y')(t)
# dydt = y.diff(t)
# expr = sym.Eq(dydt, -l*y)
# d_expr = sym.dsolve(expr)

# sym.pprint(expr)
# sym.pprint(d_expr)


# Ca influx model
t, I, V, F = sym.symbols('t I V F')
# F = 96485
# V = 1
C_ca = sym.Function('C')(t)
dc_dt = C_ca.diff(t)
ca_influx = sym.Eq(dc_dt, -I/(2*F*V))
ca_influx_i = sym.dsolve(ca_influx, C_ca, ics={C_ca.subs(t, 0): 0})
ca_current = sym.solve(ca_influx_i, I)
# influx_fd = sym.as_finite_diff(ca_influx)

sym.pprint(ca_influx)
sym.pprint(ca_influx_i)
sym.pprint(ca_current)
