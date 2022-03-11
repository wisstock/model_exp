#!/usr/bin/env python3
""" Copyright © 2022 Borys Olifirov
Reproduction of results Markram et al. (1998) with Sympy library.

"""
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# Variative model parameters
# t - s, time
# C_ca_f - nM, free Ca concentration
# C_ca_t - nM, total Ca concentration

# Global model parameters
r_val = 1           # um, radius of compartment
A_val = 20*np.pi    # um^2, surface area of compartment
V_val = 10*np.pi    # um^3, volume of compartment
d_val = 0.04        # um, thinkness of shells
I_val = 130         # pA, maximum inward Ca current
F_val = 96485       # C/mol, Faraday's constant
D_ca_val = 220      # um^2/s, cytoplasm diffusion constant for Ca
v_max_val = 75      # pmol cm^2)/s, maximum pump velocity
K_m_val = 1         # uM, pump Michaelis-Menten constant
K_d_val = 1         # uM, RP dissociation constant
C_ca_rest_val = 50  # nM, resing Ca level

# RP model parameters







# # simple exponential equation 
# sym.init_printing()
# t, l = sym.symbols('t lambda')
# y = sym.Function('y')(t)
# dydt = y.diff(t)
# expr = sym.Eq(dydt, -l*y)
# d_expr = sym.dsolve(expr)

# sym.pprint(expr)
# sym.pprint(d_expr)




# # 2nd order RP kinetic (Markram, 1998)

# t, I, F, V = sp.symbols('t I_ion F V')
# C_Ca = sp.Function('[Ca]')(t)
# dCa_dt = C_Ca.diff(t)
# influx_eq = sp.Eq(dCa_dt, -I/(2*F*V))
# influx_solv = sp.dsolve(influx_eq, C_Ca)
# fin_solv = influx_solv.subs(sp.Symbol('C1'), C_Ca.subs(t, 0))
# i_solv = sp.solve(fin_solv, I)


# sp.pprint(influx_eq, use_unicode=True)
# sp.pprint(influx_solv, use_unicode=True)
# sp.pprint(fin_solv, use_unicode=True)
# sp.pprint(fin_solv.subs({C_Ca.subs(t, 0) : 50, F : 96485, V : 10}), use_unicode=True)

# ca_influx = sp.lambdify(t, fin_solv.subs({C_Ca.subs(t, 0) : 50, F : 96485, V : 10, I : 130}).rhs, "numpy")
# t = np.arange(1, 10, 1)
# fig, ax = plt.subplots(figsize=(4, 4))
# plt.subplot()
# plt.plot(t, ca_influx(t))
# plt.show()