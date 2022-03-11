#!/usr/bin/env python3
""" Copyright © 2022 Borys Olifirov
Replication of results Markram et al. (1998) with Sympy library.

"""
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# Variative model parameters
# t - s, time
# C_ca_f - nM, free Ca concentration
# C_ca_t - nM, total Ca concentration
# C_rp_b - uM, Ca-bound form concentration of individual RP

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
C_rp_t_val = 120    # uM, total concentration on individual RP
k_on_2_val = 1e8    # 1/(mol*s), association rate constant of intermediate RP
k_off_2_val = 1e2   # 1/s, dissociation rate constant of intermediate RP







# # simple exponential equation 
# sym.init_printing()
# t, l = sym.symbols('t lambda')
# y = sym.Function('y')(t)
# dydt = y.diff(t)
# expr = sym.Eq(dydt, -l*y)
# d_expr = sym.dsolve(expr)

# sym.pprint(expr)
# sym.pprint(d_expr)


# simple model
#   k_i
# A ---> B
#   <---
#   k_o
t, k_i, k_o, A_0, B_0 = sp.symbols('t k_i k_o [A]_0 [B]_0')
At = sp.Function('[A]')(t)
dA_dt = At.diff(t)
A_eq = sp.Eq(dA_dt, -k_i*At + k_o * (A_0 - At))
A_dsolv = sp.dsolve(A_eq)

sp.pprint(A_eq)
sp.pprint(A_dsolv)
sp.pprint(A_dsolv.subs(t, 0))


# # 2nd order RP kinetic (Markram, 1998)

# t, k_on, k_off, C_ca, C_p, C_t = sp.symbols('t k_i k_o C_ion C_p C_t')
# Cb = sp.Function('C_m')(t)
# dCb_dt = Cb.diff(t)
# rp_eq = sp.Eq(- dCb_dt, -k_on * C_ca * (C_t - Cb) + k_off * Cb)
# rp_dsolv = sp.dsolve(rp_eq, Cb)

# t_0 = rp_dsolv.subs(t, 0)
# # C1_eq = sp.Eq(sp.solve(t_0, sp.Symbol('C1'))[0], C1)
# fin_solv = rp_dsolv.subs(sp.Symbol('C1'), C_t)
# # i_solv = sp.solve(fin_solv, I)


# sp.pprint(rp_eq, use_unicode=True)
# sp.pprint(rp_dsolv, use_unicode=True)
# # sp.pprint(fin_solv, use_unicode=True)
# # # sp.pprint(fin_solv.subs({C_Ca.subs(t, 0) : 50, F : 96485, V : 10}), use_unicode=True)

# rp_num = sp.lambdify(t, fin_solv.subs({C_t : C_rp_t_val, k_on : k_on_2_val, k_off : k_off_2_val}).rhs, "numpy")
# t = np.arange(1, 10, 1)
# fig, ax = plt.subplots(figsize=(4, 4))
# plt.subplot()
# plt.plot(t, rp_num(t))
# plt.show()