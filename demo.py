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

HPCA parameters (Furuta 1999, O'Callaghan 2003):
C_rp - total reaction partner concentration (35 uM)
C_rp_f - free reaction partner concentration (uM)
C_rp_b - Ca-bound form of reaction partner concentration (uM)
D_rp - cytoplasm diffusion constant for reaction partner (30 um^2/s)
n - reaction partner Hill coefficient (1.4)
K_d - reaction partner dissociation constant (0.32 uM)
t_on - insertion in membrane time (5e10^-2 s)
t_off - leaving membrane time (3 s)
k_on - 1/t_on, reaction partner association rate constant (1/s)
k_off - 1/t_off, reaction partner dissociation rate constant (1/s)


Useful links:
https://www.cfm.brown.edu/people/dobrush/am33/SymPy/index.html
https://www.sympy.org/scipy-2017-codegen-tutorial/notebooks/20-ordinary-differential-equations.html
https://www.sympy.org/scipy-2017-codegen-tutorial/notebooks/25-chemical-kinetics-intro.html
https://dabane-ghassan.github.io/ModNeuro/
https://dabane-ghassan.github.io/ModNeuro/images/1_Enzyme_Kinetics.html

"""

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


def hill_fam(K_val, n_list, c=np.arange(0, 1, 0.01)):
	""" Hill equation solutions family.

	"""
	n, K_A, C_ca = sp.symbols('n  K_a C_ion')
	Y = sp.Function('Y')(C_ca)
	hill_eq = sp.Eq(Y, C_ca**n/(K_A**n + C_ca**n))
	hill_sol_assemblage = [sp.lambdify(C_ca, hill_eq.subs({n : n_val, K_A : K_val}).rhs, "numpy") for n_val in n_list]

	sp.pprint(hill_eq, use_unicode=True)

	fig, ax = plt.subplots(figsize=(4, 4))
	plt.subplot()
	for i in range(0, len(hill_sol_assemblage)):
		plt.plot(c, hill_sol_assemblage[i](c), label=f'n={n_list[i]}')
	plt.xlabel('C')
	plt.ylabel('Y')
	plt.legend()
	plt.show()
	
# hill_fam(K_val=0.5, n_list=[-2, -1, 1, 2])  # np.arange(1, 6, 1)


def rp_fam(K_val, K_list, a_list, n_val=1.5, c=np.arange(0, 1, 0.01)):
	""" Shurik's HPCA Kd estimation model with diffusion distance variation vs. Hill model with Kd variation.
   Use simple Y estimation, isn't correct.

	"""
	n, K_d, a, D, t_off, C_ca, C_tot = sp.symbols('n K a D tau_aut [Ca] C_tot')
	C_b = sp.Function('C_m')(C_ca)

    # compartment model
	rp_eq = sp.Eq(C_b, C_tot*(C_ca**2 / (C_ca**2 + ((a**2/D)/t_off) **(1/n) * K_d)) )
	sp.pprint(rp_eq, use_unicode=True)
	rp_sol_assemblage = [sp.lambdify(C_ca, rp_eq.subs({n :n_val, K_d : K_val, a : a_val, D : 30, C_tot : 35, t_off : 3}).rhs, "numpy")
                     for a_val in a_list]

    # Hill model
	hill_eq = sp.Eq(C_b, C_tot * (C_ca**n/(K_d + C_ca**n)))
	sp.pprint(hill_eq, use_unicode=True)
	# hill_sol = sp.lambdify(C_ca, hill_eq.subs({n : n_val, K_d : K_val, C_tot : 35}).rhs, "numpy")
	hill_sol_assemblage = [sp.lambdify(C_ca, hill_eq.subs({n : n_val, K_d : K_val, C_tot : 35}).rhs, "numpy") for K_val in K_list]

	fig, ax = plt.subplots(figsize=(4, 4))
	plt.subplot()
	for i in range(0, len(rp_sol_assemblage)):
		plt.plot(c, rp_sol_assemblage[i](c), label=f'a={a_list[i]}um')
	for h in range(0, len(hill_sol_assemblage)):
		plt.plot(c, hill_sol_assemblage[h](c), linestyle='--', label=f'Kd={round(K_list[h], 1)}mM')
	plt.xlabel('[Ca] (uM)')
	plt.ylabel('Cm (uM)')
	plt.legend()
	plt.show()

# rp_fam(K_val=0.3, K_list=np.arange(0.1, 0.8, 0.2), a_list=np.arange(1, 9, 2), n_val=1.5, c=np.arange(0, 1, 0.01))


def sh_mod(K_val, a_val, K_list, a_list, n_val=1.5, C_ca_range=np.arange(0, 1, 0.01)):
	""" Pure Shutik's kinetic model.
	Series with various a (diffusion distance) and cytoplasmic Kd. 

	"""
	K_d, C_ca, n, t_out, a, D = sp.symbols('K [Ca] n tau_-1 a D')
	Y = sp.Function('Y')(C_ca)
	y_eq = sp.Eq(Y, (t_out/(a**2/D))*(C_ca/K_d)**n/(1 + (C_ca/K_d)**n + (t_out/(a**2/D))*(C_ca/K_d)**n))  # full model
	sp.pprint(y_eq, use_unicode=True)

	a_var_assemblage = [sp.lambdify(C_ca, y_eq.subs({K_d : K_val, n : n_val, t_out : 1, a : one_a, D : 30}).rhs, "numpy") for one_a in a_list] 
	K_var_assemblage = [sp.lambdify(C_ca, y_eq.subs({K_d : one_K, n : n_val, t_out : 1, a : a_val, D : 30}).rhs, "numpy") for one_K in K_list]

	fig, ax = plt.subplots(figsize=(4, 4))
	plt.subplot()
	for h in range(0, len(K_var_assemblage)):
		plt.plot(C_ca_range, K_var_assemblage[h](C_ca_range), linestyle=':', linewidth=2.5, label=f'Kd={round(K_list[h], 1)}mM')
	for i in range(0, len(a_var_assemblage)):
		plt.plot(C_ca_range, a_var_assemblage[i](C_ca_range), linestyle='--', label=f'a={a_list[i]}um')
	plt.xlabel('[Ca] (uM)')
	plt.ylabel('Y')
	plt.legend()
	plt.show()

# sh_mod(K_val=0.3, a_val=1, K_list=[0.1, 0.3, 0.6], a_list=[0.2, 0.5, 1], n_val=1.3, C_ca_range=np.arange(0, 0.5, 0.001))


def sh_k_d(K_d_lim, a_lim, n_val=1.5, t_out_val=1, D_val=30):
	""" Membrane Kd estimation

	"""
	K_d, n, t_out, a, D = sp.symbols('K n tau_-1 a D')
	K_m = sp.Function('K_m')(K_d, a)
	K_m_eq = sp.Eq(K_m, (t_out/(a**2/D))**(1/n) * K_d)
	K_m_lambda = sp.lambdify([K_d, a], K_m_eq.subs({t_out : t_out_val, n : n_val, D : D_val}).rhs, 'numpy')
	sp.pprint(K_m_eq)

	K_d_X, a_Y = np.meshgrid(np.linspace(K_d_lim[0], K_d_lim[1]), np.linspace(a_lim[0], a_lim[1]))
	K_m_Z = K_m_lambda(K_d_X, a_Y)
	print(K_m_Z.shape)

	levels = np.linspace(np.min(K_m_Z), np.max(K_m_Z))
	fig, ax = plt.subplots()
	ax.contourf(K_d_X, a_Y, K_m_Z, levels=levels)
	plt.show()

sh_k_d(K_d_lim=[0.3, 1.5], a_lim=[0.01, 0.1]) 