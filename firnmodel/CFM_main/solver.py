import numpy as np
from scipy import interpolate
from scipy.sparse import spdiags
import scipy.sparse.linalg as splin
from constants import *


def solver(a_U, a_D, a_P, b):
	'''
	function for solving matrix problem

	:param a_U:
	:param a_D:
	:param a_P:
	:param b:

	:return phi_t:
	'''

	nz = np.size(b)

	diags = (np.append([a_U, -a_P], [a_D], axis = 0))
	cols = np.array([1, 0, -1])

	big_A = spdiags(diags, cols, nz, nz, format = 'csc')
	big_A = big_A.T

	rhs = -b
	phi_t = splin.spsolve(big_A, rhs)

	return phi_t

def transient_solve_TR(z_edges_vec, z_P_vec, nt, dt, Gamma_P, phi_0, nz_P, nz_fv, phi_s, tot_rho, airdict=None):
	'''
	transient 1-d diffusion finite volume method

	:param z_edges_vec:
	:param z_P_vec:
	:param nt:
	:param dt:
	:param Gamma_P:
	:param phi_0:
	:param nz_P:
	:param nz_fv:
	:param phi_s:

	:return phi_t:
	'''

	phi_t = phi_0

	for i_time in range(nt):
		Z_P = z_P_vec

		dZ = np.concatenate(([1], np.diff(z_edges_vec), [1]))

		dZ_u = np.diff(Z_P)
		dZ_u = np.append(dZ_u[0], dZ_u)
		
		dZ_d = np.diff(Z_P)
		dZ_d = np.append(dZ_d, dZ_d[-1])
		f_u = np.append(0, (1 - (z_P_vec[1:] - z_edges_vec) / dZ_u[1:]))
		f_d = np.append(1 - (z_edges_vec - z_P_vec[0: -1]) / dZ_d[0: -1], 0)

		# Gamma_U = np.append(Gamma_P[0], Gamma_P[0: -1] )
		# Gamma_D = np.append(Gamma_P[1:], Gamma_P[-1])

		# Gamma_u =  1 / ((1 - f_u) / Gamma_P + f_u / Gamma_U)
		# Gamma_d =  1 / ((1 - f_d) / Gamma_P + f_d / Gamma_D)

		##################
		if airdict!=None: # gas diffusion takes a bit more physics
			Gamma_Po    = Gamma_P * airdict['por_op'] #? not here

			Gamma_U 	= np.append(Gamma_Po[0], Gamma_Po[0: -1] )
			Gamma_D 	= np.append(Gamma_Po[1:], Gamma_Po[-1])

			Gamma_u 	=  1 / ((1 - f_u) / Gamma_Po + f_u / Gamma_U)
			Gamma_d 	=  1 / ((1 - f_d) / Gamma_Po + f_d / Gamma_D)

			d_eddy_P    = airdict['d_eddy'] * airdict['por_op']

			d_eddy_U    = np.append(d_eddy_P[0], d_eddy_P[0:-1] )
			d_eddy_D    = np.append(d_eddy_P[1:], d_eddy_P[-1])
		
			d_eddy_u    =  1/ ( (1 - f_u)/d_eddy_P + f_u/d_eddy_U )
			d_eddy_d    =  1/ ( (1 - f_d)/d_eddy_P + f_d/d_eddy_D )
			
			if airdict['gravity']=="off" and airdict['thermal']=="off":
				S_C_0   = 0.0

			elif airdict['gravity']=='on' and airdict['thermal']=='off':
				S_C_0   = (-Gamma_d + Gamma_u) * (airdict['deltaM'] * GRAVITY / (R * airdict['Tz'])) / airdict['dz'] #S_C is independent source term in Patankar

			elif airdict['gravity']=='on' and airdict['thermal']=='on':
				dTdz    = np.gradient(airdict['Tz'])/airdict['dz']
				S_C_0   = (Gamma_d-Gamma_u) * ((-airdict['deltaM'] * GRAVITY / (R * airdict['Tz'])) + (airdict['omega'] * dTdz)) / airdict['dz'] # Check thermal - should it still work in LIZ? if so use d_eddy+diffu
			
			S_C         = S_C_0 * phi_t #this line might be the troublesome one! Should it be phi_0 instead?
			b_0         = S_C * dZ

			rho_interface = np.interp(z_edges_vec,Z_P,airdict['rho'])
			
			w_edges = w(airdict,z_edges_vec)
			w_edges[z_edges_vec>airdict['z_co']] = 0.0
			
			w_u = np.append(w_edges[0], w_edges )
			w_d = np.append(w_edges, w_edges[-1])
			
			D_u = ((Gamma_u+d_eddy_u) / dZ_u) #check signs
			D_d = ((Gamma_d+d_eddy_d) / dZ_d)
		 
			F_u =  w_u * airdict['por_op']
			F_d =  w_d * airdict['por_op']

			# F_u = 0.0 * F_u 
			# F_d = 0.0 * F_d 
			
			P_u = F_u/ D_u
			P_d = F_d/ D_d
			
			a_U = D_u * A( P_u ) + F_upwind(  F_u )
			a_D = D_d * A( P_d ) + F_upwind( -F_d )

			# a_U = D_u # 8/14/17: use this for now - check on Lagrangian need for upwinding.
			# a_D = D_d 
		
			a_P_0 = airdict['por_op'] * dZ/dt
		#######################################

		else: #just for heat, isotope diffusion
			Gamma_U = np.append(Gamma_P[0], Gamma_P[0: -1] )
			Gamma_D = np.append(Gamma_P[1:], Gamma_P[-1])

			Gamma_u =  1 / ((1 - f_u) / Gamma_P + f_u / Gamma_U) # Patankar eq. 4.9
			Gamma_d =  1 / ((1 - f_d) / Gamma_P + f_d / Gamma_D)

			S_C = 0
			S_C = S_C * np.ones(nz_P)

			D_u = (Gamma_u / dZ_u)
			D_d = (Gamma_d / dZ_d)

			b_0 = S_C * dZ

			a_U = D_u 
			a_D = D_d 

			a_P_0 = tot_rho * dZ / dt
			# a_P_0 = tot_rho * dZ / dt
			# a_P_0 = RHO_I * c_firn * dZ / dt

		S_P = 0.0
		a_P =  a_U + a_D + a_P_0 - S_P*dZ

		bc_u_0 = phi_s # need to pay attention for gas
		bc_type = 1
		bc_u   = np.concatenate(([ bc_u_0], [bc_type]))

		bc_d_0 = 0
		bc_type = 2
		bc_d   = np.concatenate(([ bc_d_0 ], [ bc_type ]))
		b = b_0 + a_P_0 * phi_t

		#Upper boundary
		a_P[0] = 1
		a_U[0] = 0
		a_D[0] = 0
		b[0] = bc_u[0]

		#Down boundary
		a_P[-1] = 1
		a_D[-1] = 0
		a_U[-1] = 1
		b[-1] = dZ_u[-1] * bc_d[0]

		phi_t = solver(a_U, a_D, a_P, b)
		a_P = a_U + a_D + a_P_0
		# phi_s = phi_t[0]

	return phi_t

'''
Functions below are for firn air
'''

def w(airdict,z_edges_vec): # Function for downward advection of air and also calculates total air content. 
# def w(z_edges,Accu,rho_interface,por_op,T,p_a,por_tot,por_cl,z_nodes,ad_method,dz): # Function for downward advection of air and also calculates total air content. 

	# global teller_co, por_cl_interface
	
	# por_tot_interface=np.interp(z_edges_vec,z_nodes,airdict['por_tot'])
	# por_cl_interface=np.interp(z_edges_vec,z_nodes,por_cl)
	por_op_interface=np.interp(z_edges_vec,airdict['z'],airdict['por_op'])
	dPdz = np.gradient(airdict['air_pressure'],airdict['dz'])
	dPdz_interface=np.interp(z_edges_vec,airdict['z'],dPdz)
	# teller_co=np.argmax(por_cl_interface)
	# w_ice=Accu*rho_i/rho_interface #Check this - is there a better way?
	# w_ice=1*w_ice
	
	# elif ad_method = 'UW':
	perm = 10.0**(-7.29) * por_op_interface**3.71 # Adolph and Albert, 2014, eq. 5
	# print('perm',perm[100])

	visc = 1.5e-5 #kg m^-1 s^-1, dynamic viscosity

	
	flux = -1.0 * perm / visc * dPdz_interface
	w_ad = flux / airdict['dt']  / por_op_interface
	# print(por_op_interface[np.where(z_edges_vec>43.0)[0][0]])
	# print(w_ad[np.where(z_edges_vec>5.0)[0][0]])
	return w_ad #, bubble_pres

def A(P): # Power-law scheme, Patankar eq. 5.34
	A = np.maximum( (1 - 0.1 * np.abs( P ) )**5, np.zeros(np.size(P) ) )
	return A    

def F_upwind(F): # Upwinding scheme
	F_upwind = np.maximum( F, 0 )
	return F_upwind