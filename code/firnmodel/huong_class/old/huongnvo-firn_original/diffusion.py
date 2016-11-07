from solver import transient_solve_TR
from constants import *
import numpy as np

class Diffusion:
    def __init__(self, z, stp, gridLen, init_Tz):
        '''
        Initializes a diffusion class -- can run heat diffusion or isotope diffusion

        :param stp:
        :param gridLen:
        :param init_Tz: initial temperature vector along the depth of the firn column
        :param init_T10m:
        '''

        self.Tz = init_Tz
        self.del_s = -30.0 * np.ones(stp)
        self.del_z = del_s[0] * np.ones(gridLen)

        fT10m      = interpolate.interp1d(z, self.Tz) #temp at 10m depth
        self.T10m  = fT10m(10)

    def heatDiff(self, z, dz, Ts, rho):
        '''
        Heat diffusion function

        :param z:
        :param dz:
        :param Ts:
        :param rho:

        :returns self.Tz:
        :returns self.T10m:
        '''

    	nz_P = len(z)
        nz_fv = nz_P - 2
        nt = 1

        z_edges_vec = z[1:-2] + dz[2:-1] / 2
        z_edges_vec = np.concatenate(([z[0]], z_edges_vec, [z[-1]]))
        z_P_vec = z
        phi_s = Ts
        phi_0 = self.Tz

        K_ice = 9.828 * np.exp(-0.0057 * phi_0)
        K_firn = K_ice * (rho / 1000) ** (2 - 0.5 * (rho / 1000))
        c_firn = 152.5 + 7.122 * phi_0
        Gamma_P = K_firn / (c_firn * rho)

        self.Tz = transient_solve_TR(z_edges_vec, z_P_vec, nt, dt, Gamma_P, phi_0, nz_P, nz_fv, phi_s)
        self.Tz = np.concatenate(([Ts], self.Tz[:-1]))

        fT10m = interpolate.interp1d(z, self.Tz) 								# temp at 10m depth
        self.T10m = fT10m(10)

        return self.Tz, self.T10m

    def isoDiff(self, iter, z, dz, rho, iso):
        '''
        Isotope diffusion function

        :param iter:
        :param z:
        :param dz:
        :param rho:
        :param iso:

        :returns self.phi_t:
        '''
    	nz_P = len(z)       												        # number of nodes in z
        nz_fv = nz_P - 2      												        # number of finite volumes in z
        nt = 1             													        # number of time steps

        z_edges_vec = z[1:-2] + dz[2:-1] / 2        						        # uniform edge spacing of volume edges
        z_edges_vec = np.concatenate(([z[0]], z_edges_vec, [z[-1]]))
        z_P_vec = z

        # Node positions
        phi_s = self.del_s[iter]    												# isotope value at surface
        phi_0 = self.del_z       												    # initial isotope profile

        # Define diffusivity for each isotopic species
        # Establish values needed for diffusivity calculation
        m = 0.018           												        # kg/mol; molar mass of water
        pz = 3.454 * np.power(10, 12) * np.exp(-6133 / self.Tz)    				    # Pa; saturation vapor pressure over ice

        alpha_18_z = 0.9722 * np.exp(11.839 / self.Tz)           				    # fractionation factor for 18_O
        # alpha_18_z = np.exp(11.839/Tz-28.224*np.power(10,-3)) 			        # alternate formulation from Eric's python code
        alpha_D_z = 0.9098 * np.exp(16288 / np.power(self.Tz, 2)) 				    # fractionation factor for D
        Po = 1.0              												        # reference pressure in atm
        P = 1.0

        # Set diffusivity in air (units of m^2/s)
        Da = 2.1 * np.power(10.0, -5.) * np.power(self.Tz / 273.15, 1.94) * (Po / P)
        Da_18 = Da / 1.0251     											        # account for fractionation factor for 18_O
        Da_D = Da / 1.0285      											        # account for fractionation factor for D

        # Calculate tortuosity
        invtau = np.zeros(gridLen)
        b = 0.25            												        # Tortuosity parameter
        invtau[rho < RHO_I / np.sqrt(b)] = 1.0 - (b * (rho[rho < RHO_I / np.sqrt(b)] / RHO_I)) ** 2
        invtau[rho >= RHO_I / np.sqrt(b)] = 0.0

        # Set diffusivity for each isotope
        if iso == '18':
            D = m * pz * invtau * Da_18 * (1 / rho - 1 / RHO_I) / (R * self.Tz * alpha_18_z)
        elif iso == 'D':
            D = m * pz * invtau * Da_D * (1 / rho - 1 / RHO_I) / (R * self.Tz * alpha_D_z)

        # Solve for vertical isotope profile at this time step i
        self.del_z = transient_solve_TR(z_edges_vec, z_P_vec, nt, dt, D, phi_0, nz_P, nz_fv, phi_s)
        self.del_z = np.concatenate(([del_s[ii]], del_z[:-1]))

        return self.del_z