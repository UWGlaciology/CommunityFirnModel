import numpy as np
from scipy import interpolate
from scipy.sparse import spdiags
import scipy.sparse.linalg as splin

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

def transient_solve_TR(z_edges_vec, z_P_vec, nt, dt, Gamma_P, phi_0, nz_P, nz_fv, phi_s):
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

    for i_time in xrange(nt):
        Z_P = z_P_vec

        dZ = np.concatenate(([1], np.diff(z_edges_vec), [1]))

        dZ_u = np.diff(Z_P)
        dZ_u = np.append(dZ_u[0], dZ_u)

        dZ_d = np.diff(Z_P)
        dZ_d = np.append(dZ_d, dZ_d[-1])

        f_u = np.append(0, (1 - (z_P_vec[1:] - z_edges_vec) / dZ_u[1:]))
        f_d = np.append(1 - (z_edges_vec - z_P_vec[0: -1]) / dZ_d[0: -1], 0)


        Gamma_U = np.append(Gamma_P[0], Gamma_P[0: -1] )
        Gamma_D = np.append(Gamma_P[1:], Gamma_P[-1])


        Gamma_u =  1 / ((1 - f_u) / Gamma_P + f_u / Gamma_U)
        Gamma_d =  1 / ((1 - f_d) / Gamma_P + f_d / Gamma_D)

        S_C = 0
        S_C = S_C * np.ones(nz_P)

        b_0 = S_C * dZ

        D_u = (Gamma_u / dZ_u)
        D_d = (Gamma_d / dZ_d)
        a_U = D_u
        a_D = D_d

        a_P_0 = dZ / dt

        a_P =  a_U + a_D + a_P_0

        bc_u_0 = phi_s
        bc_type = 1.
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

    return phi_t