import numpy as np
import time
import scipy.special as special
from scipy.integrate import quad, dblquad, fixed_quad

#constants
epsilon0  = 1./(4. * np.pi)
e = 4.80320680e-10
m = 9.1093826e-28
c = 2.99792458e10

epsilon  = -1.    #sign of electron charge

#parameters
B       = 1.         #background B strength
n_e     = 1.         #electron number density cm^-3
theta_e = 0.5         #dimensionless electron temp
theta   = np.pi/3.    #observer angle

#derived quantities
omega_p = np.sqrt(n_e * e**2. / (m * epsilon0))     # plasma frequency    (=1 in these units)
omega_c = e * B / (m * c)                           # cyclotron frequency (=1 in these units)


# In[53]:

def K_12_prefactor(omega):
    prefactor = - 1. * omega_p**2. / (omega * omega) * 1./(4. * theta_e**2. * special.kn(2, 1./theta_e))
    return prefactor

def K_12_integrand(gamma, cos_xi, tau_prime, omega):    
    prefactor  = 1j
    beta       = np.sqrt(1. - 1./gamma**2.)
    sin_xi     = np.sqrt(1. - cos_xi**2.)
    p_perp_bar = gamma * beta * sin_xi
    p_z_bar    = gamma * beta * cos_xi
    beta_perp = p_perp_bar / gamma
    beta_z    = p_z_bar    / gamma
    b         = omega/(epsilon * omega_c) * np.sin(theta) * gamma * beta_perp
    
    term1 = p_perp_bar**2. / gamma**2. * np.exp(-gamma/theta_e)
    term2 = np.exp(1j * (1. - beta_z * np.cos(theta)) * tau_prime)
    term3 = np.sin((epsilon * omega_c / omega) * tau_prime / gamma) * special.j0(2. * b * np.sin((epsilon * omega_c / omega) * tau_prime / (2. * gamma)))
    ans   = prefactor * term1 * term2 * term3
    return ans * gamma**2. * beta 

def K_32_prefactor(omega):
    prefactor = omega_p**2. / (omega * omega) * 1./(2. * theta_e**2. * special.kn(2, 1./theta_e))
    return prefactor

def K_32_integrand(gamma, cos_xi, tau_prime, omega):
    prefactor = 1.
    beta       = np.sqrt(1. - 1./gamma**2.)
    sin_xi     = np.sqrt(1. - cos_xi**2.)
    p_perp_bar = gamma * beta * sin_xi
    p_z_bar    = gamma * beta * cos_xi
    beta_perp = p_perp_bar / gamma
    beta_z    = p_z_bar    / gamma
    b         = omega/(epsilon * omega_c) * np.sin(theta) * gamma * beta_perp
    
    term1 = p_perp_bar * p_z_bar / gamma**2. * np.exp(-gamma/theta_e)
    term2 = np.exp(1j * (1. - beta_z * np.cos(theta)) * tau_prime)
    term3 = np.sin((epsilon * omega_c / omega) * tau_prime / (2. * gamma)) * special.jn(1, 2. * b * np.sin((epsilon * omega_c / omega) * tau_prime / (2. * gamma)))
    ans   = prefactor * term1 * term2 * term3
    return ans * gamma**2. * beta 

# In[54]:

def K_12_p_integrated(tau_prime, omega):
    real_part = dblquad(lambda gamma, cos_xi: K_12_integrand(gamma, cos_xi, tau_prime, omega).real, 
                       -1, 1, lambda x: 1., lambda x: np.inf)
    ans = real_part[0]
    return ans


def K_32_p_integrated(tau_prime, omega):
    real_part = dblquad(lambda gamma, cos_xi: K_32_integrand(gamma, cos_xi, tau_prime, omega).real, 
                       -1, 1, lambda x: 1., lambda x: np.inf)
    ans = real_part[0]
    return ans

#---------------------------------------------------------------------------------------------------#

def n_loc(n, tau_prime, gamma):
    if(gamma == 1.):
        return np.sign(1. - n*np.pi/tau_prime)
    
    beta   = np.sqrt(1. - 1./gamma**2.)
    ans = (1. - n*np.pi/tau_prime)/(beta * np.cos(theta))    
    
    if(ans > 1.):
        return 1.
    if(ans < -1.):
        return -1.
    
    return ans

def K_12_xi_integrator(gamma, tau_prime, omega):
    if(tau_prime == 0.):
        return 0.
    
    beta   = np.sqrt(1. - 1./gamma**2.)
    n_min  = - int(tau_prime/np.pi * (-1. + beta * np.cos(theta)))
    n_max  =   int(tau_prime/np.pi * (1. + beta * np.cos(theta)))
    
    ans = 0.
    
    ans    = quad(lambda cos_xi: K_12_integrand(gamma, cos_xi, tau_prime, omega).real, 
                      n_loc(n_max, tau_prime, gamma), -1.)[0]
    
    ans   += quad(lambda cos_xi: K_12_integrand(gamma, cos_xi, tau_prime, omega).real, 
                      1., n_loc(n_min, tau_prime, gamma))[0]
    
    for i in range(n_min, n_max):
        ans += quad(lambda cos_xi: K_12_integrand(gamma, cos_xi, tau_prime, omega).real, 
                    n_loc(i, tau_prime, gamma), n_loc(i+1, tau_prime, gamma))[0]
    return -ans

def K_32_xi_integrator(gamma, tau_prime, omega):
    if(tau_prime == 0.):
        return 0.

    beta   = np.sqrt(1. - 1./gamma**2.)
    n_min  = - int(tau_prime/np.pi * (-1. + beta * np.cos(theta)))
    n_max  =   int(tau_prime/np.pi * (1. + beta * np.cos(theta)))

    ans = 0.

    ans    = quad(lambda cos_xi: K_32_integrand(gamma, cos_xi, tau_prime, omega).real,
                      n_loc(n_max, tau_prime, gamma), -1.)[0]

    ans   += quad(lambda cos_xi: K_32_integrand(gamma, cos_xi, tau_prime, omega).real,
                      1., n_loc(n_min, tau_prime, gamma))[0]

    for i in range(n_min, n_max):
        ans += quad(lambda cos_xi: K_32_integrand(gamma, cos_xi, tau_prime, omega).real,
                    n_loc(i, tau_prime, gamma), n_loc(i+1, tau_prime, gamma))[0]
    return -ans


def K_12_p_integrator(tau_prime, omega):
    real_part = quad(lambda gamma: K_12_xi_integrator(gamma, tau_prime, omega), 1., np.inf)[0]
    return real_part

def K_32_p_integrator(tau_prime, omega):
    real_part = quad(lambda gamma: K_32_xi_integrator(gamma, tau_prime, omega), 1., np.inf)[0]
    return real_part
#-----------------------------------------------------------------------------#
time_before = time.time()

#omega_mult = np.logspace(0., 3., 10) * omega_c
#omega_mult = [200. * omega_c]
omega = 1. * omega_c

K_12_ans = 0.
K_32_ans = 0.
i        = 0.
j        = 0.
step     = 30.
num_pts  = 45
tolerance = 0.05

while(i == 0. or np.abs(step_ans_12 / K_12_ans) > tolerance):
	time_before_12 = time.time()
	step_ans_12    = fixed_quad(lambda tau_prime: np.vectorize(K_12_p_integrator)(tau_prime, omega), 
		  						   i*step, (i+1)*step, n = num_pts)[0]
	time_after_12  = time.time()
	K_12_ans      += step_ans_12
	print 'K_12 i  =', i*step, (i+1)*step, np.abs(step_ans_12/K_12_ans), time_after_12 - time_before_12
	i             += 1.

while(j == 0. or np.abs(step_ans_32 / K_32_ans) > tolerance):
	time_before_32 = time.time()
        step_ans_32    = fixed_quad(lambda tau_prime: np.vectorize(K_32_p_integrator)(tau_prime, omega), 
								   j*step, (j+1)*step, n = num_pts)[0]
        time_after_32  = time.time()
	K_32_ans    += step_ans_32
	print 'K_32 j =', j*step, (j+1)*step, np.abs(step_ans_32/K_32_ans), time_after_32 - time_before_32
        j += 1.

ans2 = 4. * np.pi * epsilon0 * omega * (K_12_prefactor(omega) * K_12_ans * np.cos(theta)
	                              + K_32_prefactor(omega) * K_32_ans * np.sin(theta))

print 'ans2/c                  ', ans2 / c


time_after = time.time()
print 'total time elapsed:', time_after - time_before, 'sec'
