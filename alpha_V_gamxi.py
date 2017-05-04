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

#-----------------------------------------------------------------------------#
time_before = time.time()

#omega_mult = np.logspace(0., 3., 10) * omega_c
omega_mult = [1. * omega_c]
omega = omega_mult[0]

#print omega_mult / omega_c

for omega in omega_mult:
	K_12_ans  = fixed_quad(lambda tau_prime: np.vectorize(K_12_p_integrated)(tau_prime, omega), 0., 30., n = 45)[0]
	print 'p1 done'
#	K_12_ans += fixed_quad(lambda tau_prime: np.vectorize(K_12_p_integrated)(tau_prime, omega), 30., 60., n = 45)[0]
#	print 'p2 done'
#	K_12_ans += fixed_quad(lambda tau_prime: np.vectorize(K_12_p_integrated)(tau_prime, omega), 60., 90., n = 45)[0]
#	print 'p3 done'
#	K_12_ans += fixed_quad(lambda tau_prime: np.vectorize(K_12_p_integrated)(tau_prime, omega), 90., 120., n = 45)[0]
#	print 'p4 done'
#	K_12_ans += fixed_quad(lambda tau_prime: np.vectorize(K_12_p_integrated)(tau_prime, omega), 120., 150., n = 45)[0]
#	print 'p5 done'
#	K_12_ans += fixed_quad(lambda tau_prime: np.vectorize(K_12_p_integrated)(tau_prime, omega), 150., 180., n = 45)[0]
#	print 'p6 done'
#	K_12_ans += fixed_quad(lambda tau_prime: np.vectorize(K_12_p_integrated)(tau_prime, omega), 180., 210., n = 45)[0]
	print 'K_12 is done'

	K_32_ans  = fixed_quad(lambda tau_prime: np.vectorize(K_32_p_integrated)(tau_prime, omega), 0., 30., n = 45)[0]
	print 'p1 done'
#	K_32_ans += fixed_quad(lambda tau_prime: np.vectorize(K_32_p_integrated)(tau_prime, omega), 30., 60., n = 45)[0]
	print 'p2 done'
#	K_32_ans += fixed_quad(lambda tau_prime: np.vectorize(K_32_p_integrated)(tau_prime, omega), 60., 90., n = 45)[0]
#	print 'p3 done'
#	K_32_ans += fixed_quad(lambda tau_prime: np.vectorize(K_32_p_integrated)(tau_prime, omega), 90., 120., n = 45)[0]
#	print 'p4 done'
#	K_32_ans += fixed_quad(lambda tau_prime: np.vectorize(K_32_p_integrated)(tau_prime, omega), 120., 150., n = 45)[0]
#	print 'p5 done'
#	K_32_ans += fixed_quad(lambda tau_prime: np.vectorize(K_32_p_integrated)(tau_prime, omega), 150., 180., n = 45)[0]
#	print 'p6 done'
#	K_32_ans += fixed_quad(lambda tau_prime: np.vectorize(K_32_p_integrated)(tau_prime, omega), 180., 210., n = 45)[0]

	ans2 = 4. * np.pi * epsilon0 * omega * (K_12_prefactor(omega) * K_12_ans * np.cos(theta)
        	                              + K_32_prefactor(omega) * K_32_ans * np.sin(theta))

	print '-ans2/c                  ', - ans2 / c


time_after = time.time()
print 'time elapsed:', time_after - time_before, 'sec'


