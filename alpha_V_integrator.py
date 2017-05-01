
# coding: utf-8

# In[1]:

#get_ipython().magic(u'matplotlib inline')
import numpy as np
import pylab as pl
import time
import scipy.special as special
from scipy.integrate import quad, dblquad, fixed_quad

# Set plot parameters to make beautiful plots
pl.rcParams['figure.figsize']  = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 20  
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex']     = True
pl.rcParams['axes.linewidth']  = 1.5
pl.rcParams['axes.titlesize']  = 'medium'
pl.rcParams['axes.labelsize']  = 'large'

pl.rcParams['xtick.major.size'] = 8     
pl.rcParams['xtick.minor.size'] = 4     
pl.rcParams['xtick.major.pad']  = 8     
pl.rcParams['xtick.minor.pad']  = 8     
pl.rcParams['xtick.color']      = 'k'     
pl.rcParams['xtick.labelsize']  = 'large'
pl.rcParams['xtick.direction']  = 'in'    

pl.rcParams['ytick.major.size'] = 8     
pl.rcParams['ytick.minor.size'] = 4     
pl.rcParams['ytick.major.pad']  = 8     
pl.rcParams['ytick.minor.pad']  = 8     
pl.rcParams['ytick.color']      = 'k'     
pl.rcParams['ytick.labelsize']  = 'large'
pl.rcParams['ytick.direction']  = 'in'


# In[52]:

# EVERYTHING IS 1 units

#constants
#e        = 1.     #electron charge
#m        = 1.     #electron mass
#c        = 1.     #speed of light
#epsilon0 = 1.     #permittivity of free space
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
    prefactor = - 1. * omega_p**2. / (omega * epsilon * omega_c) * 1./(4. * theta_e**2. * special.kn(2, 1./theta_e))
    return prefactor

def K_12_integrand(p_perp_bar, p_z_bar, tau_prime, omega):    
    prefactor = 1j
    k_perp    = omega / c * np.sin(theta)
    k_z       = omega / c * np.cos(theta)
    gamma     = np.sqrt(1. + p_perp_bar**2. + p_z_bar**2.)
    beta_perp = p_perp_bar / gamma
    beta_z    = p_z_bar    / gamma
    b         = omega/(epsilon * omega_c) * np.sin(theta) * gamma * beta_perp
    
    term1 = p_perp_bar**3. / gamma**2. * np.exp(-gamma/theta_e)
    term2 = np.exp(1j * (1. - beta_z * np.cos(theta)) * omega / (epsilon * omega_c) * tau_prime)
    term3 = np.sin(tau_prime / gamma) * special.j0(2. * b * np.sin(tau_prime / (2. * gamma)))
    ans   = prefactor * term1 * term2 * term3
    return ans

def K_32_prefactor(omega):
    prefactor = omega_p**2. / (omega * epsilon * omega_c) * 1./(2. * theta_e**2. * special.kn(2, 1./theta_e))
    return prefactor

def K_32_integrand(p_perp_bar, p_z_bar, tau_prime, omega):
    prefactor = 1.
    k_perp    = omega / c * np.sin(theta)
    k_z       = omega / c * np.cos(theta)
    gamma     = np.sqrt(1. + p_perp_bar**2. + p_z_bar**2.)
    beta_perp = p_perp_bar / gamma
    beta_z    = p_z_bar    / gamma
    b         = omega/(epsilon * omega_c) * np.sin(theta) * gamma * beta_perp
    
    term1 = p_perp_bar**2. * p_z_bar / gamma**2. * np.exp(-gamma/theta_e)
    term2 = np.exp(1j * (1. - beta_z * np.cos(theta)) * omega / (omega_c * epsilon) * tau_prime)
    term3 = np.sin(tau_prime / (2. * gamma)) * special.jn(1, 2. * b * np.sin(tau_prime / (2. * gamma)))
    ans   = prefactor * term1 * term2 * term3
    return ans


# In[54]:

def K_12_p_integrated(tau_prime, omega):
    real_part = dblquad(lambda p_perp, p_z: K_12_integrand(p_perp, p_z, tau_prime, omega).real, 
                       -np.inf, np.inf, lambda x: 0., lambda x: np.inf)
    ans = real_part[0]
    return ans


def K_32_p_integrated(tau_prime, omega):
    real_part = dblquad(lambda p_perp, p_z: K_32_integrand(p_perp, p_z, tau_prime, omega).real, 
                       -np.inf, np.inf, lambda x: 0., lambda x: np.inf)
    ans = real_part[0]
    return ans

#-----------------------------------------------------------------------------#
time_before = time.time()

#omega_mult = np.linspace(1., 2., 10) * omega_c
omega_mult = [1. * omega_c]

for omega in omega_mult:
	K_12_ans  = fixed_quad(lambda tau_prime: np.vectorize(K_12_p_integrated)(tau_prime, omega), 0., 30., n = 45)[0]
	print 'p1 done'
	K_12_ans += fixed_quad(lambda tau_prime: np.vectorize(K_12_p_integrated)(tau_prime, omega), 30., 60., n = 45)[0]
	print 'p2 done'
	K_12_ans += fixed_quad(lambda tau_prime: np.vectorize(K_12_p_integrated)(tau_prime, omega), 60., 90., n = 45)[0]
	print 'p3 done'

	print 'K_12 is done'

	K_32_ans  = fixed_quad(lambda tau_prime: np.vectorize(K_32_p_integrated)(tau_prime, omega), 0., 30., n = 45)[0]
	print 'p1 done'
	K_32_ans += fixed_quad(lambda tau_prime: np.vectorize(K_32_p_integrated)(tau_prime, omega), 30., 60., n = 45)[0]
	print 'p2 done'
	K_32_ans += fixed_quad(lambda tau_prime: np.vectorize(K_32_p_integrated)(tau_prime, omega), 60., 90., n = 45)[0]
	print 'p3 done'

	ans1 = 4. * np.pi * epsilon0 * omega * (K_12_prefactor(omega_c) * K_12_ans * np.cos(theta) 
      	                               - K_32_prefactor(omega_c) * K_32_ans * np.sin(theta))

	ans2 = 4. * np.pi * epsilon0 * omega * (K_12_prefactor(omega_c) * K_12_ans * np.cos(theta)
        	                             + K_32_prefactor(omega_c) * K_32_ans * np.sin(theta))

	print '-plus * omega_c/omega / c', ans2 * omega_c/omega / c * -1.

print omega_mult / omega_c

time_after = time.time()
print 'time elapsed:', time_after - time_before, 'sec'


