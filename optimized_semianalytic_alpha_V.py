import numpy as np
import pylab as pl
import time
import scipy.special as special
from scipy.integrate import quad, dblquad, fixed_quad
from scipy.optimize import root

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

def K_12_prefactor(omega):
    prefactor = - 1. * omega_p**2. / (omega * omega) * 1./(4. * theta_e**2. * special.kn(2, 1./theta_e))
    return prefactor

def K_32_prefactor(omega):
    prefactor = omega_p**2. / (omega * omega) * 1./(2. * theta_e**2. * special.kn(2, 1./theta_e))
    return prefactor

def guess(n, tau_prime):
    if(tau_prime == 0):
	print 'TAU IS ZERO'
	return 0.
    ans = 1./np.sqrt(1. - (n * np.pi / (2. * tau_prime))**2.)
    return ans

def I_1_analytic(alpha, delta):
    A     = np.sqrt(alpha**2. + delta**2.)
    minus = np.sqrt(A - alpha)
    plus  = np.sqrt(A + alpha)

    if(alpha == 0. or delta == 0. or minus == 0. or plus == 0.):
	return 0. #TODO: find out why minus can equal zero when delta is nonzero

    ds    = -np.sign(delta)
    outer = 2. * delta**6. / (alpha**3. * A**5. * minus**7. * plus**7.)
    inner = ((ds*2. * alpha * A**5. * (delta + ds*minus * plus)) * np.cos(alpha)
             + alpha**3. * delta * ds*(2. * alpha**2. - delta**2.) * A * np.cos(A)
             + ds*2. * alpha**4. * delta * A * np.sin(alpha)
             + ds*4. * alpha**2. * delta**3. * A * np.sin(alpha)
             + ds*2. * delta**5. * A * np.sin(alpha)
             + 2. * alpha**4. * A * minus * plus * np.sin(alpha)
             + 4. * alpha**2. * delta**2. * A * minus * plus * np.sin(alpha)
             + 2. * delta**4. * A * minus * plus * np.sin(alpha)
             - ds*alpha**3. * delta * (2. * alpha**2. + (-1. + alpha**2.) * delta**2. + delta**4.) * np.sin(A))
    ans   = outer * inner
    return ans

def I_2_analytic(alpha, delta):
    A     = np.sqrt(alpha**2. + delta**2.)
    plus  = np.sqrt(A + alpha)
    minus = np.sqrt(A - alpha)

    if(alpha == 0. or delta == 0. or minus == 0. or plus == 0.):
	return 0.

    num   = alpha * delta**8. * (3. * A * np.cos(A) + (-3. + alpha**2. + delta**2.) * np.sin(A))
    denom = A**5. * minus**7. * plus**7.
    ans   = -num / denom
    return ans

def K_12_xi_integrateda(gamma, tau_prime, omega):    

    if(gamma < 1.):
	return 1000. #wall to stop root finder from trying gamma < 1

    prefactor  = 1j
    beta       = np.sqrt(1. - 1./gamma**2.)
    alpha      = beta * np.cos(theta) * tau_prime
    delta      = 2. * omega/(epsilon * omega_c) * np.sin(theta) * gamma * beta * np.sin(
                (epsilon * omega_c / omega) * tau_prime / (2. * gamma))
    
    gamma_term = beta**2. * np.exp(-gamma/theta_e) * np.sin((epsilon * omega_c / omega) * tau_prime / gamma)
    tau_term   = np.exp(1j * tau_prime)
    
    if(alpha == 0 or delta == 0):
        return 0.
    
    xi_term    = I_1_analytic(alpha, delta)
    ans        = prefactor * gamma_term * xi_term * tau_term
    return ans * gamma**2. * beta 

def K_32_xi_integrateda(gamma, tau_prime, omega):

    if(gamma < 1.):
        return 1000. #wall to stop root finder from trying gamma < 1

    prefactor = 1.
    beta       = np.sqrt(1. - 1./gamma**2.)
    alpha      = beta * np.cos(theta) * tau_prime
    delta      = 2. * omega/(epsilon * omega_c) * np.sin(theta) * gamma * beta * np.sin(
                (epsilon * omega_c / omega) * tau_prime / (2. * gamma))
    gamma_term = beta**2. * np.exp(-gamma/theta_e) * np.sin((epsilon * omega_c / omega) * tau_prime / (2. * gamma))
    tau_term   = np.exp(1j * tau_prime)
    
    if(alpha == 0 or delta == 0):
        return 0.
    
    xi_term    = -2j * I_2_analytic(alpha, delta)
    ans        = prefactor * gamma_term * tau_term * xi_term
    
    if(delta < 0):
        ans = ans * -1.
    
    return ans * gamma**2. * beta

def roots_32(tau_prime, omega):
    n_min = 1
    n_max = int(2. * tau_prime / np.pi)
    
    solutions = np.empty(0)
    
    
    for n in range(n_min, n_max, 2):
        rootfind = root(lambda gamma: K_32_xi_integrateda(gamma, tau_prime, omega).real, 
                        guess(n, tau_prime))
        solutions = np.append(solutions, rootfind.x[0])
    
    solutions = np.unique(solutions)
    
    tolerance = 1e-5
    duplicates = np.empty(0)
    for i in range(solutions.size - 1):
        if(np.abs(solutions[i+1] - solutions[i]) < tolerance):
            duplicates = np.append(duplicates, i)

    solutions = np.delete(solutions, duplicates)
    
    return solutions

def K_32_improved_gam(tau_prime, omega):
    solns = roots_32(tau_prime, omega)
    
    if(solns.size == 0):
        return quad(lambda gamma: K_32_xi_integrateda(gamma, tau_prime, omega).real, 1., np.inf)[0]
    
    ans = 0.
    
    first_part = quad(lambda gamma: K_32_xi_integrateda(gamma, tau_prime, omega).real, 1., solns[0])[0]
    #print 1., solns[0], solns[0] - 1.
    
    if(np.isfinite(first_part)):
        ans += first_part
        #print ans
    
    step = 0.
    for i in range(solns.size-1):
        step = quad(lambda gamma: K_32_xi_integrateda(gamma, tau_prime, omega).real, solns[i], solns[i+1])[0]
        #print 'ans = ', ans
        if(np.isfinite(step) == False):
            #print 'STEP WAS NOT FINITE', step
            step = 0.
        ans += step
    
    last_part = quad(lambda gamma: K_32_xi_integrateda(gamma, tau_prime, omega).real, 
                     solns[solns.size-1], np.inf)[0]
    if(np.isfinite(last_part)):
        ans += last_part
        #print 'last part', last_part 
    
    return ans


def roots_12(tau_prime, omega):
    n_min = 1
    n_max = int(2. * tau_prime / np.pi)
    
    solutions = np.empty(0)
    
    
    for n in range(n_min, n_max, 2):
        rootfind = root(lambda gamma: K_12_xi_integrateda(gamma, tau_prime, omega).real, 
                        guess(n, tau_prime))
        solutions = np.append(solutions, rootfind.x[0])
    
    solutions = np.unique(solutions)
    
    tolerance = 1e-5
    duplicates = np.empty(0)
    for i in range(solutions.size - 1):
        if(np.abs(solutions[i+1] - solutions[i]) < tolerance):
            duplicates = np.append(duplicates, i)

    solutions = np.delete(solutions, duplicates)
    
    return solutions

def K_12_improved_gam(tau_prime, omega):
    solns = roots_12(tau_prime, omega)
    
    if(solns.size == 0):
        return quad(lambda gamma: K_12_xi_integrateda(gamma, tau_prime, omega).real, 1., np.inf)[0]
    
    ans = 0.
    
    first_part = quad(lambda gamma: K_12_xi_integrateda(gamma, tau_prime, omega).real, 1., solns[0])[0]
    #print 1., solns[0], solns[0] - 1.
    
    if(np.isfinite(first_part)):
        ans += first_part
        #print ans
    
    step = 0.
    for i in range(solns.size-1):
        step = quad(lambda gamma: K_12_xi_integrateda(gamma, tau_prime, omega).real, solns[i], solns[i+1])[0]
        #print 'ans = ', ans
        if(np.isfinite(step) == False):
            #print 'STEP WAS NOT FINITE', step
            step = 0.
        ans += step
    
    last_part = quad(lambda gamma: K_12_xi_integrateda(gamma, tau_prime, omega).real, 
                     solns[solns.size-1], np.inf)[0]
    if(np.isfinite(last_part)):
        ans += last_part
        #print 'last part', last_part 
    
    return ans


def K_32_p_integrateda(tau_prime, omega):
    real_part = quad(lambda gamma: K_32_xi_integrateda(gamma, tau_prime, omega).real, 1., np.inf)
    ans       = real_part[0]
    return ans

def K_12_p_integrateda(tau_prime, omega):
    real_part = quad(lambda gamma: K_12_xi_integrateda(gamma, tau_prime, omega).real, 1., np.inf)
    ans       = real_part[0]
    return ans

#def alpha_V(omega, tau_max):
#	K_12 = quad(lambda tau: np.vectorize(K_12_p_integrateda)(tau, omega), 0., np.inf)[0]
#	K_32 = quad(lambda tau: np.vectorize(K_32_p_integrateda)(tau, omega), 0., np.inf)[0]
#	K_12 = quad(lambda tau: np.vectorize(K_12_improved_gam)(tau, omega), 0., tau_max)[0]
#        K_32 = quad(lambda tau: np.vectorize(K_32_improved_gam)(tau, omega), 0., tau_max)[0]
#	ans = (K_12 * K_12_prefactor(omega) * np.cos(theta) 
#		 + K_32 * K_32_prefactor(omega) * np.sin(theta)) / c * omega
#
#	return ans

def alpha_V(omega, num_steps, step):
	K_12 = 0.
	K_32 = 0.
	for i in range(num_steps):
		print 'step', i+1, '/', num_steps
#		if(i * step < 200):
		K_12 += quad(lambda tau: np.vectorize(K_12_improved_gam)(tau, omega), 
			     i*step, (i+1)*step)[0]
		K_32 += quad(lambda tau: np.vectorize(K_32_improved_gam)(tau, omega), 
			     i*step, (i+1)*step)[0]
#		else:
#			K_12 += quad(lambda tau: np.vectorize(K_12_p_integrated1)(tau, omega, 0.01, 15.), 
#                       	i*step, (i+1)*step)[0]
#                	K_32 += quad(lambda tau: np.vectorize(K_32_p_integrated1)(tau, omega, 0.01, 15), 
#                            	i*step, (i+1)*step)[0]

	ans = (K_12 * K_12_prefactor(omega) * np.cos(theta) 
		 + K_32 * K_32_prefactor(omega) * np.sin(theta)) / c * omega

	return ans

time_before = time.time()
print alpha_V(200. * omega_c, 11, 30)
time_after  = time.time()
print 'time elapsed: ', time_after - time_before

