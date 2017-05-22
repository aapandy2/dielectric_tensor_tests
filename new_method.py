import numpy as np

import matplotlib #prevent SSH from killing plots
matplotlib.use('Agg')

import pylab as pl
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


# In[10]:

def K_12_prefactor(omega):
    prefactor = - 1. * omega_p**2. / (omega * omega) * 1./(4. * theta_e**2. * special.kn(2, 1./theta_e))
    return prefactor

def K_32_prefactor(omega):
    prefactor = omega_p**2. / (omega * omega) * 1./(2. * theta_e**2. * special.kn(2, 1./theta_e))
    return prefactor

def K_12_xi_integrateda(gamma, tau_prime, omega):    
    prefactor  = 1j
    beta       = np.sqrt(1. - 1./gamma**2.)
    alpha      = beta * np.cos(theta) * tau_prime * gamma
    delta      = 2. * omega/(epsilon * omega_c) * np.sin(theta) * gamma * beta * np.sin(
                (epsilon * omega_c / omega) * tau_prime / (2.))
    
    gamma_term = beta**2. * gamma * np.exp(-gamma/theta_e) * np.sin((epsilon * omega_c / omega) * tau_prime)
    tau_term   = np.exp(1j * tau_prime * gamma)    
    xi_term    = I_1_analytic(alpha, delta)
    ans        = prefactor * gamma_term * xi_term * tau_term   
    return ans * gamma**2. * beta 

def I_1_analytic(alpha, delta):
    A     = np.sqrt(alpha**2. + delta**2.)
    minus = np.sqrt(A - alpha)
    plus  = np.sqrt(A + alpha)
    if(alpha == 0 or delta == 0 or minus == 0.):
        return 0.
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
    if(alpha == 0 or delta == 0):
        return 0.
    A     = np.sqrt(alpha**2. + delta**2.)
    plus  = np.sqrt(A + alpha)
    minus = np.sqrt(A - alpha)
    num   = alpha * delta**8. * (3. * A * np.cos(A) + (-3. + alpha**2. + delta**2.) * np.sin(A))
    denom = A**5. * minus**7. * plus**7.
    delta_correction = delta / np.abs(delta)
    ans   = -num / denom * delta_correction
    return ans

def K_32_xi_integrateda(gamma, tau_prime, omega):
    prefactor = 1.
    beta       = np.sqrt(1. - 1./gamma**2.)
    alpha      = beta * np.cos(theta) * tau_prime * gamma
    delta      = 2. * omega/(epsilon * omega_c) * np.sin(theta) * gamma * beta * np.sin(
                (epsilon * omega_c / omega) * tau_prime / (2.))
    gamma_term = beta**2. * gamma * np.exp(-gamma/theta_e) * np.sin((epsilon * omega_c / omega) 
                                                                    * tau_prime / (2.))
    tau_term   = np.exp(1j * tau_prime * gamma)
    xi_term    = -2j * I_2_analytic(alpha, delta)
    ans        = prefactor * gamma_term * tau_term * xi_term
    
    if(np.isnan(ans)):
        print 'NAN in K_32_xi_integrateda'
    
    return ans * gamma**2. * beta


def K_12_p_integrateda(tau_prime, omega):
    real_part = quad(lambda gamma: K_12_xi_integrateda(gamma, tau_prime, omega).real, 1., np.inf)
    ans       = real_part[0]
    return ans

def K_32_p_integrateda(tau_prime, omega):
    real_part = quad(lambda gamma: K_32_xi_integrateda(gamma, tau_prime, omega).real, 1., np.inf)
    ans       = real_part[0]
    return ans


def K_12_zero(omega):
    ans = np.pi * omega / omega_c
    return ans

def K_32_zero(omega):
    ans = 2. * np.pi * omega / omega_c
    return ans

def tau_first_12(gamma, omega):
    #hardcode only to second zero
    ans   = 0.
    step  = 30.
    start = 0.
    end   = K_12_zero(omega) * 2. #specific to K_12
    while(end - start > step):
        ans   += quad(lambda tau: K_12_xi_integrateda(gamma, tau, omega).real, start, start + step, epsabs=0.)[0]
        start += step
    ans += quad(lambda tau: K_12_xi_integrateda(gamma, tau, omega).real, start, end, epsabs=0.)[0]
    return ans

def tau_first_32(gamma, omega):
    #hardcode only to first zero
    ans   = 0.
    step  = 30.
    start = 0.
    end   = K_32_zero(omega) * 1. #specific to K_32
    while(end - start > step):
        ans   += quad(lambda tau: K_32_xi_integrateda(gamma, tau, omega).real, start, start + step, epsabs=0.)[0]
        start += step
    ans += quad(lambda tau: K_32_xi_integrateda(gamma, tau, omega).real, start, end, epsabs=0.)[0]
    return ans

def K_12_xi_integrated_real(gamma, tau_prime, omega):    
    prefactor  = 1.
    beta       = np.sqrt(1. - 1./gamma**2.)
    alpha      = beta * np.cos(theta) * tau_prime * gamma
    delta      = 2. * omega/(epsilon * omega_c) * np.sin(theta) * gamma * beta * np.sin(
                (epsilon * omega_c / omega) * tau_prime / (2.))
    
    gamma_term = beta**2. * gamma * np.exp(-gamma/theta_e) 
    #tau_term   = np.exp(1j * tau_prime * gamma) * np.sin((epsilon * omega_c / omega) * tau_prime)
    #tau_term   = -np.sin(tau_prime * gamma) * np.sin((epsilon * omega_c / omega) * tau_prime)
    tau_term   = -1. * np.sin((epsilon * omega_c / omega) * tau_prime)
    xi_term    = I_1_analytic(alpha, delta)
    ans        = prefactor * gamma_term * xi_term * tau_term   
    return ans * gamma**2. * beta 

def K_32_xi_integrated_real(gamma, tau_prime, omega):
    prefactor = 1.
    beta       = np.sqrt(1. - 1./gamma**2.)
    alpha      = beta * np.cos(theta) * tau_prime * gamma
    delta      = 2. * omega/(epsilon * omega_c) * np.sin(theta) * gamma * beta * np.sin(
                (epsilon * omega_c / omega) * tau_prime / (2.))
    gamma_term = beta**2. * gamma * np.exp(-gamma/theta_e) 
#    tau_term   = -np.sin(tau_prime * gamma) * np.sin((epsilon * omega_c / omega) * tau_prime / (2.))
    tau_term   = -1. * np.sin((epsilon * omega_c / omega) * tau_prime / (2.))
    xi_term    = -2. * I_2_analytic(alpha, delta)
    ans        = prefactor * gamma_term * tau_term * xi_term
    
    return ans * gamma**2. * beta

def tau_first_12_mod(gamma, omega):
    #hardcode only to second zero
    ans   = 0.
    step  = 60.
    start = 0.
    end   = K_12_zero(omega) * 2. #specific to K_12, should be 2. *
    while(end - start > step):
        ans   += quad(lambda tau: K_12_xi_integrated_real(gamma, tau, omega).real, 
                      start, start + step, weight='sin', wvar=gamma, epsabs=0.)[0]
        start += step
    ans += quad(lambda tau: K_12_xi_integrated_real(gamma, tau, omega) * np.sin(gamma * tau), 
                start, end, epsabs=0.)[0]
    return ans

def tau_first_32_mod(gamma, omega):
    #hardcode only to first zero
    ans   = 0.
    step  = 60.
    start = 0.
    end   = K_32_zero(omega) * 1. #specific to K_32, should be 1. *
    while(end - start > step):
        ans   += quad(lambda tau: K_32_xi_integrated_real(gamma, tau, omega).real, 
                      start, start + step, weight='sin', wvar=gamma, epsabs=0.)[0]
        start += step
    ans += quad(lambda tau: K_32_xi_integrated_real(gamma, tau, omega) * np.sin(gamma * tau), 
                start, end, epsabs=0.)[0]
    return ans


def alpha_V_new(omega):
    K_12 = fixed_quad(lambda gamma: np.vectorize(tau_first_12_mod)(gamma, omega), 1., 20., n=45)[0]
    K_32 = fixed_quad(lambda gamma: np.vectorize(tau_first_32_mod)(gamma, omega), 1., 20., n=45)[0]
    ans = (K_12 * K_12_prefactor(omega) * np.cos(theta) 
           + K_32 * K_32_prefactor(omega) * np.sin(theta)) * omega / c
    return ans


omega = 10. * omega_c
time_before = time.time()
print alpha_V_new(omega)
time_after  = time.time()
print 'time elapsed:', time_after - time_before
