#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

//constants -- GLOBAL VARIABLES, for now
double epsilon0  = 1./(4. * M_PI);
double e         = 4.80320680e-10;
double m         = 9.1093826e-28;
double c        = 2.99792458e10;
double epsilon  = -1.;        //sign of electron charge

//parameters
double B       = 1.;          //background B strength
double n_e     = 1.;          //electron number density cm^-3
double theta_e = 0.5;         //dimensionless electron temp
double theta   = M_PI/3.;     //observer angle


struct params
{
	double omega_c;
	double omega;
	double gamma;
	int resolution_factor;
};


double I_1_analytic(double alpha, double delta)
{
	double A     = sqrt(alpha*alpha + delta*delta);
	double minus = sqrt(A - alpha);
	double plus  = sqrt(A + alpha);

	if(alpha == 0. || delta == 0. || minus == 0.)
	{
		return 0.;
	}

	double ds    = -delta / fabs(delta);
	double outer = 2. * pow(delta, 6.) / (pow(alpha, 3.) * pow(A, 5.) * pow(minus, 7.) * pow(plus, 7.));
	double inner = ((ds*2. * alpha * pow(A, 5.) * (delta + ds*minus * plus)) * cos(alpha)
                     + pow(alpha, 3.) * delta * ds*(2. * pow(alpha, 2.) - pow(delta, 2.)) * A * cos(A)
                     + ds*2. * pow(alpha, 4.) * delta * A * sin(alpha)
                     + ds*4. * pow(alpha, 2.) * pow(delta, 3.) * A * sin(alpha)
                     + ds*2. * pow(delta, 5.) * A * sin(alpha)
                     + 2. * pow(alpha, 4.) * A * minus * plus * sin(alpha)
                     + 4. * pow(alpha, 2.) * pow(delta, 2.) * A * minus * plus * sin(alpha) 
                     + 2. * pow(delta, 4.) * A * minus * plus * sin(alpha)
                     - ds*pow(alpha, 3.) * delta * (2. * pow(alpha, 2.) 
						    + (-1. + pow(alpha, 2.)) * pow(delta, 2.) 
						    + pow(delta, 4.)) * sin(A));
	double ans = outer * inner;
	return ans;
}

double K_12_integrand_real(double tau_prime, void * parameters)
{
	struct params * params = (struct params*) parameters;

	double prefactor  = 1.; //should be 1j
	double beta       = sqrt(1. - pow(params->gamma, -2.));
	double alpha      = beta * cos(theta) * tau_prime * params->gamma;
	double delta      = 2. * params->omega/(epsilon * params->omega_c) * sin(theta) * params->gamma * beta 
			   * sin((epsilon * params->omega_c / params->omega) * tau_prime / (2.));

	double gamma_term = beta*beta * params->gamma * exp(-params->gamma/theta_e);
//	double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//	double tau_term   = -sin(tau_prime * params->gamma) 
//			    * sin((epsilon * params->omega_c / params->omega) * tau_prime);
	double tau_term   = -sin((epsilon * params->omega_c / params->omega) * tau_prime);
	double xi_term    = I_1_analytic(alpha, delta);
	double ans        = prefactor * gamma_term * xi_term * tau_term * params->gamma*params->gamma * beta;
	
//	printf("K_12_value	%e	%e	%e\n", ans, params->omega/params->omega_c, params->gamma);
	return ans;
}

double K_12_integrand_imag(double tau_prime, struct params* params)
{
        double prefactor  = 1.; //should be 1j
        double beta       = sqrt(1. - pow(params->gamma, -2.));
        double alpha      = beta * cos(theta) * tau_prime * params->gamma;
        double delta      = 2. * params->omega/(epsilon * params->omega_c) * sin(theta) * params->gamma * beta
                           * sin((epsilon * params->omega_c / params->omega) * tau_prime / (2.));

        double gamma_term = beta*beta * params->gamma * exp(-params->gamma/theta_e);
//      double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
        double tau_term   = cos(tau_prime * params->gamma) 
			    * sin((epsilon * params->omega_c / params->omega) * tau_prime);
        double xi_term    = I_1_analytic(alpha, delta);
        double ans        = prefactor * gamma_term * xi_term * tau_term * params->gamma*params->gamma * beta;
        return ans;
}

/*
def tau_first_12_mod(gamma, omega, factor):
            #hardcode only to second zero
            ans   = 0.
            step  = 30.
            start = 0.
            end   = K_12_zero(omega) * 2. * factor #specific to K_12, should be 2. *
            while(end - start > step):
                ans   += quad(lambda tau: K_12_xi_integrated_real(gamma, tau, omega).real, 
                              start, start + step, weight='sin', wvar=gamma, epsabs=0.)[0]
                start += step
            ans += quad(lambda tau: K_12_xi_integrated_real(gamma, tau, omega) * np.sin(gamma * tau), 
                        start, end, epsabs=0.)[0]
            return ans
*/

double tau_integrator_12(double gamma, struct params* parameters)
{
	struct params * params = (struct params*) parameters;

        double ans_tot  = 0.;
	double ans_step = 0.;
	double error    = 0.;
        double step     = 30.;
        double start    = 0.;
        double end      = M_PI * params->omega / params->omega_c * 2. * params->resolution_factor;
	size_t n        = 50;

	/*set up GSL QAWO integrator.  Do we need a new table w every call to tau_integrator_12?*/
	/*we should also try QAWF; it fits the integrals we need, and may be faster than QAWO.  */
	gsl_integration_qawo_table * table = 
				gsl_integration_qawo_table_alloc(params->gamma, step, GSL_INTEG_SINE, n);
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
	//gsl_set_error_handler_off();
	gsl_function F;
	F.function = &K_12_integrand_real;
	F.params   = params;

        while(end - start >= step)
        {
                gsl_integration_qawo(&F, start, 0., 1e-8, 50, w, table, &ans_step, &error);
		ans_tot += ans_step;
                start   += step;
        }

	gsl_integration_qawo_table_set_length(table, end - start);
        gsl_integration_qawo(&F, start, 0., 1e-8, 50, w, table, &ans_step, &error); //INTEGRATE END STEP HERE
	ans_tot += ans_step; 

	gsl_integration_qawo_table_free(table);
	gsl_integration_workspace_free(w);

	return ans_tot;
}


double start_search_12(struct params * params)
{
	double tolerance = 0.1;
	double step      = 0.1;
	double gamma     = 1.;
	double diff      = tolerance + 10.;

	double fac1 = 0.;
	double fac2 = 0.;
	while(diff > tolerance)
	{
		params->resolution_factor = 1;
		fac1 = tau_integrator_12(gamma, params); //need to set res_factor to 1 here
		params->resolution_factor = 2;
		fac2 = tau_integrator_12(gamma, params); //need to set res_factor to 2 here
		if(fac1 != 0. && fac2 != 0.)
		{
			diff = fabs((fac2 - fac1)/fac2); 
		}
		gamma += step;
	}
	//last iteration of while loop takes 1 step after finding correct answer
	//so we subtract that off
	return (gamma - step);
}


int main(void)
{ 
	//derived quantities
	double omega_p = sqrt(n_e * e*e / (m * epsilon0));       //plasma frequency    
	double omega_c = e * B / (m * c);                        //cyclotron frequency

	double omega = 1. * omega_c;

	struct params p;
	p.gamma             = 1.5;
	p.omega             = omega;
	p.omega_c           = omega_c;
	p.resolution_factor = 1;

	printf("\n%e\n", tau_integrator_12(1.5, &p));
	
	
	double prefactor = - 1. * omega_p*omega_p / (omega * omega) 
			    * 1./(4. * theta_e*theta_e * gsl_sf_bessel_Kn(2, 1./theta_e));

	double K_12_zero = M_PI * omega / omega_c;
	/*
	def start_search_12(omega):
	    tolerance = 0.1
	    step      = 0.1
	    gamma     = 1.
	    diff      = tolerance + 10.
	
	    while(diff > tolerance):
	        fac1 = tau_first_12_mod(gamma, omega, 1)
	        fac2 = tau_first_12_mod(gamma, omega, 2)
	        if(fac1 != 0 and fac2 != 0):
	            diff = np.abs((fac2 - fac1)/fac2)
	        gamma += step
	        
	    return gamma - step
	
	def tau_first_12_mod(gamma, omega, factor):
	    #hardcode only to second zero
	    ans   = 0.
	    step  = 30.
	    start = 0.
	    end   = K_12_zero(omega) * 2. * factor #specific to K_12, should be 2. *
	    while(end - start > step):
	        ans   += quad(lambda tau: K_12_xi_integrated_real(gamma, tau, omega).real, 
	                      start, start + step, weight='sin', wvar=gamma, epsabs=0.)[0]
	        start += step
	    ans += quad(lambda tau: K_12_xi_integrated_real(gamma, tau, omega) * np.sin(gamma * tau), 
	                start, end, epsabs=0.)[0]
	    return ans
	
	K_12 = quad(lambda gamma: np.vectorize(tau_first_12_mod)(gamma, omega, 1), 
		    start_search_12(omega), np.inf)[0]
	*/
}


//int
//main (void)
//{
//  gsl_integration_workspace * w 
//    = gsl_integration_workspace_alloc (1000);
//
//  double result, error;
//  double lower = 0.;             //lower bound of integral
//  double upper = 10.;            //upper bound of integral
//  double alpha = 1.;             //oscillation frequency of sine/cosine
//  double length = upper - lower; //length of integral interval
//  int intervals = 25;
//
//  gsl_integration_qawo_table * QAWO_table
//    = gsl_integration_qawo_table_alloc (alpha, length, 
//                                        GSL_INTEG_SINE, intervals);
//
//  /*I think integral of f(x) = 1 with QAWO and parameter GSL_INTEG_SINE
//    will be the integral of sin(alpha*x), written analytically below*/
//  double expected = -cos(alpha*upper) / alpha + cos(alpha*lower) / alpha;
//
//  gsl_function F;
//  F.function = &f;
//  F.params = &alpha;
//
//  gsl_integration_qawo(&F, lower, 0, 1e-7, 1000, w, QAWO_table, 
//                       &result, &error);
//
//
//  printf ("result          = % .18f\n", result);
//  printf ("exact result    = % .18f\n", expected);
//  printf ("estimated error = % .18f\n", error);
//  printf ("actual error    = % .18f\n", result - expected);
//  printf ("intervals       = %zu\n", w->size);
//
//
//  gsl_integration_qawo_table_free (QAWO_table);
//
//  gsl_integration_workspace_free (w);
//
//  return 0;
//}
