#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include "dielectric_tensor.h"

double K_13_integrand(double tau_prime, void * parameters)
{
	struct params * params = (struct params*) parameters;

	double prefactor  = 1.; //should be 1j
	double beta       = sqrt(1. - pow(params->gamma, -2.));
	double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
	double delta      = 2. * params->omega/(params->epsilon * params->omega_c) 
			   * sin(params->theta) * params->gamma * beta 
			   * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

	double gamma_term = beta*beta * params->gamma * MJ(params);//* exp(-params->gamma/params->theta_e);
				//explicit imag part
//	double tau_term   = cos(tau_prime * params->gamma) * cos((params->epsilon * params->omega_c / params->omega) * tau_prime/2.);

//	double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//	double tau_term   = sin(tau_prime * params->gamma) 
//			    * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / 2.);

	/*K_13 is K_32 except -S_2 -> C_2*/
	double tau_term   = cos((params->epsilon * params->omega_c / params->omega) * tau_prime / 2.);
	double xi_term    = I_2_analytic(alpha, delta); 
	double ans        = prefactor * gamma_term * xi_term * tau_term * params->gamma*params->gamma * beta;

	return ans;
}

double tau_integrator_13(double gamma, void * parameters)
{
	struct params * params = (struct params*) parameters;

	/* does this only work for low theta_e? */
//	if(params->omega/params->omega_c < 10. && params->theta_e < 1.)
//	{
//		params->resolution_factor = 8;
//	}


        double ans_tot  = 0.;
	double ans_step = 0.;
	double error    = 0.;
        double step     = 5.;
        double start    = 0.;
        double end      = M_PI * params->omega / params->omega_c * 2. * params->resolution_factor;
	size_t n        = 50;
	size_t limit    = 5000;
	double epsabs   = 0.;
	double epsrel   = 1e-8;
	enum gsl_integration_qawo_enum gsl_weight;
        double sign_correction;

	//need to update value of gamma
	params-> gamma = gamma;

	/*set up GSL QAWO integrator.  Do we need a new table w every call to tau_integrator_12?*/
	/*we should also try QAWF; it fits the integrals we need, and may be faster than QAWO.  */

	if(params->real == 1)
        {
                gsl_weight      = GSL_INTEG_SINE;
                sign_correction = -1.;
        }
        else
        {
                gsl_weight      = GSL_INTEG_COSINE;
                sign_correction = 1.;
        }

	gsl_integration_qawo_table * table = 
				gsl_integration_qawo_table_alloc(gamma, step, gsl_weight, n);
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
	gsl_set_error_handler_off();
	gsl_function F;
	F.function = &K_13_integrand;
	F.params   = params;

        while(end - start >= step)
        {
	        gsl_integration_qawo(&F, start, epsabs, epsrel, limit, w, table, &ans_step, &error);
		ans_tot += ans_step;
                start   += step;
        }

	gsl_integration_qawo_table_set_length(table, end - start);
        gsl_integration_qawo(&F, start, epsabs, epsrel, limit, w, table, &ans_step, &error);
	ans_tot += ans_step; 

	gsl_integration_qawo_table_free(table);
	gsl_integration_workspace_free(w);

	return ans_tot * sign_correction;
}

double start_search_13(struct params * params)
{
	double tolerance = 0.1;
	double step      = 0.1;
	double gamma     = 1.0;
	double diff      = tolerance + 10.;

	/* describe this later */
	if(params->omega/params->omega_c < 10.)
	{
		return 1.;
	}

	double fac1 = 0.;
	double fac2 = 0.;
	while(diff > tolerance)
	{
		params->resolution_factor = 1;
		fac1 = tau_integrator_13(gamma, params); //need to set res_factor to 1 here
		params->resolution_factor = 2;
		fac2 = tau_integrator_13(gamma, params); //need to set res_factor to 2 here
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

double K_13(struct params * p)
{

//	double prefactor = p->omega_p*p->omega_p / (p->omega * p->omega) 
//			   * 1./(2. * p->theta_e*p->theta_e * gsl_sf_bessel_Kn(2, 1./p->theta_e));

	double prefactor = 2. * M_PI * p->omega_p*p->omega_p / (p->omega * p->omega);
	
	gsl_function F;
        F.function = &tau_integrator_13;
        F.params   = p;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(5000);

	double start  = start_search_13(p);
//	double end    = 150.;
	double ans    = 0.;
	double error  = 0.;
	size_t limit  = 50;
	double epsabs = 0.;
        double epsrel = 1e-8;

//        gsl_integration_qng(&F, start, end, epsabs, epsrel, &ans, &error, &limit);

	gsl_set_error_handler_off();
	gsl_integration_qagiu(&F, start, 0., 1e-8, limit, w, &ans, &error);
	gsl_integration_workspace_free(w);

	return prefactor * ans;
//	return ans;
}

double K_31(struct params * p)
{
	return K_13(p);
}