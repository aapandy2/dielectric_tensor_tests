#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include "dielectric_tensor.h"

double I_3_analytic(double alpha, double delta)
{
	if(alpha == 0. || delta == 0.)
	{
		return 0.;
	}	

	double A        = sqrt(alpha*alpha + delta*delta);
	double term1    = 6. * alpha*alpha * cos(A) / pow(A, 4.);
	double term2    = -2. * cos(A) / (A*A);
	double term3    = 6. * delta*delta * sin(A) / pow(A, 5.);
	double term4    = -4. * sin(A) / pow(A, 3.);
	double term5    = 2. * alpha*alpha * sin(A) / pow(A, 3.);
	double ans      = term1 + term2 + term3 + term4 + term5;
	return ans;
}

double I_3_limit(double alpha, double delta)
{
	if(alpha == 0.)
	{
		return 0.;
	}

        double A        = fabs(alpha);
        double term5    = 2. * alpha*alpha * sin(A) / pow(A, 3.);
        double ans      = term5;
        return ans;
}

double K_33_integrand(double tau_prime, void * parameters)
{
	struct params * params = (struct params*) parameters;

	double prefactor  = 1.; //should be 1j
	double beta       = sqrt(1. - pow(params->gamma, -2.));
	double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
	double delta      = 2. * params->omega/(params->epsilon * params->omega_c) 
			   * sin(params->theta) * params->gamma * beta 
			   * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

	double gamma_term = beta*beta * params->gamma * MJ(params);
//	double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//	double tau_term   = -sin(tau_prime * params->gamma) 
//			    * sin((epsilon * params->omega_c / params->omega) * tau_prime);
	double tau_term   = 1.;
	double xi_term    = I_3_analytic(alpha, delta) - I_3_limit(alpha, delta);
	double ans        = prefactor * gamma_term * xi_term * tau_term * params->gamma*params->gamma * beta;
	
	return ans;
}

/*note: answer is totally imaginary*/
double I_3_limit_integrated(struct params * p)
{
	double beta     = sqrt(1. - 1./pow(p->gamma, 2.));
	double ans      = 1. / (p->gamma * beta * cos(p->theta)) * 2. * atanh(beta * cos(p->theta));
	return ans;
}

/*note: answer is totally real*/
double limit_term_integrated(struct params * params)
{
        double prefactor  = 1.; //should be 1j
        double beta       = sqrt(1. - pow(params->gamma, -2.));

        double gamma_term = beta*beta * params->gamma * MJ(params);
        double tau_term   = 1.;
        double xi_term    = I_3_limit_integrated(params);
        double ans        = prefactor * gamma_term * xi_term * tau_term * params->gamma*params->gamma * beta;

        return -ans; //note: -1 factor comes from I_3_limit_integrated being imaginary, and the 1j prefactor here
}

//double tau_integrator_33(double gamma, void * parameters)
//{
//	struct params * params = (struct params*) parameters;
//
//	if(gamma == 1.)
//	{
//		return 0.;
//	}
//
//
//	/* does this only work for low theta_e? */
////	if(params->omega/params->omega_c < 10. && params-> theta_e < 1.)
////	{
////		params->resolution_factor = 8;
////	}
//
//
//        double ans_tot  = 0.;
//	double ans_step = 0.;
//	double error    = 0.;
//        double step     = 5.;
//        double start    = 0.;
//        double end      = M_PI * params->omega / params->omega_c * 2. * params->resolution_factor;
//	size_t n        = 50;
//	size_t limit    = 5000;
//	double epsabs   = 0.;
//	double epsrel   = 1e-8;
//	enum gsl_integration_qawo_enum gsl_weight;
//	double sign_correction;
//	//need to update value of gamma
//	params-> gamma = gamma;
//
//	/*set up GSL QAWO integrator.  Do we need a new table w every call to tau_integrator_12?*/
//	/*we should also try QAWF; it fits the integrals we need, and may be faster than QAWO.  */
//
//	if(params->real == 1)
//	{
//		gsl_weight      = GSL_INTEG_SINE;
//		sign_correction = -1.;
//	}
//	else
//	{
//		gsl_weight      = GSL_INTEG_COSINE;
//		sign_correction = 1.;
//	}
//
//	gsl_integration_qawo_table * table = 
//				gsl_integration_qawo_table_alloc(gamma, step, gsl_weight, n);
//	
//	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
//	gsl_set_error_handler_off();
//	gsl_function F;
//	F.function = &K_33_integrand;
//	F.params   = params;
//
//        while(end - start >= step)
//        {
//                gsl_integration_qawo(&F, start, epsabs, epsrel, limit, w, table, &ans_step, &error);
//		ans_tot += ans_step;
//                start   += step;
//        }
//
//	gsl_integration_qawo_table_set_length(table, end - start);
//        gsl_integration_qawo(&F, start, epsabs, epsrel, limit, w, table, &ans_step, &error);
//	ans_tot += ans_step; 
//
//	gsl_integration_qawo_table_free(table);
//	gsl_integration_workspace_free(w);
//
//	return ans_tot * sign_correction;
//}

double tau_integrator_33(double gamma, void * parameters)
{
	struct params * params = (struct params*) parameters;

	if(gamma == 1.)
	{
		return 0.;
	}


        double ans_tot  = 0.;
	double ans_step = 0.;
	double error    = 0.;
        double step     = 5.;
        double start    = 0.;
//        double end      = M_PI * params->omega / params->omega_c * 2. * params->resolution_factor;
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
	F.function = &K_33_integrand;
	F.params   = params;

	int i            = 0;
	int max_counter  = 200;
	double tolerance = 1e-3;
	int counts       = 0;

	int i_max        = 1000;
	double small_tol = 1e-8;
	while(i == 0 || counts < max_counter)
	{
		gsl_integration_qawo(&F, i*step, epsabs, epsrel, limit, w, table, &ans_step, &error);
		ans_tot += ans_step;
		i       += 1;

		if(fabs(ans_step / ans_tot) < tolerance)
		{
			counts += 1;
		}

		if(i >= i_max && fabs(ans_tot) < small_tol)
		{
			counts = max_counter;
		}
	}

//	printf("\n%e	%e\n", (i+1.)*step, 2. * M_PI * params->omega/params->omega_c);

	gsl_integration_qawo_table_free(table);
	gsl_integration_workspace_free(w);

	double I_3_limit_integral;
	if(params->real == 1)
	{
		I_3_limit_integral = limit_term_integrated(params);
	}
	else
	{
		I_3_limit_integral = 0.;
	}

	return ans_tot * sign_correction; //+ I_3_limit_integral;
}

double K_33(struct params * p)
{
//	double prefactor = 1. * p->omega_p*p->omega_p / (p->omega * p->omega)
//                           * 1./(2. * p->theta_e*p->theta_e * gsl_sf_bessel_Kn(2, 1./p->theta_e));

	double prefactor = 2. * M_PI * p->omega_p*p->omega_p / (p->omega * p->omega);
	
	gsl_function F;
        F.function = &tau_integrator_33;
        F.params   = p;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(5000);


	double start  = 1.; //start_search_12(p);
//	double end    = 150.; //figure out better way to do this
	double ans    = 0.;
	double error  = 0.;
	size_t limit  = 50;
	double epsabs = 0.;
	double epsrel = 1e-8;

	
	gsl_set_error_handler_off();

//	gsl_integration_qng(&F, start, end, epsabs, epsrel, &ans, &error, &limit);

	gsl_integration_qagiu(&F, start, 0., 1e-8, limit, w, &ans, &error);
	gsl_integration_workspace_free(w);

	return prefactor * ans;
//	return ans;
}
