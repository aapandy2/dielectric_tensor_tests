#include <stdio.h>
#include <math.h>
#include "dielectric_tensor.h"

int set_params(struct params *p)
{
	p->epsilon0  = 1./(4. * M_PI);
	p->e         = 4.80320680e-10;
	p->m         = 9.1093826e-28;
	p->c        = 2.99792458e10;
	p->epsilon  = -1.;        //sign of electron charge
	
	//parameters
	p->B       = 1.;          //background B strength
	p->n_e     = 1.;          //electron number density cm^-3
	p->theta_e = 10.;         //dimensionless electron temp
	p->theta   = M_PI/3.;     //observer angle

	//derived quantities
	p->omega_p = sqrt(p->n_e * p->e*p->e / (p->m * p->epsilon0));       //plasma frequency    
        p->omega_c = p->e * p->B / (p->m * p->c);                        //cyclotron frequency

	//integrator parameters
	p->gamma             = 1.; //will get reset later in integration
	p->resolution_factor = 1;

	return 1;
}

double alpha_V(struct params *p)
{
	double prefactor = 4. * M_PI * p->epsilon0 * p->omega / p->c;
	double term1     = (K_12(p) * cos(p->theta) + K_32(p) * sin(p->theta));
	double ans       = prefactor * term1;
	return ans;
}

int main(void)
{
        struct params p;
	set_params(&p);
	p.omega = 10. * p.omega_c;
//	double i = 1.;
//	while(i < 100)
//	{
//		printf("\n%e	%e", i, tau_integrator_32(i, &p));
//		i = i + 0.1;
//	}
//	printf("\n");

//	printf("\n%e\n", tau_integrator_12(10., &p));
//	printf("\n%e\n", tau_integrator_12(11., &p));
//	printf("\n%e	%e\n", 	start_search_12(&p), start_search_32(&p));	
	printf("\n%e	%e\n", p.omega/p.omega_c, alpha_V(&p));

}

