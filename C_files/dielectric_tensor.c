#include <stdio.h>
#include <math.h>
#include "dielectric_tensor.h"

int set_params(struct params *p)
{
	p->epsilon0  = 1./(4. * M_PI);
	p->e         = 4.80320680e-10;
	p->m         = 9.1093826e-28;
	p->c         = 2.99792458e10;
	p->epsilon   = -1.;        //sign of electron charge
	
	//parameters
	p->B       = 1.;          //background B strength
	p->n_e     = 1.;          //electron number density cm^-3
	p->theta_e = 0.5;         //dimensionless electron temp
	p->theta   = M_PI/3.;     //observer angle

	//derived quantities
	p->omega_p = sqrt(p->n_e * p->e*p->e / (p->m * p->epsilon0));       //plasma frequency    
        p->omega_c = p->e * p->B / (p->m * p->c);                        //cyclotron frequency

	//integrator parameters
	p->gamma             = 1.5; //will get reset later in integration
	p->resolution_factor = 1;
	p->real              = 1;

	return 1;
}

double alpha_I(struct params *p)
{
	p->real          = 0;
        double prefactor = 2. * M_PI * p->epsilon0 * p->omega / p->c;
        double term11    = (K_11(p) * pow(cos(p->theta), 2.)  
			  + K_33(p) * pow(sin(p->theta), 2.)
			  + 2. * K_13(p) * sin(p->theta) * cos(p->theta));
	double term22    = K_22(p);
        double ans       = prefactor * (term11 + term22);
        return ans;
}

double alpha_Q(struct params *p)
{
        p->real          = 0;
        double prefactor = 2. * M_PI * p->epsilon0 * p->omega / p->c;
        double term11    = (K_11(p) * pow(cos(p->theta), 2.)
                          + K_33(p) * pow(sin(p->theta), 2.)
                          + 2. * K_13(p) * sin(p->theta) * cos(p->theta));
        double term22    = K_22(p);
        double ans       = prefactor * (term11 - term22);
        return ans;
}

double alpha_V(struct params *p)
{
	p->real          = 1;
	double prefactor = 4. * M_PI * p->epsilon0 * p->omega / p->c;
	double term1     = (K_12(p) * cos(p->theta) + K_32(p) * sin(p->theta));
	double ans       = prefactor * term1;
	return ans;
}

int main(void)
{
        struct params p;
	set_params(&p);
	p.omega = 5. * p.omega_c;
	p.gamma = 1.5;
	p.real  = 0;
//	double i = 1.;
//	while(i < 10)
//	{
//		printf("\n%e	%e", i, tau_integrator_33(i, &p));
//		i = i + 0.05;
//	}
//	printf("\n");
//
//	printf("\n%e\n", tau_integrator_33(2.01, &p));
//	printf("\nK_11:	%e\n", K_11(&p));
//	printf("\nK_12: %e\n", K_12(&p));
//	printf("\nK_13: %e\n", K_13(&p));
//	printf("\nK_21: %e\n", K_21(&p));
//	printf("\nK_22: %e\n", K_22(&p));
//	printf("\nK_23: %e\n", K_23(&p));
//	printf("\nK_31: %e\n", K_31(&p));
//	printf("\nK_32: %e\n", K_32(&p));
//	printf("\nK_33: %e\n", K_33(&p));
//	printf("\n%e	%e\n", 	start_search_12(&p), start_search_32(&p));	
	printf("\n%e	%e\n", p.omega/p.omega_c, alpha_I(&p));

}

