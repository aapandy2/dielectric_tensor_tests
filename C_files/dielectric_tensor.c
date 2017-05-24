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
	p->theta_e = 0.5;         //dimensionless electron temp
	p->theta   = M_PI/3.;     //observer angle

	//derived quantities
	p->omega_p = sqrt(p->n_e * p->e*p->e / (p->m * p->epsilon0));       //plasma frequency    
        p->omega_c = p->e * p->B / (p->m * p->c);                        //cyclotron frequency

	//integrator parameters
	p->gamma             = 1.; //will get reset later in integration
	p->resolution_factor = 1;

	return 1;
}

int main(void)
{
        struct params p;
	set_params(&p);
	p.omega = 30. * p.omega_c;	

//        double prefactor = - 1. * omega_p*omega_p / (omega * omega)
//                            * 1./(4. * theta_e*theta_e * gsl_sf_bessel_Kn(2, 1./theta_e));

        printf("\n%e\n", start_search_12(&p));

}

