#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "chi_spline_calc.h"

int set_params(struct params *p)
{
  p->epsilon0  = 1./(4. * M_PI); //permittivity of free space, CGS units
  p->e         = 4.80320680e-10; //electron charge
  p->m         = 9.1093826e-28;  //electron mass
  p->c         = 2.99792458e10;  //speed of light
  p->epsilon   = -1.;            //sign of electron charge
  
  //parameters
  p->B       = 1.;          //background B strength
  p->n_e     = 1.;          //electron number density cm^-3
  p->theta   = M_PI/3.;     //observer angle
  
  //derived quantities
  p->omega_p = sqrt(p->n_e*p->e*p->e / (p->m * p->epsilon0));//plasma frequency    
  p->omega_c = p->e * p->B / (p->m * p->c);               //cyclotron frequency
  
  //integrator parameters
  p->gamma             = 1.5; //will get reset later in integration
  p->real              = 1;   //real part = 1, imag part = 0
  
  //distribution function
  p->dist              = 0; //MJ=0, PL=1, kappa=2
  
  //distribution function parameters
  p->theta_e     = 10.;         //dimensionless electron temp
  p->pl_p        = 3.;          //power-law index, p
  p->gamma_min   = 1.;          //power-law gamma_min
  p->gamma_max   = 1000.;       //power-law gamma_max
  p->kappa       = 3.5;         //kappa index
  p->kappa_width = 10.;         //kappa width, like theta_e
  p->gamma_cutoff = 1e10;       //currently unused
  
  return 1;
}

int main(void)
{
  struct params p;
  
  /*set parameters*/
  set_params(&p);
  p.omega = 1. * p.omega_c;
  p.real  = 1;

//  printf("\n%e", chi_12(&p));
//  printf("\n%e", chi_32(&p));
//  printf("\n%e", chi_22(&p));
  printf("\n%e", alpha_V(&p));
  printf("\n");

  return 0;
}

/*alpha_I: returns the absorption coefficient alpha_I, for the total intensity
 *         of light along the ray in question, for the given values of 
 *         parameters within the struct p.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: absorption coefficient for total intensity (Stokes I) 
 */
double alpha_I(struct params *p)
{
  p->real          = 0;
  double prefactor = 2. * M_PI * p->epsilon0 * p->omega / p->c;
  double term11    = (chi_11(p) * pow(cos(p->theta), 2.)  
  		  + chi_33(p) * pow(sin(p->theta), 2.)
  		  - 2. * chi_13(p) * sin(p->theta) * cos(p->theta));
  double term22    = chi_22(p);
  double ans       = prefactor * (term11 + term22);
  return ans;
}

/*alpha_Q: returns the absorption coefficient alpha_Q, for linearly polarized
 *         light along the ray in question, for the given values of parameters 
 *         within the struct p.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: absorption coefficient for linearly polarized light (Stokes Q) 
 */
double alpha_Q(struct params *p)
{
  p->real          = 0;
  double prefactor = 2. * M_PI * p->epsilon0 * p->omega / p->c;
  double term11    = (chi_11(p) * pow(cos(p->theta), 2.)
                    + chi_33(p) * pow(sin(p->theta), 2.)
                    - 2. * chi_13(p) * sin(p->theta) * cos(p->theta));
  double term22    = chi_22(p);
  double ans       = prefactor * (term11 - term22);
  return ans;
}

/*rho_Q: returns the Faraday conversion coefficient rho_Q, which corresponds
 *       to the conversion between linearly polarized and circularly
 *       polarized light by the medium.  The coefficient is calculated for
 *       the given values of parameters within the struct p.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: Faraday conversion coefficient rho_Q 
 */
double rho_Q(struct params *p)
{
  p->real          = 1;
  double prefactor = 2. * M_PI * p->epsilon0 * p->omega / p->c;
  double term11    = (chi_11(p) * pow(cos(p->theta), 2.)
                    + chi_33(p) * pow(sin(p->theta), 2.)
                    - 2. * chi_13(p) * sin(p->theta) * cos(p->theta));
  double term22    = chi_22(p);
  double ans       = prefactor * (term22 - term11);
  return ans;
}

/*alpha_V: returns the absorption coefficient alpha_V, for the circularly
 *         polarized light along the ray in question, for the given values of 
 *         parameters within the struct p.  Uses the IEEE/IAU convention for
 *         the sign of Stokes V.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: absorption coefficient for circularly polarized light (Stokes V) 
 */
double alpha_V(struct params *p)
{
  p->real            = 1;
  double prefactor   = 4. * M_PI * p->epsilon0 * p->omega / p->c;
  double term1     = (chi_12(p) * cos(p->theta) - chi_32(p) * sin(p->theta));
  double ans       = prefactor * term1;
  return ans;
}

/*rho_V: returns the Faraday rotation coefficient rho_V, which rotates the
 *       plane of polarization (EVPA) for linearly polarized light,
 *       for the given values of parameters within the struct p.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: Faraday rotation coefficient rho_V 
 */
double rho_V(struct params *p)
{
  p->real          = 0;
  double prefactor = 4. * M_PI * p->epsilon0 * p->omega / p->c;
  double term1     = (chi_12(p) * cos(p->theta) - chi_32(p) * sin(p->theta));
  double ans       = prefactor * term1;
  return ans;
}

double end_approx(struct params * params)
{
  double end;
  
  double MJ_max        = 0.5 * (3. * params->theta_e 
                         + sqrt(4. + 9. * params->theta_e * params->theta_e));
  double PL_max_real   = sqrt((1. + params->pl_p)/params->pl_p);
  double PL_max_moving = 50./sqrt(params->omega/params->omega_c) 
                         + 9. * pow(params->omega/params->omega_c, 1./3.);
  double kappa_max     = (-3. + 3. * params->kappa_width * params->kappa 
  		          + sqrt(1. - 4. * params->kappa 
  		  	         - 18. * params->kappa_width * params->kappa 
                                 + 4. * pow(params->kappa, 2.) 
                                 + 9. * pow(params->kappa_width
                                            * params->kappa, 2.))) 
  		         / (2. * (params->kappa - 2.));
  
  if(params->dist == 0)
  {
    end = 7. * MJ_max;
  }
  else if(params->dist == 1)
  {
    end = PL_max_moving;
  }
  else if(params->dist == 2)
  {
    end = 7. * kappa_max;
  }
  else
  {
    printf("\ndistribution or real/imag is set incorrectly");
    return 1.;
  }
  
  return end;
}

/*MJ: returns the value of the Maxwell-Juettner (relativistic thermal)
 *    distribution function for the parameters in the struct p.
 *    TODO: I think a minus sign and other parts of Df are absorbed
 *          into the integrands themselves, and should be moved here.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: MJ(params->gamma, params->theta_e)
 */
double MJ(struct params * params)
{
  double ans = exp(-params->gamma/params->theta_e) 
  	   / (4. * M_PI * params->theta_e*params->theta_e 
    	      * gsl_sf_bessel_Kn(2, 1./params->theta_e));
  return ans;
}

/*PL: returns the value of the power-law distribution function for 
 *    the parameters in the struct p.
 *    TODO: I think a minus sign and other parts of Df are absorbed
 *          into the integrands themselves, and should be moved here.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: PL(params->gamma, params->power_law_p, params->gamma_min,
 *             params->gamma_max)
 */
double PL(struct params * params)
{
  if(params->gamma > params->gamma_max || params->gamma < params->gamma_min)
  {
  	return 0.;
  }
  
  double beta = sqrt(1. - 1./pow(params->gamma, 2.));
  
  double ans = (params->pl_p - 1.) 
             * (-1 + 2. * params->gamma * params->gamma + params->pl_p 
                                        * (params->gamma*params->gamma - 1.))
  	     / (4. * M_PI * (pow(params->gamma_min, -1. - params->pl_p) 
                - pow(params->gamma_max, -1. - params->pl_p))
  	        * beta * (params->gamma*params->gamma - 1.)) 
             * pow(params->gamma, -3. - params->pl_p);

  return ans;	
}

/*kappa_to_be_normalized: returns the value of the unnormalized relativistic
 *                        kappa distribution for the parameters in the struct p.
 *                        This function is integrated over gamma to determine
 *                        the normalization constant in normalize_f(), below.
 *    TODO: I think a minus sign and other parts of Df are absorbed
 *          into the integrands themselves, and should be moved here.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: kappa_to_be_normalized(gamma, void * parameters)
 */
double kappa_to_be_normalized(double gamma, void * parameters)
{
  struct params * params = (struct params*) parameters;
  
  double beta = sqrt(1. - 1./pow(gamma, 2.));
  
  double body = pow((1. + (gamma - 1.)/(params->kappa * params->kappa_width)), 
                    -1. - params->kappa);
  
  double d3p = 4. * M_PI * gamma*gamma * beta;
  
  double ans = body * d3p;
  
  return ans;
}

/*normalize_f: normalizes the distribution function using GSL's 
 *             QAGIU integrator.   
 *
 *@params: struct of parameters to pass to distribution function
 *@returns: 1 over the normalization constant for the chosen distribution
 */
double normalize_f(double (*distribution)(double, void *),
                   struct params * params
                  )
{
  /*set GSL QAGIU integrator parameters */
  double lower_bound = 1.;
  double absolute_error = 0.;
  double relative_error = 1e-8;
  int limit  = 1000;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
  double result, error;
  gsl_function F;
  
  F.function = distribution;
  
  F.params = params;
  
  gsl_integration_qagiu(&F, lower_bound, absolute_error, 
                        relative_error, limit, w, &result, &error
                       );
  
  
  gsl_integration_workspace_free(w);
  
  return result;
}

/*kappa: returns the value of the relativistic kappa distribution function for 
 *    the parameters in the struct p, with the normalization handled
 *    numerically by normalize_f().
 *    TODO: I think a minus sign and other parts of Df are absorbed
 *          into the integrands themselves, and should be moved here.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: PL(params->gamma, params->power_law_p, params->gamma_min,
 *             params->gamma_max)
 */
double kappa(struct params * params)
{
  static double norm                  = 0.;
  static double previous_kappa        = 0.;
  static double previous_kappa_width  = 0.;
  static double previous_gamma_cutoff = 0.;
  if(norm == 0. || previous_kappa_width != params->kappa_width
                || previous_kappa       != params->kappa)
  {
    norm                  = 1./normalize_f(&kappa_to_be_normalized, params);
    previous_kappa        = params->kappa;
    previous_kappa_width  = params->kappa_width;
    previous_gamma_cutoff = params->gamma_cutoff;
  }
  
  double beta = sqrt(1. - 1./pow(params->gamma, 2.));
  
  double body = -pow((1. + (params->gamma - 1.)
                      /(params->kappa * params->kappa_width)), 
                    -2. - params->kappa)
  	        *(-1. - params->kappa) / (params->kappa_width * params->kappa);
  
  double ans = norm * body;
  
  return ans;
}

/*Df: chooses the distribution function given by the parameter params->dist,
 *    and returns that distribution function evaluated with the given
 *    parameters in the struct params.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: MJ(params), PL(params), or kappa(params), depending on params->dist
 */
double Df(struct params * params)
{
  if(params->dist == 0)
  {
  	return MJ(params);
  }
  else if(params->dist == 1)
  {
  	return PL(params);
  }
  else if(params->dist == 2)
  {
  	return kappa(params);
  }
  
  return 0.;
}

double chi_11(struct params * params)
{
  params->component = 11;
  return chi_ij(params);
}

double chi_12(struct params * params)
{
  params->component = 12;
  return chi_ij(params);
}

double chi_13(struct params * params)
{
  params->component = 13;
  return chi_ij(params);
}

double chi_22(struct params * params)
{
  params->component = 22;
  return chi_ij(params);
}

double chi_32(struct params * params)
{
  params->component = 32;
  return chi_ij(params);
}

double chi_33(struct params * params)
{
  params->component = 33;
  return chi_ij(params);
}

double chi_ij(struct params * params)
{
  double prefactor = 2. * M_PI * params->omega_p*params->omega_p 
                     / (params->omega * params->omega);
  double start = 1.;
  double end   = end_approx(params);

  double ans = prefactor * gauss_legendre(1., 150., params);

  return ans;
}

double chi_ij_integrand(double gamma, struct params * params)
{
  params->gamma = gamma;
  double beta = sqrt(1. - pow(gamma, -2.));

  double dist = Df(params);
  double gam_term = dist * pow(gamma, 3.) * pow(beta, 3.);
  double ans = gam_term * spline_integrand(gamma, 
                                           params->omega/params->omega_c, 
                                           params);

  return ans;
}

#define NUM_QUAD 21
double gauss_legendre(double start, double end, struct params * params)
{
  double quadPts[NUM_QUAD] = \
        {-9.93752171e-01,  -9.67226839e-01,  -9.20099334e-01,
         -8.53363365e-01,  -7.68439963e-01,  -6.67138804e-01,
         -5.51618836e-01,  -4.24342120e-01,  -2.88021317e-01,
         -1.45561854e-01,   1.98918497e-16,   1.45561854e-01,
          2.88021317e-01,   4.24342120e-01,   5.51618836e-01,
          6.67138804e-01,   7.68439963e-01,   8.53363365e-01,
          9.20099334e-01,   9.67226839e-01,   9.93752171e-01};

  double weights[NUM_QUAD] = \
        {0.01601723,  0.03695379,  0.05713443,  0.07610011,  0.09344442,
         0.1087973 ,  0.12183142,  0.13226894,  0.13988739,  0.1445244 ,
         0.14608113,  0.1445244 ,  0.13988739,  0.13226894,  0.12183142,
         0.1087973 ,  0.09344442,  0.07610011,  0.05713443,  0.03695379,
         0.01601723};

  /*first we change integration variables to x = 1/(gamma^2 - 1), where the
    upper integration bound b = 1 and the lower bound a = 0.  We then apply
    the transformation: 
    \int_a^b f(x) dx = (b-a)/2 \int_{-1}^1 f((b-a)x/2 + (a+b)/2) dx
                     =     1/2 \int_{-1}^1 f(     x/2 +     1/2) dx */
  double x      = 0.;
  int i         = 0;
  double weight = 0.;
  double sum    = 0.;
  int n = NUM_QUAD;

//  # pragma omp parallel for private(i , x, weight) reduction ( + : sum )
  for(i = 0; i < n; i++)
  {
    x        = quadPts[i];
    weight   = weights[i];

    sum = sum + (end - start)/2. 
                * chi_ij_integrand((end - start)/2. * x + (end + start)/2., params) * weight;
  }

  return sum;
}

#define GAM_ARRAY_SIZE 200
#define OM_ARRAY_SIZE 200
double spline_integrand(double gamma, double omratio, struct params * params)
{
  static int loaded_file = 0;

  //TODO: do this for all possible changeable parameters
  static int past_component = 0;

  if(past_component == 0 ||
     past_component != params->component)
  {
    loaded_file = 0;
    past_component = params->component;
  }

  double myvariable;
  int i;
  int j;

  double gamstart = 1.;
  double gamend   = 1000.;
  double omstart  = 1.;
  double omend    = 1000.;
  double gamstep  = 5.;
  double omstep   = 5.;

  static double gamma_vals[GAM_ARRAY_SIZE];
  static double om_vals[OM_ARRAY_SIZE];
  static double z_vals[GAM_ARRAY_SIZE * OM_ARRAY_SIZE];

  char fname[22]; //= "chi_12_real_step_n.txt";

  /*choose file*/
  if(params->component == 11)
  {
    if(params->real == 1)
    {
      strcpy(fname, "chi_11_real_step_5.txt");
    }
    else
    {
      strcpy(fname, "chi_11_imag_step_5.txt");
    }
  }
  else if(params->component == 12)
  {
    if(params->real == 1)
    {
      strcpy(fname, "chi_12_real_step_5.txt");
    }
    else
    {
      strcpy(fname, "chi_12_imag_step_5.txt");
    }
  }
  else if(params->component == 13)
  {
    if(params->real == 1)
    {
      strcpy(fname, "chi_13_real_step_5.txt");
    }
    else
    {
      strcpy(fname, "chi_13_imag_step_5.txt");
    }
  }
  else if(params->component == 22)
  {
    if(params->real == 1)
    {
      strcpy(fname, "chi_22_real_step_5.txt");
    }
    else
    {
      strcpy(fname, "chi_22_imag_step_5.txt");
    }
  }
  else if(params->component == 32)
  {
    if(params->real == 1)
    {
      strcpy(fname, "chi_32_real_step_5.txt");
    }
    else
    {
      strcpy(fname, "chi_32_imag_step_5.txt");
    }
  }
  else if(params->component == 33)
  {
    if(params->real == 1)
    {
      strcpy(fname, "chi_33_real_step_5.txt");
    }
    else
    {
      strcpy(fname, "chi_33_imag_step_5.txt");
    }
  }
  else
  {
    printf("\ncorresponding file not found\n");
  }

    if(loaded_file == 0)
  {
    for(i = 0; i < GAM_ARRAY_SIZE; i++)
    {
      gamma_vals[i] = gamstart + i*gamstep;
      om_vals[i]    = omstart  + i*omstep;
    }

    FILE *myfile;
    myfile=fopen(fname, "r");

    for(i = 0; i < GAM_ARRAY_SIZE; i++)
    {
      for (j = 0 ; j < OM_ARRAY_SIZE; j++)
      {
        fscanf(myfile,"%lf", &myvariable);
        z_vals[j*GAM_ARRAY_SIZE + i] = myvariable;
      }
    }

    fclose(myfile);

    loaded_file = 1;
  }


//  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  const gsl_interp2d_type *T = gsl_interp2d_bicubic;
  gsl_spline2d *spline = gsl_spline2d_alloc(T, GAM_ARRAY_SIZE, OM_ARRAY_SIZE);
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();

  /* initialize interpolation */
  gsl_spline2d_init(spline, gamma_vals, om_vals, z_vals, GAM_ARRAY_SIZE, OM_ARRAY_SIZE);

  double ans = gsl_spline2d_eval(spline, omratio, gamma, xacc, yacc);

//printf("\n%e\n", ans);

  gsl_spline2d_free(spline);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);

  return ans;
}
