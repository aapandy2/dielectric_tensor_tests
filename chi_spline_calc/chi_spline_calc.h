struct params
{
        double epsilon0;
	double epsilon;
	double e;
	double m;
	double c;
	double B;
	double n_e;
	double theta;
	double theta_e;
	double pl_p;
	double gamma_min;
	double gamma_max;
	double kappa;
	double kappa_width;
	double gamma_cutoff;
	double omega_c;
	double omega_p;
        double omega;
        double gamma;
	int real;
	int dist;
	double (*tau_integrand)(double, void * parameters);
	double (*gamma_integrand)(double, void * parameters);

        int component;
};

double spline_integrand(double gamma, double omratio, struct params * params);
double chi_ij_spline_integrand(double gamma, struct params * params);
double gauss_legendre(double start, double end, struct params * params);

double chi_11(struct params * params);
double chi_12(struct params * params);
double chi_13(struct params * params);
double chi_22(struct params * params);
double chi_32(struct params * params);
double chi_33(struct params * params);

double chi_ij(struct params * params);
double end_approx(struct params * params);
double MJ(struct params * params);
double PL(struct params * params);
double kappa(struct params * params);
double kappa_to_be_normalized(double gamma, void * parameters);
double Df(struct params * params);
double normalize_f(double (*distribution)(double, void *), struct params * params);

double alpha_I(struct params *p);
double alpha_Q(struct params *p);
double alpha_V(struct params *p);
double rho_Q(struct params *p);
double rho_V(struct params *p);

