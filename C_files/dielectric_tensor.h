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
	double omega_c;
	double omega_p;
        double omega;
        double gamma;
        int resolution_factor;
	int real;
};

double I_2_analytic(double alpha, double delta);
double MJ(struct params * params);


double K_13_integrand(double tau_prime, void * parameters);

double K_12(struct params * p);
double K_32(struct params * p);
double K_13(struct params * p);
