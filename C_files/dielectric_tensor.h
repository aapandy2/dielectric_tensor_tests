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
};

double K_12(struct params * p);
double K_32(struct params * p);
double start_search_12(struct params * p);
double start_search_32(struct params * p);
double tau_integrator_12(double tau_prime, void * parameters);
double tau_integrator_32(double tau_prime, void * parameters);
