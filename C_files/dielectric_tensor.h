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

double I_1_of_2(double alpha, double delta);
double I_1_analytic(double alpha, double delta);
double I_2_analytic(double alpha, double delta);
double MJ(struct params * params);

double K_11_integrand(double tau_prime, void * parameters);
double K_13_integrand(double tau_prime, void * parameters);
double K_33_integrand(double tau_prime, void * parameters);

double tau_integrator_11(double gamma, void * parameters);
double tau_integrator_13(double gamma, void * parameters);
double tau_integrator_33(double gamma, void * parameters);
double tau_integrator_22(double gamma, void * parameters);

double start_search_13(struct params * params);
double start_search_11(struct params * params);

double K_11(struct params * p);
double K_22(struct params * p);
double K_33(struct params * p);
double K_12(struct params * p);
double K_21(struct params * p);
double K_23(struct params * p);
double K_32(struct params * p);
double K_13(struct params * p);
double K_31(struct params * p);
