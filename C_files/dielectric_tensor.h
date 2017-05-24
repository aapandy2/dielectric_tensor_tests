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
