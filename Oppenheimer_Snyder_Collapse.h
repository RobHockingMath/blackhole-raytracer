#ifndef OPPENHEIMER_SNYDER_COLLAPSE_H
#define OPPENHEIMER_SNYDER_COLLAPSE_H

#include <vector>
#include <string>
#include "Spacetime.h"
#include <math.h>

using namespace std;

class Oppenheimer_Snyder_Collapse : public Spacetime {
public:

    Oppenheimer_Snyder_Collapse(double pos_[3], double tol_, double mass, double R0);
	void integrate_a(double eps,double h);
	double get_a(double tau);
	double get_a_old(double tau);
	
	double tau_of_eta(double eta);
	double a_of_eta(double eta);
	double eta_from_tau(double tau);
	double get_a_old2(double tau);

	double eta_from_tau_newton(double tau);
	
	double compute_energy(double p_loc[]);

    bool accel_1st_order(double*& yp, double* y, double E);

    

    void get_null_vector_in_direction(double e1[],double e2[],double e3[],double P_obs_loc[],double V_obs_loc[],double V[],double V_loc[],double &psi);
	void local_to_cartesian(double o[],double e1[],double e2[], double e3[], double p_local[],double p_cart[],double psi);
	double error_norm(double *p,double *x);
	double min_step_size(double* y);
	int inside_universe(double *p);
	
	double dust_time(double r);
	bool inside_dust_cloud(double r,double t);
	bool exiting_dust_cloud(double X,double Xdot);

	double findDustBoundaryCrossingBisection(double y[4],double yp[3],double y_new[4],double yp_new[3],double tol,int maxIters);
	double findDustBoundaryCrossing(double y[4],double y_new[4],double tol,int maxIters);
	double findDustExitCrossingBisection(double y[4],double yp[3],double y_new[4],double yp_new[3],double tol,int maxIters);
	double findTauZeroCrossing(double y[4],double y_new[4],double tol,int maxIters);
	double findEtaZeroCrossing(double y[4],double y_new[4],double tol,int maxIters);
	void to_interior_coordinates(double y[],double yp[],double y_new[],double yp_new[]);
	void to_interior_coordinates_eta(double y[],double yp[],double y_new[],double yp_new[]);
	void to_exterior_coordinates(double y[],double yp[],double y_new[],double yp_new[]);
	void to_exterior_coordinates_eta(double y[],double yp[],double y_new[],double yp_new[]);
	void to_exterior_coordinates_through_tau_equals_zero(double y[],double yp[],double y_new[],double yp_new[]);
	void to_exterior_coordinates_through_tau_equals_zero_eta(double y[],double yp[],double y_new[],double yp_new[]);
	
	double get_redshift(double r);
	
	double M;
private:
    /// Example of storing initial position
    double pos[3];

    /// Gravitational mass of the dust (exterior Schwarzschild)

    /// Initial star radius at t=0 (matching or user-chosen)
    double R0;
	double X0;
	double am;
	
	double rho0; // initial density.
	
	double C;
	
	int t_dir = -1;
	double h;
	bool inside=false;
	bool inside_prev=false; // we need to know not only about the current step, but the previous one.
	
	bool use_eta = false;

	vector<double> m_tau;
	vector<double> m_a;

};

#endif // OPPENHEIMER_SNYDER_COLLAPSE_H