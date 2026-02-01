#ifndef SPACETIME_H
#define SPACETIME_H

#include <vector>
#include <string>
#include <math.h>
#include <iostream>

using namespace std;


class Spacetime{
public:
	Spacetime(double tol_);
	virtual bool accel_1st_order(double*& yp,double* y, double E){return false;}; // acceleration as a 1st order system - returns false if we attempt to compute at a point outside of the universe
	virtual double error_norm(double *p,double *x){return 0;};  
	virtual int inside_universe(double *p){return 0;}; 
	virtual double get_h0(double *p){return 1.0;}; 
	virtual double compute_energy(double *p){return 1;}; 
	virtual void scale_velocity(double v[],double *p){return;};
	virtual double get_a(double u,double v){return 0;}
	virtual double get_r(double u,double v){return 0;}
	virtual void teleport(int i,int j,double *p){return;};
	
	virtual void get_null_vector_in_direction(double e1[],double e2[],double e3[],double P_obs_loc[],double V_obs_loc[],double V[],double V_loc[],double &psi){return;};
	virtual void cartesian_to_local(double e1[],double e2[], double e3[], double p_cart[],double v_cart[],double p_local[],double v_local[],double &psi){p_local[0]=p_cart[0];p_local[1]=p_cart[1];p_local[2]=p_cart[2];v_local[0]=v_cart[0];v_local[1]=v_cart[1];v_local[2]=v_cart[2];}; // default is identity transform.
	virtual void local_to_cartesian(double o[],double e1[],double e2[], double e3[], double p_local[],double p_cart[],double psi){p_cart[0]=p_local[0];p_cart[1]=p_local[1];p_cart[2]=p_local[2];};  // default is identity transform.
	bool take_step(double P[],double V[],double P2[],double V2[],double &E,double &H,bool &step_taken, bool force_step); // assumed to be in local manifold coordinates.
	bool Adaptive_Step(double*& y,bool& y_changed,int n,double E,double& h,double tol, bool force_step);
	virtual double min_step_size(double* y){return 1.0;};
	int encode(int n, int m) { return (n << 16) | m;};
	void decode(int N, int& n, int& m) {n = N >> 16;m = N & 0xFFFF;};
	void recompute_orthonormal_frame(double o[],double p[],double u[],double e1_new[],double e2_new[],double e3_new[]);
	void recompute_orthonormal_frame2(double o[],double p[],double e1_new[],double e2_new[],double e3_new[]);
	void orthonormalBasis(double e1[3], double e2[3], double e3[3]);
	
	bool is_wormhole;
	bool is_collapse;
	bool has_dustcloud;
	bool inside=false; // only relevant for spacetimes with a dust cloud.
	bool inside_prev=false; // only relevant for spacetimes with a dust cloud.
	virtual void to_interior_coordinates(double P[],double V[],double P_interior[],double V_interior[]){P_interior[0]=P[0];P_interior[1]=P[1];P_interior[2]=P[2];P_interior[3]=P[3];V_interior[0]=V[0];V_interior[1]=V[1];V_interior[2]=V[2];};
	virtual void to_interior_coordinates_eta(double P[],double V[],double P_interior[],double V_interior[]){P_interior[0]=P[0];P_interior[1]=P[1];P_interior[2]=P[2];P_interior[3]=P[3];V_interior[0]=V[0];V_interior[1]=V[1];V_interior[2]=V[2];};
	virtual void to_exterior_coordinates(double P[],double V[],double P_exterior[],double V_exterior[]){P_exterior[0]=P[0];P_exterior[1]=P[1];P_exterior[2]=P[2];P_exterior[3]=P[3];V_exterior[0]=V[0];V_exterior[1]=V[1];V_exterior[2]=V[2];};
	virtual void to_exterior_coordinates_eta(double P[],double V[],double P_exterior[],double V_exterior[]){P_exterior[0]=P[0];P_exterior[1]=P[1];P_exterior[2]=P[2];P_exterior[3]=P[3];V_exterior[0]=V[0];V_exterior[1]=V[1];V_exterior[2]=V[2];};
	virtual void to_exterior_coordinates_through_tau_equals_zero(double P[],double V[],double P_exterior[],double V_exterior[]){P_exterior[0]=P[0];P_exterior[1]=P[1];P_exterior[2]=P[2];P_exterior[3]=P[3];V_exterior[0]=V[0];V_exterior[1]=V[1];V_exterior[2]=V[2];};
	virtual void to_exterior_coordinates_through_tau_equals_zero_eta(double P[],double V[],double P_exterior[],double V_exterior[]){P_exterior[0]=P[0];P_exterior[1]=P[1];P_exterior[2]=P[2];P_exterior[3]=P[3];V_exterior[0]=V[0];V_exterior[1]=V[1];V_exterior[2]=V[2];};
	virtual bool inside_dust_cloud(double r,double t){return false;};
	virtual bool exiting_dust_cloud(double X,double Xdot){return false;};
	virtual double findDustBoundaryCrossing(double y[4],double y_new[4],double tol,int maxIters){return -1;};
	virtual double findDustBoundaryCrossingBisection(double y[4],double yp[3],double y_new[4],double yp_new[3],double tol,int maxIters){return -1;};
	virtual double findDustExitCrossingBisection(double y[4],double yp[3],double y_new[4],double yp_new[3],double tol,int maxIters){return -1;};
	virtual double findTauZeroCrossing(double y[4],double y_new[4],double tol,int maxIters){return -1;};
	virtual double findEtaZeroCrossing(double y[4],double y_new[4],double tol,int maxIters){return -1;};
	virtual double dust_time(double r){return 0;};
	
	virtual double get_redshift(double r){return 1.0;};
	virtual double get_a(double tau){return 1.0;};
	double M=0.075;
protected:
	double e1[3];
	double e2[3];
	double radius;
	double phi;
	bool has_plane_coordinates;
private:
	double tolerance;
};

#endif