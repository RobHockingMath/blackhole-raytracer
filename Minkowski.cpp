#include "Minkowski.h"

Minkowski::Minkowski(double tol_)
	: Spacetime(tol_)
{
}

int Minkowski::inside_universe(double *p){
	return 1;
}

double Minkowski::get_h0(double *p){
	return 1; // take a big step.
}

// We think of p as being a point and x a vector.  
// This format is useful as it can be used directly in adaptive Runge Kutta routines.

double Minkowski::error_norm(double *p,double *x){
	
	return fabs(x[0])+fabs(x[1])+fabs(x[2]);
	
}

bool Minkowski::accel_1st_order(double*& yp,double* y, double E){
	
	// y and yp are giant vectors of the form y=[p,v] where p is position and v is velocity.
	
	// position derivative is velocity.
	yp[0]=y[3]; 
	yp[1]=y[4]; 
	yp[2]=y[5]; 
	
	// acceleration is zero.
	yp[3]=0;
	yp[4]=0;
	yp[5]=0;
	return true;
	
}

double Minkowski::compute_energy(double *p){
	return 0;
}

void Minkowski::scale_velocity(double v[],double *p){
	return;
}




