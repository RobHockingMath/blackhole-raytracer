#ifndef MINKOWSKI_H
#define MINKOWSKI_H

#include <vector>
#include <string>

#include "Spacetime.h"

using namespace std;


class Minkowski : public Spacetime{
public:
	Minkowski(double tol_);
	bool accel_1st_order(double*& yp,double* y,double E); // acceleration as a 1st order system - returns false if we attempt to compute at a point outside of the universe
	double error_norm(double *p,double *x);  
	int inside_universe(double *p);
	double get_h0(double *p);
	double compute_energy(double *p);
	void scale_velocity(double v[],double *p);
};

#endif