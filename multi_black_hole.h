#ifndef MULTI_BLACKHOLE_H
#define MULTI_BLACKHOLE_H

#include <vector>
#include <string>
#include <array>
#include "Spacetime.h"

using namespace std;


class Multi_black_hole : public Spacetime{
public:
	Multi_black_hole(vector<array<double, 3>> positions_,double tol_, vector<double> masses_,vector<int> stitching_);
	bool accel_1st_order(double*& yp,double* y,double E); // acceleration as a 1st order system - returns false if we attempt to compute at a point outside of the universe
	double error_norm(double *p,double *x);  
	int inside_universe(double *p);
	double get_h0(double *p);
	double compute_energy(double *p);
	void scale_velocity(double v[],double *p);
	void build_plane_coordinates(double P[],double V[]);
	void teleport(int i,int j,double *p);
	double min_step_size(double* y);
	void rotate(double x[], double v[], double theta);
private:
	vector<array<double, 3>> positions;
	vector<double> masses;
	vector<int> stitching;
	int num_black_holes;
};

#endif