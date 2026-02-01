#ifndef MULTI_WORMHOLE_H
#define MULTI_WORMHOLE_H

#include <vector>
#include <string>
#include <array>
#include "Spacetime.h"

using namespace std;


class Multi_worm_hole : public Spacetime{
public:
	Multi_worm_hole(vector<array<double, 3>> positions_,double tol_, vector<double> as_,vector<int> stitching_);
	bool accel_1st_order(double*& yp,double* y,double E); // acceleration as a 1st order system - returns false if we attempt to compute at a point outside of the universe
	double error_norm(double *p,double *x);  
	int inside_universe(double *p);
	double get_h0(double *p);
	double compute_energy(double *p);
	void scale_velocity(double v[],double *p);
	void build_plane_coordinates(double P[],double V[]);
	void teleport(int i,int j,double *p);
	double min_step_size(double* y);
private:
	double scale_a = 1.05; 
	vector<array<double, 3>> positions;
	vector<double> as;
	vector<int> stitching;
	int num_wormholes;
};

#endif