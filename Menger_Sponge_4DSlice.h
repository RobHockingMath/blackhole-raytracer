#ifndef MENGER_4DSLICE_H
#define MENGER_4DSLICE_H

#include "Distance_fractal.h"


class Menger_Sponge_4DSlice : public Distance_fractal{
public:
	//first 3 are x,y,z of center, next is the scale, then the depth,then type, 
	// then epsilon,then offset, then three ints for the color
	
	Menger_Sponge_4DSlice(double,double,double,double,int,int,double,double,int,int,int);
	double Distance(double[]);
	bool wireframe_condition(double p[]);
	double is_Outside(double p[]);
	void rot(double p[]);
	void get_Texture_Coords(double[],double[]);
	void get_Bump_Perturbation(double x[], vector<Texture*> textures, double dn[]);
	//void get_Normal(double[],double[]);

private:
	int depth;
	double scale;
	double offset;
	int type;
};

#endif
