#ifndef MENGER_4DPROJECTION_H
#define MENGER_4DPROJECTION_H

#include "Distance_fractal.h"


class Menger_Sponge_4DProjection : public Distance_fractal{
public:
	//first 3 are x,y,z of center, next is the scale, then the depth,then type, 
	// then epsilon,then offset, then three ints for the color
	
	Menger_Sponge_4DProjection(double,double,double,double,int,double,double,double,double,double,double,int,int,int);
	double Distance(double[]);
	bool wireframe_condition(double p[]);
	//double is_Outside(double p[]);
	//void rot(double p[]);
	//void get_Normal(double[],double[]);

private:
	int depth;
	double scale;
	double offset;
	int type;
	
	double e1[4],e2[4],e3[4],N[4];
	
};

#endif
