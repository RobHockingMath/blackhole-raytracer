#ifndef JULIABULB_H
#define JULIABULB_H

#include "Distance_fractal.h"


class Juliabulb : public Distance_fractal{
public:
	//first 3 are x,y,z of center, then epsilon, then power (i.e. is it a z^2, z^8 juliabulb or what?),
	//then c1 c2 c3, then r g b values of color
	
	
	Juliabulb(double,double,double,double,int,double,double,double,int,int,int);
	double Distance(double[]);
	bool wireframe_condition(double[]);

private:
	//the C in Qn+1=Qn^2+C is split into cS and cV,
	//the scalar and vector parts respectively
	int power;
	double c[3];
};

#endif
