#ifndef MANDELBULB_H
#define MANDELBULB_H

#include "Distance_fractal.h"


class Mandelbulb : public Distance_fractal{
public:
	//first 3 are x,y,z of center, then epsilon, then power (i.e. is it a z^2, z^8 mandelbulb or what?),
	//then r g b values of color
	
	
	Mandelbulb(double,double,double,double,int,int,int,int);
	void Mult(double[],double[]);
	void raise_to_power(double[],int,double[]);
	void raise_to_power_overwrite(double[],int);
	void Next(double[],double[],double[]);
	double Distance(double[]);
	bool wireframe_condition(double[]);


private:
	//the C in Qn+1=Qn^2+C is split into cS and cV,
	//the scalar and vector parts respectively
	int power;
};

#endif
