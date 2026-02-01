#ifndef QUATERNION_H
#define QUATERNION_H

#include "Distance_fractal.h"


class Quaternion_Julia_Set : public Distance_fractal{
public:
	//first 3 are x,y,z of center, next 3 are R,I,and theta, then epsilon, 
	//then r g b values of color
	
	//even though defining a set by 4 numbers c1,c2,c3,c4 is more natural, 
	//it seems to be the convention to use 3 parameters R,I, and theta
	//so in Feb 2010 I changed the code to fit the convention.  This way
	//one can look up sets online, read of the R,I,theta values, and reproduce
	//the same sets with this code.
	//Why 3 parameters?  Because of the symmetry relation J_Rq=RJ_q, the space of possible
	//quat julia sets is really 2d.  The third parameter is somehow related to the angle
	//of the 3d hyperplane you intersect the set with.
	
	Quaternion_Julia_Set(double,double,double,double,double,double,double,int,int,int);
	void Mult(double[],double[]);
	void Next(double[],double[]);
	double Distance(double[]);
	bool wireframe_condition(double[]);

private:
	//the C in Qn+1=Qn^2+C is split into cS and cV,
	//the scalar and vector parts respectively
	double C[4];
	double e2[4];
	double Theta;
};

#endif
