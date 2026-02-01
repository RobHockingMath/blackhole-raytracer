#ifndef MENGER_TYPE2_H
#define MENGER_TYPE2_H

#include "Distance_fractal.h"


class Menger_Sponge_Type2 : public Distance_fractal{
public:
	//first 3 are x,y,z of center, next is the scale, then the depth, 
	// then epsilon, then three ints for the color
	
	Menger_Sponge_Type2(double c1,double c2,double c3,double scale_,int depth_,double E,int red,int green, int blue);
	double Distance(double[]);

private:
	double scale;
	double depth;
};

#endif
