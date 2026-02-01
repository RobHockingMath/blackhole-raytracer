#ifndef MENGER_4DSLICE_TYPE2_H
#define MENGER_4DSLICE_TYPE2_H

#include "Distance_fractal.h"


class Menger_Sponge_4DSlice_Type2 : public Distance_fractal{
public:
	//first 3 are x,y,z of center, next is the scale, then the depth,then type, 
	// then epsilon,then offset, then three ints for the color
	
	Menger_Sponge_4DSlice_Type2(double c1,double c2,double c3,double scale_,int depth_,int type_,double E,double offset_,int red,int green, int blue);
	double Distance(double[]);
	double Distance3D(double x,double y, double z);
	double flip = false;

private:
	int depth;
	double scale;
	double offset;
	int type;
};

#endif
