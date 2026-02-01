#ifndef MENGER_H
#define MENGER_H

#include "Distance_fractal.h"


class Menger_Sponge : public Distance_fractal{
public:
	//first 3 are x,y,z of center, next is the scale, then the depth, 
	// then epsilon, then three ints for the color
	
	Menger_Sponge(double,double,double,double,int,double,double,int,int,int);
	double Distance(double[]);
	bool wireframe_condition(double p[]);
	bool rat;

private:
	int depth;
	double shift;
	double scale;
};

class Upright_Cylinder_Infinite : public Distance_fractal{
public:
	Upright_Cylinder_Infinite(double zt, double zb,double x0, double y0, double r,double E, int red,int green, int blue);
	void inverse_texture(int i,int j,int n, int m,float query_pt[]);
	double Distance(double[]);
	bool wireframe_condition(double p[]);
	void get_Texture_Coords(double x[], double u[]);
	void get_Shading(double x[],vector<Texture*> textures,int obj_id,int C[]);
	bool has_upper_bound,has_lower_bound;
	double upper_bound,lower_bound;
private:
	double x_0,y_0,R,z_b,z_t;
};

#endif
