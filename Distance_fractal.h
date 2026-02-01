#ifndef FRACTAL_H
#define FRACTAL_H

#include <vector>
#include "Shape.h"
#include "Texture.h"
#include "Sphere.h"

using namespace std;

class Distance_fractal : public BasicShape{
public:
	
	Distance_fractal(double,double,double,double,int,int,int);
	virtual double Distance(double[]){return 0;};
	double get_Intersection(double[],double[],double[]);
	double get_Intersection2(double[],double[],double,double[]);
	//void get_Texture_Coords(double[],double[]);
	//void get_Bump_Perturbation(double x[], vector<Texture*> textures, double dn[]);
	void get_Texture_Basis(double[],double[],double[]);
	void get_Normal(double[],double[]);
	virtual bool wireframe_condition(double[]){return false;};
	double c[3];
	double epsilon;
private:
	
protected:
	double bounding_radius;
};

#endif
