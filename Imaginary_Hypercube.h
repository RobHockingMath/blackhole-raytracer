#ifndef IMAGINARY_HYPERCUBE
#define IMAGINARY_HYPERCUBE

#include <vector>
#include "Sphere.h"
#include "Shape.h"

class Imaginary_Hypercube : public BasicShape{
public:
	//first 3 are x,y,z of center, then epsilon, then power (i.e. is it a z^2, z^8 mandelbulb or what?),
	//then r g b values of color
	
	
	Imaginary_Hypercube(double N[],bool is_T_,int red, int green, int blue);
	double get_Intersection(double p[],double v[],double x[]);
	double get_Intersection_depth(double p[],double v[],double x[],int depth);
	void get_Texture_Coords(double x[], double u[]);
	void get_Normal(double x[], double n[]);
	void set_Texture_IDs(int IDs[]);
	void get_Color(double x[], vector<Texture*> textures, int col[]);
	bool wireframe_condition(double[]);

private:
	double h;
	int recursion_depth = 4;
	double s = 3;
	double Vertices_H[14][4];
	double Vertices_T[17][4];
	
	std::vector<double> X_H3, Y_H3, Z_H3;
	std::vector<double> X_T3, Y_T3, Z_T3;
	
	double Translations_H[27][4];
	double Translations_T[27][4];
	double Translations_H3[27][3];
	double Translations_T3[27][3];
	
	//normal vector
	double n[4];
	double e1[4],e2[4],e3[4];
	bool is_T;
	vector<Triangle*> Triangles;
	
};

#endif