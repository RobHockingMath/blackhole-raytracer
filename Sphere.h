#ifndef SPHERE_H
#define SPHERE_H

#include "Shape.h"


class Sphere : public BasicShape{
public:
	//first 3 are x,y,z of center, then radius, then 
	//r g b values of color
	Sphere(double,double,double,double,int,int,int);
	double get_Intersection(double[],double[],double[]);
	double get_Intersection2(double[],double[],double,double[]);
	void get_Texture_Coords(double[],double[]);
	void get_Normal(double[],double[]);
	bool wireframe_condition(double[]);

private:
	double r; //radius
	double c[3];//center
};

class PlusSign3D : public BasicShape{
public:
	//first 3 are x,y,z of center, then width of central cube, then 
	//r g b values of color
	PlusSign3D(double,double,double,double,int,int,int);
	double get_Intersection(double[],double[],double[]);
	double get_Intersection2(double[],double[],double,double[]);
	void get_Texture_Coords(double[],double[]);
	void get_Normal(double[],double[]);
	bool wireframe_condition(double[]);

private:
	double width; //radius
	double c[3];//center
	
};

class Parallelogram : public BasicShape{
public:
	Parallelogram(double[],double[],double[],int,int,int);
	double get_Intersection(double[],double[],double[]);
	double get_Intersection2(double[],double[],double,double[]);
	void cheap_hack(double[],double[]);
	void get_Texture_Coords(double[],double[]);
	void get_Normal(double[],double[]);
	bool wireframe_condition(double[]);
	int rat;
private:
	double c[3]; //center
	double o[3];//corner
	double x1[3],x2[3];//edges
	double n[3];//normal
};

class Triangle : public BasicShape{
public:
	Triangle(double p1_[],double p2_[],double p3_[], double n_[], int red,int green,int blue);
	double get_Intersection(double[],double[],double[]);
	double get_Intersection2(double[],double[],double,double[]);
	//void cheap_hack(double[],double[]);
	void get_Texture_Coords(double[],double[]);
	void get_Normal(double[],double[]);
	bool wireframe_condition(double[]);
	void scaled_translation(double x,double y, double z,double s);
	int rat;
	double p1[3],p2[3],p3[3]; //corners
	double n[3];//normal
};

class InfinitePlane: public BasicShape{
public:	
	InfinitePlane(double Z, int red,int green,int blue);
	double get_Intersection(double[],double[],double[]);
	double get_Intersection2(double[],double[],double,double[]);
	//void cheap_hack(double[],double[]);
	void get_Shading(double x[],vector<Texture*> textures,int obj_id,int C[]);
	void inverse_texture(int i,int j,int n, int m,float query_pt[]);
	void place_dot(double x[],int ***im_data,int w,int h);
	void get_Texture_Coords(double[],double[]);
	void get_Bump_Perturbation(double x[], vector<Texture*> textures, double dn[]);
	void get_Normal(double[],double[]);
	bool wireframe_condition(double[]);
	bool tiling;
	int coordinate;
	int normal_dir;
	double tile_size;
	int upper_bound_coordinate;
	int upper_bound_variable;
	double upper_bound;
	double upper_bound_coefficient;
	double texture_center[2];
private:	
	double z;
};


class Upright_Cylinder : public BasicShape{
public:
	Upright_Cylinder(double zb,double zt,double x0, double y0, double r,int red,int green, int blue);
	double get_Intersection(double[],double[],double[]);
	double get_Intersection2(double[],double[],double,double[]);
	//void cheap_hack(double[],double[]);
	void get_Texture_Coords(double[],double[]);
	void get_Normal(double[],double[]);
	bool wireframe_condition(double[]);

private:
	double z_b,z_t,x_0,y_0,R;
	
};


class Cylinder : public BasicShape{
public:
	Cylinder(double P1[3],double P2[3],double R,int red,int green, int blue);
	double get_Intersection(double[],double[],double[]);
	double get_Intersection2(double[],double[],double,double[]);
	//void cheap_hack(double[],double[]);
	void get_Texture_Coords(double[],double[]);
	void get_Normal(double[],double[]);
	bool wireframe_condition(double[]);
private:
	double p1[3];
	double p2[3];
	double r;
	double H;
	
	double e1[3];
	double e2[3];
	double e3[3];
};


class Cube : public CompoundShape{
public:
	//center, then dims, then color
	Cube(double,double,double,double,double,double,int,int,int);
	bool wireframe_condition(double[]);

private:
	double dim[3]; //dimensions
	double c[3]; //center
};

/**class ConvexHull : public CompoundShape{
public:
	ConvexHull(std::vector<double> X,std::vector<double> Y,std::vector<double> Z,int Red, int Green, int Blue);
	bool wireframe_condition(double[]);
private:
	int num_triangles;
};**/

class PlusSign3DMesh: public CompoundShape{
public:
	//first 3 are x,y,z of center, then width of central cube, then 
	//r g b values of color
	PlusSign3DMesh(double,double,double,double,int,int,int);
	bool wireframe_condition(double[]);

private:
	double width; //dimensions
	double c[3]; //center
};

#endif
