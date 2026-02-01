#ifndef SHAPE_H
#define SHAPE_H

#include <vector>
#include "Texture.h"
#include "invisibility.h"
#include <iostream>

using namespace std;

class Shape {
public: 
	//r g b values of color
	Shape(int,int,int);

	//these 2 compute intersects for normal light rays
	//and and light rays bending around the invisibility 
	//cloak respectively
	virtual double get_Intersection(double[],double[],double[]){return 0;};
	virtual double get_Intersection2(double[],double[],double,double[]){return 0;};

	virtual void get_Texture_Coords(double[],double[]){return;};
	virtual void inverse_texture(int i,int j,int n, int m,float query_pt[]){return;};
	virtual void place_dot(double x[],int ***im_data,int w,int h){return;};
	virtual void get_Normal(double[],double[]){return;};
	virtual void get_Color(double[],vector<Texture*>,int[]){return;};
	virtual void get_Shading(double y[],vector<Texture*>,int obj_id,int C[]){return;};
	virtual double get_Alpha(double x[], vector<Texture*> textures){return 0;};
	virtual void get_Bump_Perturbation(double x[], vector<Texture*> textures, double dn[]){return;};
	void get_reflection(double[],double[]);
	bool get_refracted_ray(double[],double[]);
	double reflectivity;
	double transparency;
	double ref_index; //index of refraction
	int texture_ID;
	int alpha_ID;
	int bump_ID;
	virtual bool wireframe_condition(double[]){return false;};
	virtual void scaled_translation(double x,double y, double z,double s){return;};
	bool is_wireframe;
	bool self_illuminating;
	bool reflection_only;
	bool one_sided;
protected:
	int color[3];
};

class BasicShape : public Shape{
public:
	BasicShape(int,int,int);
	virtual double get_Intersection(double[],double[],double[]){return 0;};
	virtual double get_Intersection2(double[],double[],double,double[]){return 0;};
	virtual void get_Texture_Coords(double[],double[]){return;};
	virtual void get_Normal(double[],double[]){return;};
	virtual bool wireframe_condition(double[]){return false;};
	void get_Color(double[],vector<Texture*>,int[]);
	double get_Alpha(double x[], vector<Texture*> textures);
	//void get_Bump_Perturbation(double x[], vector<Texture*> textures, double dn[]);
	//double get_Bump_Value(double u[], vector<Texture*> textures);
};

class CompoundShape : public Shape{
public:
	CompoundShape(int,int,int);
	double get_Intersection(double[],double[],double[]);
	double get_Intersection2(double[],double[],double,double[]);
	void get_Texture_Coords(double[],double[]);
	void get_Normal(double[],double[]);
	virtual bool wireframe_condition(double[]){return false;};
	void get_Color(double[],vector<Texture*>,int[]);
	double get_Alpha(double x[], vector<Texture*> textures);
	void set_Texture_IDs(int[]);
	vector<BasicShape*> parts;
};

#endif
