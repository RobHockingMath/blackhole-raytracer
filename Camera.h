#ifndef CAMERA_H
#define CAMERA_H

#include <vector>

using namespace std;

class Camera {
public:
	Camera(double[],double[],double[],double);
	void reset_position(double P[], double L[], double U[]);
	virtual void get_rays(int,int,const int, const int,vector< vector<double> >*,vector< vector<double> >*){return;};
	virtual void get_rays_geola(int,int,const int, const int,vector< vector<double> >*,vector< vector<double> >*){return;};
	virtual void get_rays_chimera(int,int,const int, const int,vector< vector<double> >*,vector< vector<double> >*){return;};
	void get_basis(double o[],double e1[],double e2[],double e3[]);
	double p[3]; //position
	double l[3]; //what it's looking at
	double up[3]; //this direction will point up in the pic
	double O[3][3]; //orthonormal frame

	double field_of_view;
	//double h_m;
	//double h_M;
	//double h_c;

};

class Pinhole_Camera : public Camera{
public:
	Pinhole_Camera(double[],double[],double[],double);
	void get_rays(int,int,const int, const int,vector< vector<double> >*,vector< vector<double> >*);
	void get_rays_geola(int i, int j, const int w, const int h, vector< vector<double> >* P,vector< vector<double> >* V);
	void get_rays_chimera(int i, int j, const int w, const int h, vector< vector<double> >* P,vector< vector<double> >* V);
};

class Lense_Camera : public Camera{
public: 
	Lense_Camera(double[],double[],double[],double,double,double,double,int);
	void get_rays(int,int,const int,const int,vector< vector<double> >*,vector< vector<double> >*);
private:
	int num_rays;
	double len_W;
	double len_H;
	double focal_length;
};

#endif

