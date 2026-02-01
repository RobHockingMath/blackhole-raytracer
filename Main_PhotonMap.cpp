#include <iostream>
#include <cstdlib>
#include <math.h>
#include "World.h"
#include "Shape.h"
#include "Sphere.h"
#include "Camera.h"
#include "Texture.h"
#include "Quaternion_Julia_Set.h"
#include "Menger_Sponge_4DSlice.h"
#include "Menger_Sponge_4DSlice_Type2.h"
#include "Menger_Sponge.h"
#include "Menger_Sponge_Type2.h"
#include "Mandelbulb.h"
#include "Minkowski.h"
#include "Schwarzschild.h"
#include <string>
#include <cstring>
#include <omp.h>
#include "PhotonMapLight.h"

#define pi 3.1415
#define SQRT3 1.7321

using namespace std;


int main(int argc,char* argv[])
{
//default values for output image filename 
//and width and height, if not entered by command line
char* filename="world.bmp";
int w=1920;
int h=1080;
int eye=1;
double offset = 0;
//attempt to read in these things
if(argc>1){filename=argv[1];}
if(argc>2){w=atoi(argv[2]);}
if(argc>3){h=atoi(argv[3]);}
//if(argc>4){eye=atoi(argv[4]);}
if(argc>4){offset=stod(argv[4]);}

int rat = 30;
string base = "rat";
string extension = to_string(rat);
cout <<(base+extension).c_str()<<endl;



double c1=0,c2=0*0.5,c3=0,scale = 1,E=0.001;
int depth = 4;
E=0.0001;
//double offset = 1.0;
	
	

char* water_bump = "./Water2/bump0.txt";
Texture* T1 = new Texture(water_bump);

			
			
World earth(0);


InfinitePlane* P0 = new InfinitePlane(-6,0,255,255);
P0->normal_dir = -1;
P0->texture_ID=4;
//P0->self_illuminating=true;
P0->tile_size=12;
//P0->bump_ID=2;

InfinitePlane* P1 = new InfinitePlane(2,0,255,255);
P1->texture_ID=2;
P1->bump_ID=0;
P1->self_illuminating=true;
P1->reflectivity=0.5;
P1->tile_size=12;
//Parallelogram *P1=new Parall

//InfinitePlane* P2 = new InfinitePlane(0,255,0,0);
//P2->texture_ID=4;
//P2->self_illuminating=true;
//P2->alpha_ID=4;
//P2->tiling=false;
//P1->transparency=0.5;

InfinitePlane* P3 = new InfinitePlane(6,0,255,255);
P3->coordinate = 0;
P3->texture_ID=0;
P3->self_illuminating=true;
P3->tile_size=12;
InfinitePlane* P4 = new InfinitePlane(-6,0,255,255);
P4->normal_dir = -1;
P4->coordinate = 0;
P4->texture_ID=0;
P4->tile_size=12;
//P4->self_illuminating=true;
InfinitePlane* P5 = new InfinitePlane(6,0,255,255);
P5->coordinate = 1;
P5->texture_ID=0;
P5->tile_size=12;
//P5->self_illuminating=true;
InfinitePlane* P6 = new InfinitePlane(-6,0,255,255);
P6->normal_dir=-1;
P6->coordinate = 1;
P6->texture_ID=0;
P6->tile_size=12;
//P6->self_illuminating=true;

double pos[3]={0,0,-2.5};
int color[3]={255,255,255};
double bb = 1000;
double radius = 0.75;

PhotonMapLight *L = new spherePhotonMapLight(pos, bb,radius, color);

double z_t=2;
double z_b=-6;
double x00 = 0.0;
double y00 = 2.3;
double rr = 0.5;
Upright_Cylinder_Infinite* C1 = new Upright_Cylinder_Infinite(z_t,z_b,x00,y00,rr,E,200,200,200);
//C1->self_illuminating=true;
C1->texture_ID=1;
Upright_Cylinder_Infinite* C2 = new Upright_Cylinder_Infinite(z_t,z_b,x00,-y00,rr,E,200,200,200);
//C2->self_illuminating=true;
C2->texture_ID=1;

Upright_Cylinder_Infinite* C3 = new Upright_Cylinder_Infinite(z_t,z_b,-3,y00,rr,E,200,200,200);
//C3->self_illuminating=true;
C3->texture_ID=1;
Upright_Cylinder_Infinite* C4 = new Upright_Cylinder_Infinite(z_t,z_b,-3,-y00,rr,E,200,200,200);
//C4->self_illuminating=true;
C4->texture_ID=1;


Minkowski* FlatSpace = new Minkowski(0.01);
double error_tol=1e-6;
double black_hole_radius=0.20;
int dimension = 4;
Schwarzschild* BlackHole = new Schwarzschild(error_tol,black_hole_radius,dimension);
earth.space_time=BlackHole;
//earth.space_time=FlatSpace;


earth.env_map_ID=3; //the texture ID of the background image - explained later

earth.w=w; //this is sloppy, but I have the world know how big the image is going to be
	  //so it can resize the background image accordingly.    
earth.h=h;

//earth.addShape(M);
//earth.addShape(M2);
//earth.addShape(M2p5);
earth.addTexture(T1);
//earth.addShape(M3);
earth.addShape(P1);
earth.addShape(P0);
earth.addShape(P3);
earth.addShape(P4);
earth.addShape(P5);
earth.addShape(P6);
//earth.addShape(P2);
earth.addShape(C1);
earth.addShape(C2);
earth.addShape(C3);
earth.addShape(C4);
earth.addPhotonMapLight(L);

earth.make_photon_map(1000,1000);
		
}
	