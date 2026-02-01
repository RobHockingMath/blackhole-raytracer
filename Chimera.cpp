#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>

using namespace std;
#include "Chimera.h"
#define PI 3.14159

ChimeraCameraArray::ChimeraCameraArray(double R, double H_m, double H_M, double UP[3],double arclength, double holo_W, double holo_H, double hogel_size, double zoom){

	// all cameras in the array look at the same point and have the same up vector.
	l[0]=0;l[1]=0;l[2]=0;
	up[0]=UP[0];up[1]=UP[1];up[2]=UP[2];
	//we align our cylinder relative to the up vector.  This will break if the up vector is parallel to e1.
	double len = sqrt(up[0]*up[0]+up[1]*up[1]+up[2]*up[2]);
	up[0]/=len;up[1]/=len;up[2]/=len;
	
	double e1[3]={1-(1*up[0])*1,0,0};
	len = sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);
	e1[0]/=len;e1[1]/=len;e1[2]/=len;
	double e2[3]={up[1]*e1[2]-e1[1]*up[2],e1[0]*up[2]-up[0]*e1[2],up[0]*e1[1]-e1[0]*up[1]}; 
	
	cout<<e1<<endl;
	cout<<e2<<endl;
	cout<<up<<endl;
	
	// zoom is a percentage, so that a zoom of 100 does nothing.
	
	double F=2*tan(holo_W/(2*R))/(zoom/100);
	// let's build the array of positions and location frames for each camera in the area.
	vector<Pinhole_Camera*> column;
	for(int i = 0 ; i<768;i++){
		cout<<"building Chimera array...column "<<i<<" of "<<768<<endl;
		double theta = arclength/2-i*arclength/767;
		column.clear();
		for(int j=45;j<90;j++){
			double z = H_m + j*(H_M-H_m)/89;
			p[0] = R*cos(theta)*e1[0]+R*sin(theta)*e2[0]+z*up[0];
			p[1] = R*cos(theta)*e1[1]+R*sin(theta)*e2[1]+z*up[1];
			p[2] = R*cos(theta)*e1[2]+R*sin(theta)*e2[2]+z*up[2];
			
			Pinhole_Camera* C = new Pinhole_Camera(p,l,up,F);
			column.push_back(C);
			
		}
		cameras.push_back(column);
	}
	holoW=holo_W;
	holoH=holo_H;
	hogelSize = hogel_size;

	// now let's compute the width and height of the images in pixels.
	// notice the coversion factor of 10000 because holoW is in cm but
	// hogelSize is in micrometers.
	
	// also notice the scale factor of 1.1 as 10% of each image is cut away.
	
	w = floor(1.1*holoW*10000/(hogel_size));
	h = floor(1.1*holoH*10000/(hogel_size));

}








/**Pinhole_Camera::Pinhole_Camera(double P[], double L[], double U[], double F)
	: Camera(P,L,U,F)
{

}

//returns a line passing from pixel i,j on the back of the camera
//through the pinhole.
//Line is returned by returning two vectors P and V, intended
//to be the vectors of the vector form of a line L(t)=P+t*V.
//P and V are taken as parameters and modified.

void Pinhole_Camera::get_rays(int i, int j, const int w, const int h, vector< vector<double> >* P,vector< vector<double> >* V){

	int num_rays = 3;

	srand( time(NULL) );
	
	double i_r,j_r; //these will hold small random perturbations in position index

	P->clear();V->clear();
	vector<double> P1,V1;
	//in this case P and V are 1 by 3

	for(int ii=0;ii<num_rays;ii++){
		for(int jj=0;jj<num_rays;jj++){
			//P is easy, it's just the camera location
			for(int k=0;k<3;k++){P1.push_back(p[k]);}
			P->push_back(P1);

			i_r=((rand() % 100)-50)/100.0; //i.e. rand(-0.5,0.5)
			j_r=((rand() % 100)-50)/100.0; //i.e. rand(-0.5,0.5)

			//cout<<"float(ii+i_r)/num_rays : "<<float(ii+i_r)/num_rays<<" float(jj+j_r)/num_rays : "<<float(jj+j_r)/num_rays<<endl;
			//V is more complicated
			double x=field_of_view*((float(i)+float(ii+i_r)/num_rays)/w-0.5);
			double y=field_of_view*((float(j)+float(jj+j_r)/num_rays)/h-0.5);
			//double x=field_of_view*((float(i))/w-0.5);
			//double y=field_of_view*((float(j))/h-0.5);
			//cout<<"x y = "<<x<<" "<<y<<endl;
			double theta=sqrt(x*x+y*y);
			double phi=atan2(y,x);
			for(int k=0;k<3;k++){
				V1.push_back(cos(theta)*O[0][k]+sin(theta)*(cos(phi)*O[1][k]+sin(phi)*O[2][k]));
			}
			V->push_back(V1);
			P1.clear();
			V1.clear();
		}
	}
	
}

Lense_Camera::Lense_Camera(double P[], double L[], double U[], double F,double W,double H,double f,int num)
	: Camera(P,L,U,F)
{
	num_rays=num;
	len_W=W;
	len_H=H;
	focal_length=f;
}

void Lense_Camera::get_rays(int i, int j, const int w, const int h,  vector< vector<double> >* P,vector< vector<double> >* V){

	vector<double> P1,V1;
	double x,y,z,xp,yp,theta,phi,l;
	
	P->clear();
	V->clear();
	
	//we'll need random numbers so things look good
	srand( time(NULL) );
	
	double i_r,j_r; //these will hold small random perturbations in position index
	
	//iterate over positions on the lense
	int count=0;
	for(int ii=0;ii<num_rays;ii++){
		for(int jj=0;jj<num_rays;jj++){
			//these are the xyz coords of point (i,j) on the lense
			
			i_r=((rand() % 100)-50)/100.0; //i.e. rand(-0.5,0.5)
			j_r=((rand() % 100)-50)/100.0; //i.e. rand(-0.5,0.5)
			
			for(int k=0;k<3;k++){
				P1.push_back(p[k]+((ii+i_r)*len_W/num_rays-0.5*len_W)*O[1][k]+((jj+j_r)*len_H/num_rays-0.5*len_H)*O[2][k]);
			}
			//add position to our list of positions
			P->push_back(P1);
			
			//now we figure out the corresponding position on the image plane
			x=2*focal_length*tan(0.5*field_of_view)*(float(i)/w-0.5);
			y=2*focal_length*tan(0.5*field_of_view)*(float(j)/h-0.5);
			//and get the displacement vector
			xp=len_W*(float(ii+i_r)/num_rays-0.5);
			yp=len_H*(float(jj+j_r)/num_rays-0.5);
			x=x-xp;
			y=y-yp;
			z=focal_length;
			l=sqrt(x*x+y*y+z*z);
			x=x/l;y=y/l;z=z/l;
			
			for(int k=0;k<3;k++){
				V1.push_back(z*O[0][k]+x*O[1][k]+y*O[2][k]);
			}
			V->push_back(V1);
			count++;
			P1.clear();
			V1.clear();
		}
	}

}**/
