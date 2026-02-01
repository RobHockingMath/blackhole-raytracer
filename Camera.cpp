#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>

using namespace std;
#include "Camera.h"
#define PI 3.14159

Camera::Camera(double P[], double L[], double U[], double F){

//copy data into private variables
p[0]=P[0];p[1]=P[1];p[2]=P[2];
l[0]=L[0];l[1]=L[1];l[2]=L[2];
up[0]=U[0];up[1]=U[1];up[2]=U[2];
field_of_view=F*PI/180;

//make orthonormal frame
O[0][0]=L[0]-P[0];O[0][1]=L[1]-P[1];O[0][2]=L[2]-P[2];
double len=sqrt(O[0][0]*O[0][0]+O[0][1]*O[0][1]+O[0][2]*O[0][2]);
O[0][0]=O[0][0]/len;O[0][1]=O[0][1]/len;O[0][2]=O[0][2]/len;

//cout<<"len = "<<len<<endl;

O[1][0]=O[0][1]*U[2]-O[0][2]*U[1];
O[1][1]=O[0][2]*U[0]-O[0][0]*U[2];
O[1][2]=O[0][0]*U[1]-O[0][1]*U[0];

len=sqrt(O[1][0]*O[1][0]+O[1][1]*O[1][1]+O[1][2]*O[1][2]);
O[1][0]=O[1][0]/len;O[1][1]=O[1][1]/len;O[1][2]=O[1][2]/len;

O[2][0]=O[0][1]*O[1][2]-O[0][2]*O[1][1];
O[2][1]=O[0][2]*O[1][0]-O[0][0]*O[1][2];
O[2][2]=O[0][0]*O[1][1]-O[0][1]*O[1][0];

cout << O[0][0]<<" "<<O[0][1]<<" "<<O[0][2]<<endl;
cout << O[1][0]<<" "<<O[1][1]<<" "<<O[1][2]<<endl;
cout << O[2][0]<<" "<<O[2][1]<<" "<<O[2][2]<<endl;

cout << O[0][0]*O[1][0]+O[0][1]*O[1][1]+O[0][2]*O[1][2] <<endl;
cout << O[0][0]*O[2][0]+O[0][1]*O[2][1]+O[0][2]*O[2][2] <<endl;
cout << O[2][0]*O[1][0]+O[2][1]*O[1][1]+O[2][2]*O[1][2] <<endl;

cout << O[0][0]*O[0][0]+O[0][1]*O[0][1]+O[0][2]*O[0][2] <<endl;
cout << O[1][0]*O[1][0]+O[1][1]*O[1][1]+O[1][2]*O[1][2] <<endl;
cout << O[2][0]*O[2][0]+O[2][1]*O[2][1]+O[2][2]*O[2][2] <<endl;

}


void Camera::reset_position(double P[], double L[], double U[]){

//copy data into private variables
p[0]=P[0];p[1]=P[1];p[2]=P[2];
l[0]=L[0];l[1]=L[1];l[2]=L[2];
up[0]=U[0];up[1]=U[1];up[2]=U[2];

//make orthonormal frame
O[0][0]=L[0]-P[0];O[0][1]=L[1]-P[1];O[0][2]=L[2]-P[2];
double len=sqrt(O[0][0]*O[0][0]+O[0][1]*O[0][1]+O[0][2]*O[0][2]);
O[0][0]=O[0][0]/len;O[0][1]=O[0][1]/len;O[0][2]=O[0][2]/len;

O[1][0]=O[0][1]*U[2]-O[0][2]*U[1];
O[1][1]=O[0][2]*U[0]-O[0][0]*U[2];
O[1][2]=O[0][0]*U[1]-O[0][1]*U[0];

len=sqrt(O[1][0]*O[1][0]+O[1][1]*O[1][1]+O[1][2]*O[1][2]);
O[1][0]=O[1][0]/len;O[1][1]=O[1][1]/len;O[1][2]=O[1][2]/len;

O[2][0]=O[0][1]*O[1][2]-O[0][2]*O[1][1];
O[2][1]=O[0][2]*O[1][0]-O[0][0]*O[1][2];
O[2][2]=O[0][0]*O[1][1]-O[0][1]*O[1][0];

//cout << O[0][0]*O[1][0]+O[0][1]*O[1][1]+O[0][2]*O[1][2] <<endl;
//cout << O[0][0]*O[2][0]+O[0][1]*O[2][1]+O[0][2]*O[2][2] <<endl;
//cout << O[2][0]*O[1][0]+O[2][1]*O[1][1]+O[2][2]*O[1][2] <<endl;

//cout << O[0][0]*O[0][0]+O[0][1]*O[0][1]+O[0][2]*O[0][2] <<endl;
//cout << O[1][0]*O[1][0]+O[1][1]*O[1][1]+O[1][2]*O[1][2] <<endl;
//cout << O[2][0]*O[2][0]+O[2][1]*O[2][1]+O[2][2]*O[2][2] <<endl;

return;

}

void Camera::get_basis(double o[],double e1[],double e2[],double e3[]){
	
	e1[0]=O[0][0]; e1[1]=O[0][1]; e1[2]=O[0][2];
	e2[0]=O[1][0]; e2[1]=O[1][1]; e2[2]=O[1][2];
	e3[0]=O[2][0]; e3[1]=O[2][1]; e3[2]=O[2][2];
	o[0]=p[0];o[1]=p[1];o[2]=p[2];

}



Pinhole_Camera::Pinhole_Camera(double P[], double L[], double U[], double F)
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

// gets the rays in a way that is useful for the geola setup.  In particular, we only apply the field of view in the horizontal direction, 
// whereas we use h_m and h_M for the vertical direction.  We also assume the camera translates along the y direction and looks straight at the xz plane.
void Pinhole_Camera::get_rays_geola(int i, int j, const int w, const int h, vector< vector<double> >* P,vector< vector<double> >* V){

	int num_rays = 3;

	srand( time(NULL) );
	
	double i_r,j_r; //these will hold small random perturbations in position index

	P->clear();V->clear();
	vector<double> P1,V1;
	//in this case P and V are 1 by 3
	
	// I'm going to assume the camera translates along the y direction and has a location with a positive x value.
	double y_max = (p[0]-2)*tan(field_of_view/2.0);
	//cout<<" p[0] "<<p[0]<<endl;
	//cout<<"field_of_view = "<<field_of_view<<endl;
	//cout<<"field_of_view/2.0 = "<<field_of_view/2.0<<endl;
	//cout<<"tan(field_of_view/2.0) = "<<tan(field_of_view/2.0)<<endl;
	//cout<<" y_max = "<<y_max<<endl;

	double h_m=-3;
	double h_M=3;
	double h_c = 0;

	for(int ii=0;ii<num_rays;ii++){
		for(int jj=0;jj<num_rays;jj++){
			//P is easy, it's just the camera location
			for(int k=0;k<3;k++){P1.push_back(p[k]);}
			P->push_back(P1);

			i_r=((rand() % 100)-50)/100.0; //i.e. rand(-0.5,0.5)
			j_r=((rand() % 100)-50)/100.0; //i.e. rand(-0.5,0.5)

			//cout<<"float(ii+i_r)/num_rays : "<<float(ii+i_r)/num_rays<<" float(jj+j_r)/num_rays : "<<float(jj+j_r)/num_rays<<endl;
			//V is more complicated
			//double x=field_of_view*((float(i)+float(ii+i_r)/num_rays)/w-0.5);
			//double y=field_of_view*((float(j)+float(jj+j_r)/num_rays)/h-0.5);
			double x = -y_max+2*y_max*((float(i)+float(ii+i_r)/num_rays)/w-0.0);
			double y = h_m+(h_M-h_m)*((float(j)+float(jj+j_r)/num_rays)/h-0.0);
			//double x=field_of_view*((float(i))/w-0.5);
			//double y=field_of_view*((float(j))/h-0.5);
			//cout<<"x y = "<<x<<" "<<y<<endl;
			double ray[3];
			for(int k=0;k<3;k++){
				ray[k]=(p[0]-2)*O[0][k]+x*O[1][k]+y*O[2][k];
			}
			double length = sqrt(ray[0]*ray[0]+ray[1]*ray[1]+ray[2]*ray[2]);
			ray[0]=ray[0]/length;
			ray[1]=ray[1]/length;
			ray[2]=ray[2]/length;
			for(int k=0;k<3;k++){
				//V1.push_back(cos(theta)*O[0][k]+sin(theta)*(cos(phi)*O[1][k]+sin(phi)*O[2][k]));
				V1.push_back(ray[k]);
			}
			V->push_back(V1);
			P1.clear();
			V1.clear();
		}
	}
}

// gets the rays in a way that is useful for the chimera setup.  
void Pinhole_Camera::get_rays_chimera(int i, int j, const int w, const int h, vector< vector<double> >* P,vector< vector<double> >* V){

	int num_rays = 3;

	srand( time(NULL) );
	
	double i_r,j_r; //these will hold small random perturbations in position index

	P->clear();V->clear();
	vector<double> P1,V1;
	//in this case P and V are 1 by 3
	
	double buffer = 1.1;
	
	// I'm going to assume the camera translates along the y direction and has a location with a positive x value.
	double y_max = 2.25*buffer;
	//cout<<" p[0] "<<p[0]<<endl;
	//cout<<"field_of_view = "<<field_of_view<<endl;
	//cout<<"field_of_view/2.0 = "<<field_of_view/2.0<<endl;
	//cout<<"tan(field_of_view/2.0) = "<<tan(field_of_view/2.0)<<endl;
	//cout<<" y_max = "<<y_max<<endl;
	
	double angle = (PI/180)*24.0;
	
	double h_m=-2.5/cos(angle);
	double h_M=2.5/cos(angle);
	
	y_max = 0.5*(h_M-h_m)*(3.0/4)*buffer;
	
	double h_c = 0;
	double h_mp = (h_m-h_c)*buffer+h_c;
	double h_Mp = (h_M-h_c)*buffer+h_c;

	//double r = 9;
	double r = 9;

	for(int ii=0;ii<num_rays;ii++){
		for(int jj=0;jj<num_rays;jj++){
			//P is easy, it's just the camera location
			for(int k=0;k<3;k++){P1.push_back(p[k]);}
			P->push_back(P1);

			i_r=((rand() % 100)-50)/100.0; //i.e. rand(-0.5,0.5)
			j_r=((rand() % 100)-50)/100.0; //i.e. rand(-0.5,0.5)

			//cout<<"float(ii+i_r)/num_rays : "<<float(ii+i_r)/num_rays<<" float(jj+j_r)/num_rays : "<<float(jj+j_r)/num_rays<<endl;
			//V is more complicated
			//double x=field_of_view*((float(i)+float(ii+i_r)/num_rays)/w-0.5);
			//double y=field_of_view*((float(j)+float(jj+j_r)/num_rays)/h-0.5);
			double x = -y_max+2*y_max*((float(i)+float(ii+i_r)/num_rays)/w-0.0);
			double y = h_mp+(h_Mp-h_mp)*((float(j)+float(jj+j_r)/num_rays)/h-0.0);
			//double x=field_of_view*((float(i))/w-0.5);
			//double y=field_of_view*((float(j))/h-0.5);
			//cout<<"x y = "<<x<<" "<<y<<endl;
			double ray[3];
			for(int k=0;k<3;k++){
				ray[k]=r*O[0][k]+x*O[1][k]+y*O[2][k];
			}
			double length = sqrt(ray[0]*ray[0]+ray[1]*ray[1]+ray[2]*ray[2]);
			ray[0]=ray[0]/length;
			ray[1]=ray[1]/length;
			ray[2]=ray[2]/length;
			for(int k=0;k<3;k++){
				//V1.push_back(cos(theta)*O[0][k]+sin(theta)*(cos(phi)*O[1][k]+sin(phi)*O[2][k]));
				V1.push_back(ray[k]);
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

}
