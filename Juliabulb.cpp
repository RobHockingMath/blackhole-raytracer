#include <math.h>
#include <iostream>
#include "Juliabulb.h"
#include "Sphere.h"
#include "invisibility.h"

#define PI 3.14159

using namespace std;

//c1,c2,c3 are xyz of center.
//c_1,c_2,c_3 are components of "c" in z'=z^n+c.

Juliabulb::Juliabulb(double c1,double c2,double c3,double E,int N,double c_1,double c_2,double c_3,int red,int green, int blue)
	: Distance_fractal(c1,c2,c3,E,red,green,blue)
{
   power=N;
   c[0]=c_1;c[1]=c_2;c[2]=c_3; 
   is_wireframe=false;
}

bool Juliabulb::wireframe_condition(double p[]){
	return false;
}

double Juliabulb::Distance(double p[]){
	double pd=power-1;
	
	double z[3]={p[0],p[1],p[2]};
	
	double r=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
	double th=atan2(p[1],p[0]);
	double ph=asin(p[2]/r);
	
	double dz[3]={0,0,0};
	double ph_dz=0;
	double th_dz=0;
	double r_dz=1;
	double powR,powRsin;
	
	
	for(int i=0;i<500;i++){
		powR=power*pow(r,pd);
		powRsin=powR*r_dz*sin(ph_dz+pd*ph);
		dz[0]=powRsin*cos(th_dz+pd*th);
		dz[1]=powRsin*sin(th_dz+pd*th);
		dz[2]=powR*r_dz*cos(ph_dz+pd*ph);

		r_dz=sqrt(dz[0]*dz[0]+dz[1]*dz[1]+dz[2]*dz[2]);
		th_dz=atan2(dz[1],dz[0]);
		ph_dz=acos(dz[2]/r_dz);
		
		powR=pow(r,power);
		powRsin=sin(power*ph);
		z[0]=powR*powRsin*cos(power*th);
		z[1]=powR*powRsin*sin(power*th);
		z[2]=powR*cos(power*ph);
		z[0]+=c[0];z[1]+=c[1];z[2]+=c[2];
		
		r=sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);
		if(r>500){break;}
		
		th=atan2(z[1],z[0]);
		ph=acos(z[2]/r);
	}
	return 0.5*r*log(r)/r_dz;
}
