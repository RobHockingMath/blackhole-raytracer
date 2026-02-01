#include <math.h>
#include <iostream>
#include "Mandelbulb.h"
#include "Sphere.h"
#include "invisibility.h"

#define PI 3.14159

using namespace std;

Mandelbulb::Mandelbulb(double c1,double c2,double c3,double E,int N,int red,int green, int blue)
	: Distance_fractal(c1,c2,c3,E,red,green,blue)
{
   power=N;
   is_wireframe=false;
}

bool Mandelbulb::wireframe_condition(double p[]){
	return false;
}

//x2 is overwritten with x1*x2.  x1 is unchanged
void Mandelbulb::Mult(double x1[],double x2[]){

	double p1=sqrt(x1[0]*x1[0]+x1[1]*x1[1]);
	double p2=sqrt(x2[0]*x2[0]+x2[1]*x2[1]);
	double gets_reused=(1-x1[2]*x2[2]/(p1*p2));
	double temp[3];
	temp[0]=(x1[0]*x2[0]-x1[1]*x2[1])*gets_reused;
	temp[1]=(x2[0]*x1[1]+x1[0]*x2[1])*gets_reused;
	temp[2]=p1*x2[2]+p2*x1[2];
	for(int i=0;i<3;i++){x2[i]=temp[i];}
}

//x unchanged, y is x^n
void Mandelbulb::raise_to_power(double x[], int n,double y[]){
	y[0]=1;
	y[1]=0;
	y[2]=0;
	for(int i=0;i<n;i++){
		Mult(x,y);
	}
}
//x is replaced with x^n
void Mandelbulb::raise_to_power_overwrite(double x[], int n){
	double y[3]={1,0,0};
	for(int i=0;i<n;i++){
		Mult(x,y);
	}
	x[0]=y[0];x[1]=y[1];x[2]=y[2];
}

//two sequences x_n and x'_n are relavent for the distance estimator
//this updates them both.
void Mandelbulb::Next(double x[],double xp[],double x0[]){
	//formula are x'_n+1=power*x_n^(power-1)*x'_n+1,
	//  x_n+1=x_n^power+x_0
	double y[3]={1,0,0};
	raise_to_power(x,power-1,y);
	raise_to_power_overwrite(x,power);
	Mult(y,xp);
	
	for(int i=0;i<3;i++){
		x[i]+=x0[i];
		xp[i]*=power;
		if(i==0){
			xp[i]+=1;
		}
	}
	

}

double Mandelbulb::Distance(double p[]){
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
	
	
	for(int i=0;i<3000;i++){
		powR=power*pow(r,pd);
		powRsin=powR*r_dz*sin(ph_dz+pd*ph);
		dz[0]=powRsin*cos(th_dz+pd*th)+1;
		dz[1]=powRsin*sin(th_dz+pd*th);
		dz[2]=powR*r_dz*cos(ph_dz+pd*ph);

		r_dz=pow(fabs(dz[0])+fabs(dz[1])+fabs(dz[2]),3);
		th_dz=atan2(dz[1],dz[0]);
		ph_dz=acos(dz[2]/r_dz);
		
		powR=pow(r,power);
		powRsin=sin(power*ph);
		z[0]=powR*powRsin*cos(power*th);
		z[1]=powR*powRsin*cos(power*th);
		z[2]=powR*sin(power*ph);
		z[0]+=p[0];z[1]+=p[1];z[2]+=p[2];
		
		r=pow(fabs(z[0])+fabs(z[1])+fabs(z[2]),3);
		if(r>500){break;}
		
		th=atan2(z[1],z[0]);
		ph=acos(z[2]/r);
	}
	return 0.5*r*log(r)/r_dz;
}

/**
//returns a lower bound for the distance from point p to the set
double Mandelbulb::Distance(double p[]){

	//use magic formula
	//d(q)=lim n->inf |qn| log |qn| /(2 |q'n|)

	double q[3]={p[0],p[1],p[2]};
	double qp[3]={1,0,0};
	for(int i=0;i<100;i++){
		Next(q,qp,p);
		if(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]>100){break;}
	}
	
	double len_q=sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
	double len_qp=sqrt(qp[0]*qp[0]+qp[1]*qp[1]+qp[2]*qp[2]);

	//cout<<"here "<<len_qp<<" "<<len_q<<endl;
	if(len_q!=0 && len_qp!=0){
		return len_q*log(len_q)/(2*len_qp);}
	else{return 10;}
	

}
**/
