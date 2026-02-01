#include <math.h>
#include <iostream>
#include "Menger_Sponge_Type2.h"
#include "Sphere.h"
#include "invisibility.h"

#define PI 3.14159

using namespace std;

Menger_Sponge_Type2::Menger_Sponge_Type2(double c1,double c2,double c3,double scale_,int depth_,double E,int red,int green, int blue)
	: Distance_fractal(c1,c2,c3,E,red,green,blue)
{
   scale = scale_;
   depth=depth_;
}


//returns a lower bound for the distance from point p to the set
double Menger_Sponge_Type2::Distance(double p[]){

	double x=p[0], y=p[1], z=p[2];
	x=x/scale;y=y/scale;z=z/scale; //center it by changing position and scale
	double cx=0.5;double cy=0.5;double cz=0.5;
	double h = 1.0/3;
	double xx=abs(x-cx)-h/2;
	double yy=abs(y-cy)-h/2;
	double zz=abs(z-cz)-h/2;
	double dist = max(xx,max(yy,zz));
	h = 1.0/9;
	for(int i=-1;i<=1;i++){
		if(i==-1){
			cx=1.0/6;
		}
		if(i==0){
			cx=1.0/2;
		}
		if(i==1){
			cx=1-1.0/6;
		}
		for(int j=-1;j<=1;j++){
			if(j==-1){
				cy=1.0/6;
			}
			if(j==0){
				cy=1.0/2;
			}
			if(j==1){
				cy=1-1.0/6;
			}
			for(int k=-1;k<=1;k++){
				if(k==-1){
					cz=1.0/6;
				}
				if(k==0){
					cz=1.0/2;
				}
				if(k==1){
					cz=1-1.0/6;
				}
				if(!(i==0 && j==0 && k==0)){
					xx=abs(x-cx)-h/2;
					yy=abs(y-cy)-h/2;
					zz=abs(z-cz)-h/2;
					dist = min(dist,max(xx,max(yy,zz)));
				}
			}
		}
	}
	return dist;
	

}
