#include <math.h>
#include <iostream>
#include "Menger_Sponge_4DSlice_Type2.h"
#include "Sphere.h"
#include "invisibility.h"

#define PI 3.14159

using namespace std;
Menger_Sponge_4DSlice_Type2::Menger_Sponge_4DSlice_Type2(double c1,double c2,double c3,double scale_,int depth_,int type_,double E,double offset_,int red,int green, int blue)
	: Distance_fractal(c1,c2,c3,E,red,green,blue)
{
   depth = depth_;
   scale = scale_;
   offset = offset_;
   type = type_;
}


//returns a lower bound for the distance from point p to the set
double Menger_Sponge_4DSlice_Type2::Distance(double p[]){

	double X=p[0], Y=p[1], Z=p[2];
	
	X=X/scale;Y=Y/scale;Z=Z/scale; //center it by changing position and scale	
	
	/**double x = (X+Y-Z)/2.0;
	double y = (-X-Y-Z)/2.0;
	double z = (X-Y+Z)/2.0;
	double w = offset+(-X+Y+Z)/2.0;**/
	
	double x = 0.25*offset+(X+Y-Z)/2.0;
	double y = 0.25*offset+(-X-Y-Z)/2.0;
	double z = 0.25*offset+(X-Y+Z)/2.0;
	double w = 0.25*offset+(-X+Y+Z)/2.0;
	
	x = (1+x)/2;y = (1+y)/2; z = (1+z)/2;w=(1+w)/2;

	double dx,dy,dz,dw;
	if(0<=x & x<=1){
		dx = 0;
	}
	if(x>1){dx = x-1;}
	if(x<0){dx = -x;}

	if(0<=y & y<=1){
		dy = 0;
	}
	if(y>1){dy = y-1;}
	if(y<0){dy = -y;}

	if(0<=z & z<=1){
		dz = 0;
	}
	if(z>1){dz = z-1;}
	if(z<0){dz = -z;}

	if(0<=w & w<=1){
		dw = 0;
	}
	if(w>1){dw = w-1;}
	if(w<0){dw = -w;}	
	

	double d1 = max(Distance3D(x,y,z),dw);
	if(flip==false){
		d1 = min(d1,max(Distance3D(y,z,w),dx));
		d1 = min(d1,max(Distance3D(z,w,x),dy));
		d1 = min(d1,max(Distance3D(w,x,y),dz));
	}
	else{
		d1 = max(d1,max(Distance3D(y,z,w),dx));
		d1 = max(d1,max(Distance3D(z,w,x),dy));
		d1 = max(d1,max(Distance3D(w,x,y),dz));		
	}
	return d1;
	
	
	

}

double Menger_Sponge_4DSlice_Type2::Distance3D(double x,double y, double z){

	x=x/scale;y=y/scale;z=z/scale; //center it by changing position and scale
	/**double cx=0.5;double cy=0.5;double cz=0.5;
	double h = 1.0/3;
	double xx=abs(x-cx)-h/2;
	double yy=abs(y-cy)-h/2;
	double zz=abs(z-cz)-h/2;
	double dist = max(xx,max(yy,zz));
	double d = dist;
	h = 1.0/9;**/
	double xx=abs(x-0.5)-0.5, yy=abs(y-0.5)-0.5, zz=abs(z-0.5)-0.5;//ww=offset_p-x-y-z;//ww=abs(offset_p-x-y-z-0.5)-0.5;//ww=abs(w-0.5)-0.5;
	double d1;
	d1=max(xx,max(yy,zz)); //distance to the box
	double dist;
	if(flip==false){
		dist = 10000000;
	}
	else{
		dist= -10000000;
	}
	/**for(int i=-1;i<=1;i++){
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
	}**/
	double pp=1.0;
	for (int i=1; i<=depth; ++i) {
		double xa = fmod(3.0*x*pp,3.0);
		double ya = fmod(3.0*y*pp,3.0);
		double za = fmod(3.0*z*pp,3.0);
		pp*=3.0;

		//we can also translate/rotate (xa,ya,za) without affecting the DE estimate
		if(flip==false){
			xx=-(0.5-abs(xa-1.5)); yy=-(0.5-abs(ya-1.5)); zz=-(0.5-abs(za-1.5));
		}
		else{
			xx=(0.5-abs(xa-1.5)); yy=(0.5-abs(ya-1.5)); zz=(0.5-abs(za-1.5));
		}
		//double xx=abs(x-0.5)-0.5, yy=abs(y-0.5)-0.5, zz=abs(z-0.5)-0.5;
		/**double dxy,dxz,dyz;
		if(xx>=0 || yy>=0){
			dxy = max(xx,yy);
		}
		else{
			dxy = min(xx,yy);
		}
		if(xx>=0 || zz>=0){
			dxz = max(xx,zz);
		}
		else{
			dxz = min(xx,zz);
		}
		if(yy>=0 || zz>=0){
			dyz = max(yy,zz);
		}
		else{
			dyz = min(yy,zz);
		}**/
		if(flip==false){
			dist=min(dist,max(max(xx,yy),zz) / pp); //distance inside the 3 axis-aligned square tubes
		}
		else{
			dist=max(dist,-max(max(-xx,-yy),-zz) / pp); //distance inside the 3 axis-aligned square tubes
		}
		//dist=min(dist,min(max(xx,zz),min(max(xx,yy),max(yy,zz))) / pp); //w ommitted
		///d1 = max(d1,min(max(xx,ww),min(max(xx,yy),max(yy,ww))) / pp); // z ommitted
		//d1 = max(d1,min(max(xx,zz),min(max(xx,ww),max(ww,zz))) / pp); // y ommitted 
		//d1 = max(d1,min(max(ww,zz),min(max(ww,yy),max(yy,zz))) / pp); // x ommitted
		//d1 = min(dxy,min(dxz,dyz));
		/**if(dxy>=0 || dxz>=0 || dyz>=0){
			d1 = min(dxy,min(dxz,dyz));
		}
		else{
			d1 = min(dxy,min(dxz,dyz));
		}**/
		//d1 = -d1;
		//d=max(d,d1); //intersection
	}
	dist=max(dist,d1); //intersection
	return dist*scale;

}