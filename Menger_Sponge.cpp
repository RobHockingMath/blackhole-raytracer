#include <math.h>
#include <iostream>
#include "Menger_Sponge.h"
#include "Sphere.h"
#include "invisibility.h"

#define PI 3.14159
#define SQRT3 1.7321

using namespace std;

Menger_Sponge::Menger_Sponge(double c1,double c2,double c3,double scale_,int depth_,double shift_,double E,int red,int green, int blue)
	: Distance_fractal(c1,c2,c3,E,red,green,blue)
{
   depth = depth_;
   scale = scale_;
   shift= shift_;
   is_wireframe=true;
   rat = false;
}


bool Menger_Sponge::wireframe_condition(double p[]){

	double x=p[0], y=p[1], z=p[2];
	
	x=x/scale;y=y/scale;z=z/scale; //center it by changing position and scale


	x=x-c[0];
	y=y-c[1];
	z=z-c[2];

	if(!rat){
		double n = 1.0/SQRT3;
		
		double d_plane = (x-shift/3)*n+(y-shift/3)*n+(z-shift/3)*n;  //(p-1*c).n = dist
		
		double h = 0.01;
		/**d_plane = d_plane-h;
		if(d_plane<0 && d_plane>-2*h){d_plane=0;}
		else if(d_plane<0){
			d_plane = 2*h+d_plane;
		}**/
		
		//double delta = 0.02*2.0/3;
		double delta = 0.04*2.0/3;
		//double delta = 0.02;
		//if(abs(d_plane)<delta && (abs(x)<delta || abs(x-1) <delta || abs(y)<delta || abs(y-1) <delta || abs(z)<delta || abs(z-1) <delta)){
		int edge_variable = 0;
		if(abs(x)<delta || abs(1-x)<delta){
			edge_variable++;
		}
		if(abs(y)<delta || abs(1-y)<delta){
			edge_variable++;
		}
		if(abs(z)<delta || abs(1-z)<delta){
			edge_variable++;
		}
		if(edge_variable>=2){
			return true;
		}
		else{
			return false;
		}
	}
	else{
		double xx=abs(x-0.5)-0.5, yy=abs(y-0.5)-0.5, zz=abs(z-0.5)-0.5;
		
		double d1 = min(max(xx,zz),min(max(xx,yy),max(yy,zz)));
		
		double n = 1.0/SQRT3;
		double d_plane = (x-shift/3)*n+(y-shift/3)*n+(z-shift/3)*n;
		if(abs(d_plane)<0.05 && abs(d1)<0.05){return true;}
		else{return false;}
		/**double dx = abs(x); double dy=abs(y); double dz=abs(z);
		double dx1 = abs(x-1); double dy1 = abs(y-1); double dz1=abs(z-1);
		double delta = 0.02*2.0/3;
		//cout<<" dx = "<<dx<<endl;
		//if((dx<delta && dy<delta) || (dx<delta && dz<delta) || (dy<delta && dz<delta)){
		if(  (dx<delta && (min(dy,dy1)<delta || min(dz,dz1)<delta)) || (dy<delta && (min(dx,dx1)<delta || min(dz,dz1)<delta)) || (dz<delta && (min(dy,dy1)<delta || min(dx,dx1)<delta))  ){
			//cout<<" rat = "<<rat<<endl;
			return true;
		}
		else{
			return false;
		}**/
		
	}

}


//returns a lower bound for the distance from point p to the set
double Menger_Sponge::Distance(double p[]){

	double x=p[0], y=p[1], z=p[2];
	
	x=x/scale;y=y/scale;z=z/scale; //center it by changing position and scale

	double n = 1.0/SQRT3;
	
	double d_plane = (x-shift/3)*n+(y-shift/3)*n+(z-shift/3)*n;  //(p-1*c).n = dist
	if(is_wireframe && !rat){
		double h = 0.01;
		/**d_plane = d_plane-h;
		if(d_plane<0 && d_plane>-2*h){d_plane=0;}
		else if(d_plane<0){
			d_plane = 2*h+d_plane;
		}**/
		d_plane=d_plane;
	}
	else{
		d_plane=-d_plane;
	}
	
	
	

	//x = (1+x)/2;y = (1+y)/2; z = (1+z)/2;
	
	

	double xx=abs(x-0.5)-0.5, yy=abs(y-0.5)-0.5, zz=abs(z-0.5)-0.5;
	double d1;
	/**if(xx>=0 || yy>=0 || zz>=0){
		d1=max(xx,max(yy,zz));
	}
	else{
		d1=min(xx,min(yy,zz));
	}**/
	d1 = max(xx,max(yy,zz));
	//d1 = min(max(xx,zz),min(max(xx,yy),max(yy,zz)));
	
	d1 = max(d1,d_plane);
	
	
	double d = d1;
	double pp=1.0;
	for (int i=1; i<=depth; ++i) {
		double xa = fmod(3.0*x*pp,3.0);
		double ya = fmod(3.0*y*pp,3.0);
		double za = fmod(3.0*z*pp,3.0);
		pp*=3.0;

		//we can also translate/rotate (xa,ya,za) without affecting the DE estimate

		double xx=0.5-abs(xa-1.5), yy=0.5-abs(ya-1.5), zz=0.5-abs(za-1.5);
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
		d1=min(max(xx,zz),min(max(xx,yy),max(yy,zz))) / pp; //distance inside the 3 axis-aligned square tubes
		//d1 = min(dxy,min(dxz,dyz));
		/**if(dxy>=0 || dxz>=0 || dyz>=0){
			d1 = min(dxy,min(dxz,dyz));
		}
		else{
			d1 = min(dxy,min(dxz,dyz));
		}**/
		//d1 = -d1;
		d=max(d,d1); //intersection
	}
	d = d;
	//d = max(d,d_plane);
	//return d*2.0; //the distance estimate. The *2 is because of the scaling we did at the beginning of the function
	return d*scale;
	

}


Upright_Cylinder_Infinite::Upright_Cylinder_Infinite(double zt, double zb,double x0, double y0, double r, double E, int red,int green, int blue)
	: Distance_fractal(0,0,0,E,red,green,blue)
{
   x_0 = x0;
   y_0 = y0;
   z_b = zb;
   z_t = zt;
   R= r;
   is_wireframe=false;
   has_upper_bound=false;
   has_lower_bound=false;
   upper_bound=0;
   lower_bound=0;
}

double Upright_Cylinder_Infinite::Distance(double p[]){

	
	double dist = fabs(sqrt((p[0]-x_0)*(p[0]-x_0)+(p[1]-y_0)*(p[1]-y_0))-R);
	if(!has_upper_bound && !has_lower_bound){
		return dist;
	}
	else if(has_upper_bound && !has_lower_bound){
		if(p[2]>upper_bound){
			return fmax(p[2]-upper_bound,dist);
		}
	}
	else if(!has_upper_bound && has_lower_bound){
		if(p[2]<lower_bound){
			return fmax(lower_bound-p[2],dist);
		}
	}
	else if(has_upper_bound && has_lower_bound){
		if(p[2]>upper_bound){
			dist = fmax(p[2]-upper_bound,dist);
		}
		if(p[2]<lower_bound){
			return fmax(lower_bound-p[2],dist);
		}
	
	}
	return dist;

}

bool Upright_Cylinder_Infinite::wireframe_condition(double p[]){
	return true;
}

void Upright_Cylinder_Infinite::inverse_texture(int i,int j,int n, int m,float query_pt[]){
	double u = i/double(n);
	double v = j/double(m);
	double x = x_0-R*cos(2*PI*u);
	double y = y_0-R*sin(2*PI*u);
	double z = z_b + (z_t-z_b)*v;
	query_pt[0]=x;
	query_pt[1]=y;
	query_pt[2]=z;
	return;
}

void Upright_Cylinder_Infinite::get_Texture_Coords(double x[], double u[]){
	double phi = atan2(x[1]-y_0,x[0]-x_0);
	u[0]=(atan2(x[1]-y_0,x[0]-x_0)+PI)/(2*PI);
	u[1]=(x[2]-z_b)/(z_t-z_b);
	u[0]=u[0]-floor(u[0]/1)*1;
	u[1]=u[1]-floor(u[1]/1)*1;
	return;
}

void Upright_Cylinder_Infinite::get_Shading(double x[],vector<Texture*> textures,int obj_id,int C[]){
	double u[2]={0,0};
	this->get_Texture_Coords(x,u);
	textures[obj_id]->bilinear_get_at(u,C);
}