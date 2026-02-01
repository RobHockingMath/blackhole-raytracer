
#include <math.h>
#include <iostream>
using namespace std;
#include "Light.h"

#define PI 3.14159265359

Light::Light(double b,int C[]){
	color[0]=C[0];
	color[1]=C[1];
	color[2]=C[2];
	brightness=b;
	illumination_map_ID=-1;
	pp=2;
	is_flat = true;
}

planeLight::planeLight(double Z, double b, int C[],Texture* my_Texture)
	: Light(b,C)
{
	z=Z;
	texture = my_Texture;
}

void planeLight::get_Contribution(double n[],double x[],double v[],vector<Texture*> textures,double R,double frac[]){

	double d;
	d=fabs(z-x[2]);
	double f=(brightness)*(n[2]);
	
	int col[3]={0,0,0};
	this->getColor(x[0],x[1],col);
	
	if(f<0){f=0;}
	if(f>1){f=1;}
	for(int i=0;i<3;i++){frac[i]=double(f*col[i])/255.0;}
}

void planeLight::getColor(double x, double y, int col[3]){
	if(use_texture){
		double u[2]={0,0};
		u[0] = (x-(-3))/(2*3.84);
		u[1] = (y-(-2.5))/(2*2.16);
		
		u[0]=u[0]-floor(u[0]/1)*(1);
		u[1]=u[1]-floor(u[1]/1)*(1);
		
		int w=texture->get_width()-1;
		int h=texture->get_height()-1;
		int i=int(u[0]*w);int j=int((1-u[1])*h);
		texture->get_at(i,j,col);
		//texture->bilinear_get_at(u, col);
	}
	else{
		col[0]=color[0];
		col[1]=color[1];
		col[2]=color[2];
	}

	
	
}

void planeLight::get_Specular_Contribution(double r[],double x[],double specular_term[]){
	
	double len=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
	r[0]/=len; r[1]/=len; r[2]/=len;
	int col[3]={0,0,0};
	this->getColor(x[0],x[1],col);
	
	
	double specular = pow(fabs(1*r[2]),25.0);
	specular_term[0]=100*specular*col[0]/255;
	specular_term[1]=100*specular*col[1]/255;
	specular_term[2]=100*specular*col[2]/255;
}


bool planeLight::is_Visible(double x[],double R,double dp[],vector<Shape*> shapes){
	//double V[3]={pos[0]-x[0],pos[1]-x[1],pos[2]-x[2]};
	double V[3]={0,0,x[2]-z};
	double dist = sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
	V[0]=V[0]/dist;V[1]=V[1]/dist;V[2]=V[2]/dist;
	double t_i1,t_i2,y[4];
	double eps = 0.01;
	double pos1[3]={x[0],x[1],z};
	double pos2[3]={x[0],x[1],z};
	pos1[0]+=V[0]*eps;
	pos1[1]+=V[1]*eps;
	pos1[2]+=V[2]*eps;
	pos2[0]-=V[0]*eps;
	pos2[1]-=V[1]*eps;
	pos2[2]-=V[2]*eps;
	//if(R>0){get_direction(pos,x,R,V);}	

	for(int i=0; i<shapes.size();i++){
		if(shapes[i]->is_wireframe==false){
			if(R==0){t_i1=shapes[i]->get_Intersection(pos1,V,y);t_i2=shapes[i]->get_Intersection(pos2,V,y);}
			//else{t_i1=shapes[i]->get_Intersection2(pos,V,R,y);}
			//things are scaled so that t=0 -> x,
			//t=1 -> pos.  Only things between x and pos
			//cast shadows, hence the t_i<1.
			if(t_i1>0.000*dist && t_i1<0.9999*dist && t_i2>0.000*dist && t_i2<0.9999*dist){return false;}
		}
	}
	return true;
}

double planeLight::percent_Visible(double x[],double R,vector<Shape*> shapes){
	double dp[3]={0,0,0};
	if(this->is_Visible(x,R,dp,shapes)){
		return 1.0;
	}
	else{
		return 0.0;
	}
	
}


pointLight::pointLight(double P[], double b, int C[])
	: Light(b,C)
{
	pos[0]=P[0];pos[1]=P[1];pos[2]=P[2];
}

void pointLight::color_from_illumination_map(double X[],double n[],vector<Texture*> textures, int C[]){
	/**if(illumination_map_ID<0){
		C[0]=color[0];
		C[1]=color[1];
		C[2]=color[2];
	}
	else{
		double r = sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
		double u[2]={0,0};
		u[0]=(atan2(V[1],V[0])+PI)/(2*PI);
		u[1]=1-acos(V[2]/r)/PI;
		
		// periodicity
		u[0]=2*u[0];
		u[0]=u[0]-floor(u[0]);
		
		if(u[0]<0){u[0]=0;}
		if(u[1]<0){u[1]=0;}
		if(u[0]>1){u[0]=1;}
		if(u[1]>1){u[1]=1;}
		
		int W,H;
		
		//# pragma omp critical
		{
			W=textures[illumination_map_ID]->get_width();
			H=textures[illumination_map_ID]->get_height();
		}
		
		int ii=floor(u[0]*float(W-1));
		int ip=min(ii+1,W-1);
		int jj=floor(u[1]*float(H-1));
		int jp=min(jj+1,H-1);
		
		int C00[3]={C[0],C[1],C[2]};
		int C01[3]={C[0],C[1],C[2]};
		int C10[3]={C[0],C[1],C[2]};
		int C11[3]={C[0],C[1],C[2]};
		
		//# pragma omp critical
		{
			textures[illumination_map_ID]->get_at(ii,jj,C00);
			textures[illumination_map_ID]->get_at(ip,jj,C10);
			textures[illumination_map_ID]->get_at(ii,jp,C01);
			textures[illumination_map_ID]->get_at(ip,jp,C11);
		}
		
		double s = u[0]*float(W-1)-ii;
		double t = u[1]*float(H-1)-jj;
		
		C[0]=int((1-s)*(1-t)*C00[0]+s*(1-t)*C10[0]+(1-s)*t*C01[0]+s*t*C11[0]);
		C[1]=int((1-s)*(1-t)*C00[1]+s*(1-t)*C10[1]+(1-s)*t*C01[1]+s*t*C11[1]);
		C[2]=int((1-s)*(1-t)*C00[2]+s*(1-t)*C10[2]+(1-s)*t*C01[2]+s*t*C11[2]);
		
		return;
	}**/
	
	if(!is_flat){
	
		double u[2]={0,0};
		
		double n12[3]={0.4714,0.8165,0.3333};
		double e12[3]={-0.8660,0.5000,0};
		double u12[3]={-0.1667,-0.2887,0.9428};
		double d12 = (n[0]-n12[0])*(n[0]-n12[0])+(n[1]-n12[1])*(n[1]-n12[1])+(n[2]-n12[2])*(n[2]-n12[2]);
		double d12m = (n[0]+n12[0])*(n[0]+n12[0])+(n[1]+n12[1])*(n[1]+n12[1])+(n[2]+n12[2])*(n[2]+n12[2]);
		
		double n23[3]={-0.9428,0,0.3333};
		double e23[3]={0,-1,0};
		double u23[3]={0.3333,0,0.9428};
		double d23 = (n[0]-n23[0])*(n[0]-n23[0])+(n[1]-n23[1])*(n[1]-n23[1])+(n[2]-n23[2])*(n[2]-n23[2]);
		double d23m = (n[0]+n23[0])*(n[0]+n23[0])+(n[1]+n23[1])*(n[1]+n23[1])+(n[2]+n23[2])*(n[2]+n23[2]);
		
		double n31[3]={0.4714,-0.8165,0.3333};
		double e31[3]={0.8660,0.5000,0};
		double u31[3]={-0.1667,0.2887,0.9428};
		double d31 = (n[0]-n31[0])*(n[0]-n31[0])+(n[1]-n31[1])*(n[1]-n31[1])+(n[2]-n31[2])*(n[2]-n31[2]);
		double d31m = (n[0]+n31[0])*(n[0]+n31[0])+(n[1]+n31[1])*(n[1]+n31[1])+(n[2]+n31[2])*(n[2]+n31[2]);
		
		double n4[3]={0,0,1};
		double e4[3]={1,0,0};
		double u4[3]={0,1,0};
		double d4 = (n[0]-n4[0])*(n[0]-n4[0])+(n[1]-n4[1])*(n[1]-n4[1])+(n[2]-n4[2])*(n[2]-n4[2]);
		double d4m = (n[0]+n4[0])*(n[0]+n4[0])+(n[1]+n4[1])*(n[1]+n4[1])+(n[2]+n4[2])*(n[2]+n4[2]);
		
		double d_min = min( min(min(d12,d23),min(d31,d4)), min(min(d12m,d23m),min(d31m,d4m)) );
		
		double L=0.5;
		
		if(d_min == d12 || d_min == d12m){
			u[0]=(X[0]*e12[0]+X[1]*e12[1]+X[2]*e12[2]-(-0.5*L))/L;
			u[1]=(X[0]*u12[0]+X[1]*u12[1]+X[2]*u12[2]-(-0.5*L))/L;
		}
		else if(d_min == d23 || d_min == d23m){
			u[0]=(X[0]*e23[0]+X[1]*e23[1]+X[2]*e23[2]-(-0.5*L))/L;
			u[1]=(X[0]*u23[0]+X[1]*u23[1]+X[2]*u23[2]-(-0.5*L))/L;
		}
		else if(d_min == d31 || d_min == d31m){
			u[0]=(X[0]*e31[0]+X[1]*e31[1]+X[2]*e31[2]-(-0.5*L))/L;
			u[1]=(X[0]*u31[0]+X[1]*u31[1]+X[2]*u31[2]-(-0.5*L))/L;
		}
		else if(d_min == d4 || d_min == d4m){
			u[0]=(X[0]-(-0.5*L))/L;
			u[1]=(X[1]-(-0.5*L))/L;
		}
		u[0]=u[0]-floor(u[0]/1)*1;
		u[1]=u[1]-floor(u[1]/1)*1;
		
		textures[illumination_map_ID]->bilinear_get_at(u,C);
	
	}
	else{
	
		double L=1;
	
		double u[2]={0,0};
		
		double V[3]={pos[0]-X[0],pos[1]-X[1],pos[2]-X[2]};
		double t,xx[3];
		if(V[2]!=0){
			t = (2-X[2])/V[2];
			xx[0]=X[0]+t*V[0];
			xx[1]=X[1]+t*V[1];
			u[0]=(xx[0]-(-6))/12;
			u[1]=(xx[1]-(-6))/12;
			u[0]=2*fabs(u[0]/2-floor(u[0]/2+0.5));
			u[1]=2*fabs(u[1]/2-floor(u[1]/2+0.5));
			textures[illumination_map_ID]->bilinear_get_at(u,C);
			C[0]=max(C[0]-40,0);
			C[1]=C[0];
			C[2]=C[0];
		}
		else{
			C[0]=color[0];
			C[1]=color[1];
			C[2]=color[2];
		}
		
		//plane at z=2
		
		//u[0]=(X[0]-0.5*L)/L;
		//u[1]=(X[1]-0.5*L)/L;
		//u[0]=u[0]-floor(u[0]/1)*1;
		//u[1]=u[1]-floor(u[1]/1)*1;
		//u[0]=2*fabs(u[0]/2-floor(u[0]/2+0.5));
		//u[1]=2*fabs(u[1]/2-floor(u[1]/2+0.5));
		
		//textures[illumination_map_ID]->bilinear_get_at(u,C);
		
		//double brightness = sqrt(C[0]*C[0]+C[1]*C[1]+C[2]*C[2]);
		//double max_brightness = sqrt(3)*255;
		//C[0]=int(floor((100*brightness/max_brightness)));
		//C[1]=int(floor((255*brightness/max_brightness)));
		//C[2]=int(floor((255*brightness/max_brightness)));
		//C[1]=255*brightness/cyan_brightness;
		
	}
	
}

void pointLight::get_Contribution(double n[],double x[],double v[],vector<Texture*> textures,double R,double frac[]){

	double d[3];
	d[0]=pos[0]-x[0];d[1]=pos[1]-x[1];d[2]=pos[2]-x[2];
	double len=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
	if(R>0){get_direction(x,pos,R,d);}
	double len2=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
	double f;
	if(fabs(pp-2)<0.1){
		f=(brightness/(len2*len2))*(n[0]*d[0]+n[1]*d[1]+n[2]*d[2])/len2;
	}
	else{
		f=brightness*(n[0]*d[0]+n[1]*d[1]+n[2]*d[2])/len2;
	}
	
	if(f<0){f=0;}
	if(f>1){f=1;}
	if(illumination_map_ID<0){
		for(int i=0;i<3;i++){frac[i]=double(f*color[i])/255.0;}
	}
	else{
		int C[3]={0,0,0};
		color_from_illumination_map(x,n,textures,C);
		for(int i=0;i<3;i++){frac[i]=double(f*C[i])/255.0;}
	}
}

void pointLight::get_Specular_Contribution(double r[],double x[],double specular_term[]){
	
	double len=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
	r[0]/=len; r[1]/=len; r[2]/=len;
	
	double d[3];
	d[0]=pos[0]-x[0];d[1]=pos[1]-x[1];d[2]=pos[2]-x[2];
	len=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
	d[0]/=len; d[1]/=len; d[2]/=len;
	double specular = pow(max(0.0,d[0]*r[0]+d[1]*r[1]+d[2]*r[2]),25.0);
	specular_term[0]=100*specular;
	specular_term[1]=100*specular;
	specular_term[2]=100*specular;
}

//void pointLight::get_Specular_Contribution(double n[], double x[], double v[]){
//	
//}

bool pointLight::is_Visible(double x[],double R,double dp[],vector<Shape*> shapes){
	//double V[3]={pos[0]-x[0],pos[1]-x[1],pos[2]-x[2]};
	double V[3]={x[0]-pos[0],x[1]-pos[1],x[2]-pos[2]};
	double dist = sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
	V[0]=V[0]/dist;V[1]=V[1]/dist;V[2]=V[2]/dist;
	double t_i1,t_i2,y[4];
	double eps = 0.01;
	double pos1[3]={pos[0],pos[1],pos[2]};
	double pos2[3]={pos[0],pos[1],pos[2]};
	pos1[0]+=V[0]*eps;
	pos1[1]+=V[1]*eps;
	pos1[2]+=V[2]*eps;
	pos2[0]-=V[0]*eps;
	pos2[1]-=V[1]*eps;
	pos2[2]-=V[2]*eps;
	if(R>0){get_direction(pos,x,R,V);}	

	dp[0]=V[0];
	dp[1]=V[1];
	dp[2]=V[2];

	for(int i=0; i<shapes.size();i++){
		if(shapes[i]->is_wireframe==false && shapes[i]->reflection_only==false){
			if(R==0){t_i1=shapes[i]->get_Intersection(pos1,V,y);t_i2=shapes[i]->get_Intersection(pos2,V,y);}
			else{t_i1=shapes[i]->get_Intersection2(pos,V,R,y);}
			//things are scaled so that t=0 -> x,
			//t=1 -> pos.  Only things between x and pos
			//cast shadows, hence the t_i<1.
			if(t_i1>0.000*dist && t_i1<0.9999*dist && t_i2>0.000*dist && t_i2<0.9999*dist){return false;}
		}
	}
	return true;
}

double pointLight::percent_Visible(double x[],double R,vector<Shape*> shapes){
	
	// we have to break the light into a series of points.
	double r0 = 0.5;
	double sum = 0.0;
	double max_sum = 0.0;
	for(int i=0;i<5;i++){
		double theta = 2*PI*i/5;
		for(int j=1;j<4;j++){
			double phi = PI*j/5;
			double dA = sin(phi);
			double dp[3] = {r0*cos(theta)*sin(phi),r0*sin(theta)*sin(phi),r0*cos(phi)};
			if(this->is_Visible(x,R,dp,shapes)){
				sum+=dA;
			}
			max_sum+=dA;
		}
	}
	return sum/max_sum;
}




ambientLight::ambientLight(double b, int C[])
	: Light(b,C)
{}

void ambientLight::get_Contribution(double n[],double x[],double v[],vector<Texture*> textures,double R,double frac[]){
	for(int i=0;i<3;i++){	
		frac[i]=brightness*color[i]/255;
		if(frac[i]>255){frac[i]=255;}
		if(frac[i]<0){frac[i]=0;}
	}
}

void ambientLight::get_Specular_Contribution(double r[],double x[],double specular_term[]){
	specular_term[0]=0;
	specular_term[1]=0;
	specular_term[2]=0;
}

bool ambientLight::is_Visible(double x[],double R,double dp[], vector<Shape*> shapes){
	return true;
}

double ambientLight::percent_Visible(double x[],double R, vector<Shape*> shapes){
	return 1;
}
