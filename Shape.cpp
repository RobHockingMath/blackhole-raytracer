#include <math.h>
#include <iostream>
#include "Shape.h"


Shape::Shape(int red,int green, int blue){

   color[0]=red;color[1]=green;color[2]=blue;
   reflectivity=0; //non-reflective by default
   texture_ID=-1; //no texture by default
   bump_ID=-1; //no bumpmap by default
   transparency=0; //opaque by default
   ref_index=1; //no refraction by default
   is_wireframe=false; //default
   self_illuminating=false;
   reflection_only = false;
   alpha_ID=-1;
   one_sided=false;

}

void Shape::get_reflection(double x[], double v[]){
	double n[3];
	this->get_Normal(x,n);
	double dot=v[0]*n[0]+v[1]*n[1]+v[2]*n[2];
	v[0]=v[0]-2*dot*n[0];
	v[1]=v[1]-2*dot*n[1];
	v[2]=v[2]-2*dot*n[2];
}

//modifies v for refraction and returns true if internal 
//reflection takes place, false otherwise
bool Shape::get_refracted_ray(double x[],double v[]){

	if(ref_index==1){return false;}//no refraction

	double n[3];
	this->get_Normal(x,n);
	double v_dot_n,V,vt,Vt,n1,n2;
	v_dot_n=v[0]*n[0]+v[1]*n[1]+v[2]*n[2];
	//if v_dot_n<0, we are outside the object going in
	//otherwise we're already inside going out
	if(v_dot_n<0){n1=1;n2=ref_index;}
	else{n1=ref_index;n2=1;}
	V=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	vt=sqrt(V*V-v_dot_n*v_dot_n);

	if(vt==0){return false;}//parallel to normal
	
	double temp=(n1/n2)*(vt/V);
	if(temp*temp<=1){
		Vt=v_dot_n*temp/sqrt(1-temp*temp);
	
		v[0]=v[0]+(v[0]-v_dot_n*n[0])*(Vt/vt-1);
		v[0]=v[0]+(v[0]-v_dot_n*n[0])*(Vt/vt-1);
		v[0]=v[0]+(v[0]-v_dot_n*n[0])*(Vt/vt-1);
		return false;
	}
	else{ //total internal reflection
		get_reflection(x,v);
		return true;
	}
}

BasicShape::BasicShape(int red,int green,int blue)
	: Shape(red,green,blue)
{

}

void BasicShape::get_Color(double x[], vector<Texture*> textures, int col[]){
	double n[3];
	
	if(this->texture_ID>-1 && this->texture_ID<textures.size()){
		double u[2];
		int w=textures[texture_ID]->get_width()-1;
		int h=textures[texture_ID]->get_height()-1;
		this->get_Texture_Coords(x,u);
		//int i=int(u[0]*w);int j=int((1-u[1])*h);
		//textures[texture_ID]->get_at(i,j,col);
		textures[texture_ID]->bilinear_get_at(u,col);
	}
	else{
		col[0]=color[0];
		col[1]=color[1];
		col[2]=color[2];
	}
}

double BasicShape::get_Alpha(double x[], vector<Texture*> textures){
	if(!(fabs(x[0])<=3 && fabs(x[1])<=3) ){
		return 1.0;
	}
	int C[3]={0,0,0};
	this->get_Color(x,textures,C);
	return double(1-max(max(C[0],C[1]),C[2])/255.0);	
}

/**void BasicShape::get_Bump_Perturbation(double x[], vector<Texture*> textures, double dn[]){
	double u[2]={0,0};
	double u1[2]={0,0};
	double u2[2]={0,0};
	double u3[2]={0,0};
	
	double h = 0.05;
	
	double x1[3]={x[0]+h,x[1],x[2]};
	double x2[3]={x[0],x[1]+h,x[2]};
	double x3[3]={x[0],x[1],x[2]+h};
	
	this->get_Texture_Coords(x,u);
	this->get_Texture_Coords(x1,u1);
	this->get_Texture_Coords(x2,u2);
	this->get_Texture_Coords(x3,u3);
	double z0 = this->get_Bump_Value(u,textures);
	double z1 = this->get_Bump_Value(u1,textures);
	double z2 = this->get_Bump_Value(u2,textures);
	double z3 = this->get_Bump_Value(u3,textures);
	
	dn[0]=(z1-z0)/h;
	dn[1]=(z2-z0)/h;
	dn[2]=(z3-z0)/h;
	
	
}**/

/**double BasicShape::get_Bump_Value(double u[], vector<Texture*> textures){
	
	int col[3]={0,0,0};
	
	if(this->bump_ID>-1 && this->bump_ID<textures.size()){
		int w=textures[bump_ID]->get_width()-1;
		int h=textures[bump_ID]->get_height()-1;
		int i=int(u[0]*w);int j=int((1-u[1])*h);
		textures[bump_ID]->get_at(i,j,col);
	}
	else{
		col[0]=color[0];
		col[1]=color[1];
		col[2]=color[2];
	}
	return sqrt(double(col[0]*col[0]+col[1]*col[1]+col[2]*col[2]))/255.0;
}**/



CompoundShape::CompoundShape(int red,int green,int blue)
	: Shape(red,green,blue)
{

}

//4th component of x records which part is hit
double CompoundShape::get_Intersection(double p[],double v[],double x[]){
	
	double t=-1;
	double t_i;
	double y[4];
	for(int i=0;i<parts.size();i++){
		
		t_i=parts[i]->get_Intersection(p,v,x);
		if(t_i>0.01 && (t_i<t || t<0) ){
			t=t_i;
			y[0]=x[0];y[1]=x[1];y[2]=x[2];y[3]=double(i);
		}
	}
	if(t!=-1){ x[0]=y[0];x[1]=y[1];x[2]=y[2];x[3]=y[3];}
	return t;
}

double CompoundShape::get_Intersection2(double p[],double v[],double R,double x[]){
	double t=-1;
	double t_i;
	double y[4];
	for(int i=0;i<parts.size();i++){
		t_i=parts[i]->get_Intersection2(p,v,R,x);
		//cout<<"rat "<<t_i<<endl;
		if(t_i>0.01 && (t_i<t || t<0) ){
			t=t_i;
			y[0]=x[0];y[1]=x[1];y[2]=x[2];y[3]=double(i);
		}
	}
	if(t!=-1){ x[0]=y[0];x[1]=y[1];x[2]=y[2];x[3]=y[3];}
	return t;
}

//4th component of x tells which part is relevant
void CompoundShape::get_Texture_Coords(double x[], double u[]){
	parts[int(x[3])]->get_Texture_Coords(x,u);
}

//4th component of x tells which part is relevant
void CompoundShape::get_Normal(double x[], double n[]){
	//std::cout<<" part number "<<int(x[3])<<std::endl;
	parts[int(x[3])]->get_Normal(x,n);
}

void CompoundShape::set_Texture_IDs(int IDs[]){
	for(int i=0;i<parts.size();i++){
		parts[i]->texture_ID=IDs[i];
	}
}

void CompoundShape::get_Color(double x[], vector<Texture*> textures, int col[]){
	parts[int(x[3])]->get_Color(x,textures,col);
}

double CompoundShape::get_Alpha(double x[], vector<Texture*> textures){
	return parts[int(x[3])]->get_Alpha(x,textures);
}
