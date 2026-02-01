#include <math.h>
#include <iostream>
#include "Sphere.h"
//#include "matrix.h"
#include "poly_root.h"

//#include "quickhull.h"

#include <array>
#include <iterator>
#include <limits>


#define PI 3.14159

using namespace std;
//using namespace math;

Sphere::Sphere(double c1,double c2,double c3,double R,int red,int green, int blue)
	: BasicShape(red,green,blue)
{
   r=R;
   c[0]=c1;
   c[1]=c2;
   c[2]=c3;
   is_wireframe=false;
}

bool Sphere::wireframe_condition(double p[]){
	return false;
}

//p and v come from the vector equation of a line
//x(t)=p+vt.  returns the value of t at which the line
//intersects the sphere.  If their are two such t, it
//returns the one closest to p.  We are only interested in 
//things in front of the eye, so if t is negative it simply
//returns -1.  If the line does not hit the sphere, it also
//returns -1.  On exit, p and v are unchanged but x[]
//holds x(t_intersect) (unless no intersection is found),
//in which case it is unchanged

double Sphere::get_Intersection(double p[],double v[],double x[]){

   double d[3];d[0]=p[0]-c[0];d[1]=p[1]-c[1];d[2]=p[2]-c[2];
   double a=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
   double b=2*(v[0]*d[0]+v[1]*d[1]+v[2]*d[2]);
   double c=d[0]*d[0]+d[1]*d[1]+d[2]*d[2]-r*r;
   double desc=b*b-4*a*c;
   if(desc<0){
	return -1;
   }
   desc=sqrt(desc);
   double t1=(-b+desc)/(2*a);
   double t2=(-b-desc)/(2*a);
   if(t1>0 && t1<t2){ x[0]=p[0]+t1*v[0]; x[1]=p[1]+t1*v[1]; x[2]=p[2]+t1*v[2]; return t1;}
   else if(t2>0 && t2<=t1){x[0]=p[0]+t2*v[0]; x[1]=p[1]+t2*v[1]; x[2]=p[2]+t2*v[2];return t2;}
   return -1;
}

double Sphere::get_Intersection2(double P[],double V[],double R,double X[]){
	return 0;
	double p[3],v[3],x[3];
	contract(P,R,p);
	contract_derivative(P,V,R,v);

	double alpha=R*R+c[0]*c[0]+c[1]*c[1]+c[2]*c[2]-r*r;
	double A[3]={p[0]*p[0]+p[1]*p[1]+p[2]*p[2],2*(p[0]*v[0]+p[1]*v[1]+p[2]*v[2]),v[0]*v[0]+v[1]*v[1]+v[2]*v[2]};
	double B[2]={c[0]*p[0]+c[1]*p[1]+c[2]*p[2],c[0]*v[0]+c[1]*v[1]+c[2]*v[2]};
	double C[5]={A[0]*A[0],2*A[0]*A[1],A[1]*A[1]+2*A[0]*A[2],2*A[1]*A[2],A[2]*A[2]};
	double D[7]={A[0]*C[0],A[0]*C[1]+A[1]*C[0],A[0]*C[2]+A[1]*C[1]+A[2]*C[0],A[0]*C[3]+A[1]*C[2]+A[2]*C[1],A[0]*C[4]+A[1]*C[3]+A[2]*C[2],A[1]*C[4]+A[2]*C[3],A[2]*C[4]};
	double E[5]={A[0]*B[0]*B[0],B[0]*B[0]*A[1]+2*B[0]*B[1]*A[0],B[0]*B[0]*A[2]+2*B[0]*B[1]*A[1]+A[0]*B[1]*B[1],B[1]*B[1]*A[1]+2*A[2]*B[0]*B[1],A[2]*B[1]*B[1]};
	
	double coef[6],roots1[6],roots2[6];
	coef[6]=D[0]+2*alpha*C[0]+alpha*alpha*A[0]-4*E[0]-4*R*R*B[0]*B[0];
	coef[5]=D[1]+2*alpha*C[1]+alpha*alpha*A[1]-4*E[1]-8*R*R*B[0]*B[1];
	coef[4]=D[2]+2*alpha*C[2]+alpha*alpha*A[2]-4*E[2]-4*R*R*B[1]*B[1];
	coef[3]=D[3]+2*alpha*C[3]-4*E[3];
	coef[2]=D[4]+2*alpha*C[4]-4*E[4];
	coef[1]=D[5];
	coef[0]=D[6];

	rpoly(coef,6,roots1,roots2);

	//cout<<"roots1 "<<roots1[0]<<" "<<roots1[1]<<" "<<roots1[2]<<" "<<roots1[3]<<" "<<roots1[4]<<" "<<roots1[5]<<endl;
	//cout<<"roots2 "<<roots2[0]<<" "<<roots2[1]<<" "<<roots2[2]<<" "<<roots2[3]<<" "<<roots2[4]<<" "<<roots2[5]<<endl;

	double t=-1;
	double dist;
	for(int i=0;i<6;i++){
		x[0]=p[0]+roots1[i]*v[0];
		x[1]=p[1]+roots1[i]*v[1];
		x[2]=p[2]+roots1[i]*v[2];
		expand(x,R,X);
		dist=sqrt((X[0]-c[0])*(X[0]-c[0])+(X[1]-c[1])*(X[1]-c[1])+(X[2]-c[2])*(X[2]-c[2]))-r;
		if(dist<0.01){
			if( roots1[i]>0.00001 && (t<0 || roots1[i]<t) ){
				t=roots1[i];
			}
		}
		//else{cout<<"hmm "<<dist<<endl;}		
	}
	if(t!=-1){
	x[0]=p[0]+t*v[0];
	x[1]=p[1]+t*v[1];
	x[2]=p[2]+t*v[2];
	expand(x,R,X);
	//expand_derivative(x,v,R,V);
	}
	return t;

}

void Sphere::get_Texture_Coords(double x[], double u[]){
	double X[3];X[0]=x[0]-c[0];X[1]=x[1]-c[1];X[2]=x[2]-c[2];
	u[0]=(atan2(X[1],X[0])+PI)/(2*PI);
	u[1]=1-acos(X[2]/r)/PI;
	if(u[0]<0){u[0]=0;}
	if(u[1]<0){u[1]=0;}
	if(u[0]>1){u[0]=1;}
	if(u[1]>1){u[1]=1;}
}

void Sphere::get_Normal(double x[], double n[]){
	for(int i=0;i<3;i++){	
		n[i]=(x[i]-c[i])/r;
	}
}



InfinitePlane::InfinitePlane(double Z, int red,int green,int blue)
	: BasicShape(red,green,blue)
{
	z=Z;
	tiling=true;
	coordinate = 2;
	normal_dir=1;
	tile_size=6;
	upper_bound_coordinate=-1;
	upper_bound=-1;
	texture_center[0]=0;
	texture_center[1]=0;
	upper_bound_coefficient = 0;
	upper_bound_variable = 0;
}

double InfinitePlane::get_Intersection(double p[],double v[],double x[]){
	double t;
	double n[3]={0,0,0};
	this->get_Normal(x,n);
	double v_dot_n = v[0]*n[0]+v[1]*n[1]+v[2]*n[2];
	if(one_sided && v_dot_n <= 0){t=-1; return t;}
	if(v[coordinate]!=0){
		t = (z-p[coordinate])/v[coordinate];
		x[0]=p[0]+t*v[0];
		x[1]=p[1]+t*v[1];
		x[2]=p[2]+t*v[2];
		if(upper_bound_coordinate>-1 && x[upper_bound_coordinate]>upper_bound+upper_bound_coefficient*x[upper_bound_variable]){
			t=-1;
		}
	}
	else{
		t=-1;
	}
	return t;
}

double InfinitePlane::get_Intersection2(double P[],double V[],double R,double X[]){
	return 0;
}
void InfinitePlane::get_Texture_Coords(double x[], double u[]){
	
	int i1,i2;
	if(coordinate == 2){
		i1=0;
		i2=1;
	}
	else if(coordinate == 1){
		i1=0;
		i2=2;
	}
	else if(coordinate == 0){
		i1=1;
		i2=2;
	}
	
	if(tiling){
	
		u[0]=(x[i1]+texture_center[i1]-(-0.5*tile_size))/tile_size;
		u[1]=(x[i2]+texture_center[i2]-(-0.5*tile_size))/tile_size;
		//u[0]=u[0]-floor(u[0]/1)*1;
		//u[1]=u[1]-floor(u[1]/1)*1;
		u[0]=2*fabs(u[0]/2-floor(u[0]/2+0.5));
		u[1]=2*fabs(u[1]/2-floor(u[1]/2+0.5));
	
	}
	else{
		u[0]=(x[i2]+texture_center[i1]-(-0.5*tile_size))/tile_size;
		u[1]=(x[i2]+texture_center[i2]-(-0.5*tile_size))/tile_size;
		u[0]=fmax(fmin(u[0],1.0),0.0);
		u[1]=fmax(fmin(u[1],1.0),0.0);
	}
	
}

void InfinitePlane::get_Shading(double x[],vector<Texture*> textures,int obj_id,int C[]){
	int i1,i2;
	if(coordinate == 2){
		i1=0;
		i2=1;
	}
	else if(coordinate == 1){
		i1=0;
		i2=2;
	}
	else if(coordinate == 0){
		i1=1;
		i2=2;
	}
	
	double tile_size2 = 8.0;
	double u[2]={0,0};

	u[0]=(x[i1]-(-0.5*tile_size2))/tile_size2;
	u[1]=(x[i2]-(-0.5*tile_size2))/tile_size2;
	u[0]=fmax(fmin(u[0],1.0),0.0);
	u[1]=fmax(fmin(u[1],1.0),0.0);

	//cout<<"u = "<<u[0]<<" "<<u[1]<<endl;

	textures[obj_id]->bilinear_get_at(u,C);
	//if(coordinate==0 && x[1]<0){
	//	C[2]=0;
	//}
	
}

void InfinitePlane::inverse_texture(int i,int j,int n, int m,float query_pt[]){
	
	float x = -0.5*tile_size+i*tile_size/double(n);
	float y = -0.5*tile_size+j*tile_size/double(m);
	
	if(coordinate == 2){
		query_pt[0]=x;
		query_pt[1]=y;
		query_pt[2]=z;
	}
	else if(coordinate == 1){
		query_pt[0]=x;
		query_pt[1]=z;
		query_pt[2]=y;
	}
	else if(coordinate == 0){
		query_pt[0]=z;
		query_pt[1]=x;
		query_pt[2]=y;
	}
	return;
	
}

void InfinitePlane::place_dot(double x[],int ***im_data,int w,int h){
	
	
	
	if(coordinate == 2){
		int i = int(min(max(0,int(floor((x[0]+0.5*tile_size)*double(w)/tile_size))),w-1));
		int j = int(min(max(0,int(floor((x[1]+0.5*tile_size)*double(h)/tile_size))),h-1));
		im_data[i][j][0]=255;
		im_data[i][j][1]=255;
		im_data[i][j][2]=255;
	}
	else if(coordinate == 1){
		int i = int(min(max(0,int(floor((x[0]+0.5*tile_size)*double(w)/tile_size))),w-1));
		int j = int(min(max(0,int(floor((x[2]+0.5*tile_size)*double(h)/tile_size))),h-1));
		im_data[i][j][0]=255;
		im_data[i][j][1]=255;
		im_data[i][j][2]=255;
	}
	else if(coordinate == 0){
		int i = int(min(max(0,int(floor((x[1]+0.5*tile_size)*double(w)/tile_size))),w-1));
		int j = int(min(max(0,int(floor((x[2]+0.5*tile_size)*double(h)/tile_size))),h-1));
		im_data[i][j][0]=255;
		im_data[i][j][1]=255;
		im_data[i][j][2]=255;
	}
}

void InfinitePlane::get_Bump_Perturbation(double x[], vector<Texture*> textures, double dn[]){


	//cout<<"inside bump perturbation"<<endl;
	int col[3]={0,0,-2};
	double u[2]={0,0};
	this->get_Texture_Coords(x,u);

	if(this->bump_ID>-1 && this->bump_ID<textures.size()){
			//int w=textures[bump_ID]->get_width()-1;
			//int h=textures[bump_ID]->get_height()-1;
			//int i=int(u[0]*w);int j=int((1-u[1])*h);
			//textures[bump_ID]->get_at(i,j,col);
			//cout<<"u[0] = "<<u[0]<<" u[1] = "<<u[1]<<endl;
			textures[bump_ID]->bilinear_get_at(u, col);
	}
	else{
		col[0]=color[0];
		col[1]=color[1];
		col[2]=color[2];
	}
	//cout<<"col[0] = "<<col[0]<<" col[1] = "<<col[1]<<endl;
	dn[0]=double((col[1]-127)/127.0);
	dn[1]=-double((col[0]-127)/127.0);
	dn[2]=0;
}

bool InfinitePlane::wireframe_condition(double p[]){
	return false;
}

void InfinitePlane::get_Normal(double x[], double n[]){
n[0]=0;n[1]=0;n[2]=0;
n[coordinate]=normal_dir;
}

PlusSign3D::PlusSign3D(double c1,double c2,double c3,double Width,int red,int green, int blue)
	: BasicShape(red,green,blue)
{
   width=Width;
   c[0]=c1;
   c[1]=c2;
   c[2]=c3;
   is_wireframe=false;
}

bool PlusSign3D::wireframe_condition(double p[]){
	
	double x=p[0], y=p[1], z=p[2];
	


	x=x-c[0];
	y=y-c[1];
	z=z-c[2];
	
	int edge_variable = 0;
	int bound_variable = 0;

	double delta = width/5.0;

	if(abs(abs(x)-0.5*width)<delta){
		edge_variable++;
	}
	if(abs(abs(y)-0.5*width)<delta){
		edge_variable++;
	}
	if(abs(abs(z)-0.5*width)<delta){
		edge_variable++;
	}
	if(edge_variable>=1){
		return true;
	}
	else{
		return false;
	}
}

double PlusSign3D::get_Intersection(double p[],double v[],double x[]){

	double tp,tm;
	double X,Y,Z;
	double t = -1;

	for(int k=0;k<3;k++){
		if(v[k]!=0){
			tp = (-p[k]+(0.5*width+c[k]))/v[k];
			tm = (-p[k]+(-0.5*width+c[k]))/v[k];
			X = p[0]+tp*v[0];
			Y = p[1]+tp*v[1];
			Z = p[2]+tp*v[2];
			


			X=X-c[0];
			Y=Y-c[1];
			Z=Z-c[2];
			
			if( abs(X)<1.5*width && abs(Y)<1.5*width && abs(Z)<1.5*width){
				int count = 0;
				if(abs(X)>0.5*width){
					count++;
				}
				if(abs(Y)>0.5*width){
					count++;
				}
				if(abs(Z)>0.5*width){
					count++;
				}
				if(count<=1 && (t==-1 || tp<t) ){
					t = tp;
				}
			}
			X = p[0]+tm*v[0];
			Y = p[1]+tm*v[1];
			Z = p[2]+tm*v[2];
			


			X=X-c[0];
			Y=Y-c[1];
			Z=Z-c[2];
			
			if( abs(X)<1.5*width && abs(Y)<1.5*width && abs(Z)<1.5*width){
				int count = 0;
				if(abs(X)>0.5*width){
					count++;
				}
				if(abs(Y)>0.5*width){
					count++;
				}
				if(abs(Z)>0.5*width){
					count++;
				}
				if(count<=1 && (t==-1 || tm<t) ){
					t = tm;
				}
			}
		}
	}
	return t;

}

double PlusSign3D::get_Intersection2(double P[],double V[],double R,double X[]){
	return 0;
}

void PlusSign3D::get_Texture_Coords(double x[], double u[]){
	
	// not implemented
	u[0]=0;
	u[1]=0;
}

void PlusSign3D::get_Normal(double x[], double n[]){
	double xx=x[0], y=x[1], z=x[2];


	xx=xx-c[0];
	y=y-c[1];
	z=z-c[2];
	
	double dxp = abs(xx-width/2);
	double dxm = abs(xx+width/2);
	
	double dyp = abs(y-width/2);
	double dym = abs(y+width/2);
	double dzp = abs(z-width/2);
	double dzm = abs(z+width/2);
	
	double eps = 0.01;
	
	double d_min = min(dxp,min(dxm,min(dyp,min(dym,min(dzp,dzm)))));
	
	if( abs(d_min-dxp)<eps ){
		n[0]=1;n[1]=0;n[2]=0;
	}
	else if( abs(d_min-dxm)<eps){
		n[0]=-1;n[1]=0;n[2]=0;
	}
	else if( abs(d_min-dyp)<eps ){
		n[0]=0;n[1]=1;n[2]=0;
	}
	else if( abs(d_min-dym)<eps ){
		n[0]=0;n[1]=-1;n[2]=0;
	}
	else if( abs(d_min-dzp)<eps ){
		n[0]=0;n[1]=0;n[2]=1;
	}
	else if( abs(d_min-dzp)<eps ){
		n[0]=0;n[1]=0;n[2]=-1;
	}
	
}

Triangle::Triangle(double p1_[],double p2_[],double p3_[], double n_[], int red,int green,int blue)
	: BasicShape(red,green,blue)
{
	p1[0]=p1_[0];p1[1]=p1_[1];p1[2]=p1_[2];
	p2[0]=p2_[0];p2[1]=p2_[1];p2[2]=p2_[2];
	p3[0]=p3_[0];p3[1]=p3_[1];p3[2]=p3_[2];
	n[0]=n_[0];n[1]=n_[1];n[2]=n_[2];
	double len=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	n[0]/=len;n[1]/=len;n[2]/=len;
	//n[0]*=-1;n[1]*=-1;n[2]*=-1;
	is_wireframe=false;
	rat = 0;
}

bool Triangle::wireframe_condition(double p[]){
	return false;
}

void Triangle::scaled_translation(double x,double y, double z,double s){
	
	p1[0]=(p1[0]+x)/s;p2[0]=(p2[0]+x)/s;p3[0]=(p3[0]+x)/s;
	p1[1]=(p1[1]+y)/s;p2[1]=(p2[1]+y)/s;p3[1]=(p3[1]+y)/s;
	p1[2]=(p1[2]+z)/s;p2[2]=(p2[2]+z)/s;p3[2]=(p3[2]+z)/s;

}

double Triangle::get_Intersection(double p[],double v[],double x[]){
	
	// We want to solve Ax=b, where A = [p1-p3,p2-p3,-v], x=[u;vv;t], and b=p-p3
	// We're going to do this by brute force, using the general formula for the inverse of a 3x3 matrix.
	// First, define A =[a,b,c;d,e,f;g,h,i]
	double a,b,c,d,e,f,h,i;
	a=p1[0]-p3[0];
	b=p2[0]-p3[0];
	c=-v[0];
	d=p1[1]-p3[1];
	e=p2[1]-p3[1];
	f=-v[1];
	g=p1[2]-p3[2];
	h=p2[2]-p3[2];
	i=-v[2];
	
	double determ = a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g);
	if(determ==0){
		return -1; //ray parallel to triangle
	}
	else{
		//A^-1 = [A,B,C;D,E,F;G,H,I]
		double A,B,C,D,E,F,G,H,I;
		A = (e*i-f*h)/determ;
		B = (c*h-b*i)/determ;
		C = (b*f-c*e)/determ;
		D = (f*g-d*i)/determ;
		E = (a*i-c*g)/determ;
		F = (c*d-a*f)/determ;
		G = (d*h-e*g)/determ;
		H = (b*g-a*h)/determ;
		I = (a*e-b*d)/determ;
		double u = A*(p[0]-p3[0])+B*(p[1]-p3[1])+C*(p[2]-p3[2]);
		double vv = D*(p[0]-p3[0])+E*(p[1]-p3[1])+F*(p[2]-p3[2]);
		double t = G*(p[0]-p3[0])+H*(p[1]-p3[1])+I*(p[2]-p3[2]);
		double w = 1-u-vv; // Barycentric Coordinates (u,vv,w)


		if(!(0<=u && u<=1 && 0<=vv && vv<=1 && 0<=w && w<=1)){
			// we hit the plane, but missed the triangle
			return -1;
		}

		x[0]=p[0]+t*v[0];x[1]=p[1]+t*v[1];x[2]=p[2]+t*v[2];
		return t;
	}
}

double Triangle::get_Intersection2(double P[],double V[],double R,double X[]){

	// not implemented
	return -1;
   
}

void Triangle::get_Texture_Coords(double x[], double u[]){
	
	// not implemented
	u[0]=0;
	u[1]=0;
}

void Triangle::get_Normal(double x[], double N[]){
	N[0]=n[0];N[1]=n[1];N[2]=n[2];
}


Parallelogram::Parallelogram(double O[3],double a[3],double b[3],int red,int green, int blue)
	: BasicShape(red,green,blue)
{
	o[0]=O[0];o[1]=O[1];o[2]=O[2];
	x1[0]=a[0];x1[1]=a[1];x1[2]=a[2];
	x2[0]=b[0];x2[1]=b[1];x2[2]=b[2];
	c[0]=o[0]+0.5*(a[0]+b[0]);
	c[1]=o[1]+0.5*(a[1]+b[1]);
	c[2]=o[2]+0.5*(a[2]+b[2]);
	n[0]=x1[1]*x2[2]-x2[1]*x1[2];
	n[1]=x1[2]*x2[0]-x1[0]*x2[2];
	n[2]=x1[0]*x2[1]-x2[0]*x1[1];
	double len=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	n[0]/=len;n[1]/=len;n[2]/=len;
	n[0]*=-1;n[1]*=-1;n[2]*=-1;
	is_wireframe=false;
	rat = 0;
}

bool Parallelogram::wireframe_condition(double p[]){
	double u[2]={0,0};
	this->get_Texture_Coords(p, u);
	double delta = 0.02*2.0/3;
	if(rat==0){
		if( abs(u[0])<delta || abs(1-u[0])<delta || abs(u[1])<delta || abs(1-u[1])<delta){
			return true;
		}
		else{
			return false;
		}
	}
	else if(rat==1){
		if( abs(u[0])<delta || abs(u[1])<delta ){
			return true;
		}
		else{
			return false;
		}
	}
	else if(rat==2){
		if( abs(u[0]-1)<delta || abs(u[1]-1)<delta ){
			return true;
		}
		else{
			return false;
		}
	}
	else{
		return false;
	}
}

double Parallelogram::get_Intersection(double p[],double v[],double x[]){

	//have o+u*x1+vv*x2=p+v*t - 3x3 linear system for u,vv,t.
	/**matrix<double> A(3,3),b(3,1),X(3,1);

	for(int i=0;i<3;i++){
	A(i,0)=x1[i];
	A(i,1)=x2[i];
	A(i,2)=-v[i];
	b(i,0)=p[i]-o[i];
	}

	//X=[u v t]
	try {
		X=A.Solve(b);
	}
	catch (...){return -1;} //ray parallel to parallelogram**/

	// We want to solve Ax=b, where A = [x1,x2,-v], x=[u;vv;t], and b=p-o
	// We're going to do this by brute force, using the general formula for the inverse of a 3x3 matrix.
	// First, define A =[a,b,c;d,e,f;g,h,i]
	double a,b,c,d,e,f,h,i;
	a=x1[0];
	b=x2[0];
	c=-v[0];
	d=x1[1];
	e=x2[1];
	f=-v[1];
	g=x1[2];
	h=x2[2];
	i=-v[2];
	double determ = a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g);
	if(determ==0){
		return -1; //ray parallel to parallelogram
	}
	else{
		//A^-1 = [A,B,C;D,E,F;G,H,I]
		double A,B,C,D,E,F,G,H,I;
		A = (e*i-f*h)/determ;
		B = (c*h-b*i)/determ;
		C = (b*f-c*e)/determ;
		D = (f*g-d*i)/determ;
		E = (a*i-c*g)/determ;
		F = (c*d-a*f)/determ;
		G = (d*h-e*g)/determ;
		H = (b*g-a*h)/determ;
		I = (a*e-b*d)/determ;
		double u = A*(p[0]-o[0])+B*(p[1]-o[1])+C*(p[2]-o[2]);
		double vv = D*(p[0]-o[0])+E*(p[1]-o[1])+F*(p[2]-o[2]);
		double t = G*(p[0]-o[0])+H*(p[1]-o[1])+I*(p[2]-o[2]);


		if(!(0<=u && u<=1 && 0<=vv && vv<=1)){
			// we hit the plane, but missed the parallelogram
			return -1;
		}

		x[0]=p[0]+t*v[0];x[1]=p[1]+t*v[1];x[2]=p[2]+t*v[2];
		return t;
	}
}

//intersection routine for the 1st invisibility model
double Parallelogram::get_Intersection2(double P[],double Vo[],double R,double X[]){
	return 0;
	double coef[5],roots1[4],roots2[4];
	double p[3],v[3],x[3];

	contract(P,R,p);
	contract_derivative(P,Vo,R,v);

	//first get the coefficients of the relevant 
	//4th degree polynomial

	double c=n[0]*o[0]+n[1]*o[1]+n[2]*o[2];
	double pdotp,vdotp,vdotv,pdotn,vdotn;
	pdotp=p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
	vdotp=v[0]*p[0]+v[1]*p[1]+v[2]*p[2];
	vdotv=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
	pdotn=p[0]*n[0]+p[1]*n[1]+p[2]*n[2];
	vdotn=v[0]*n[0]+v[1]*n[1]+v[2]*n[2];

	double alpha,beta,gamma,U,V,W;

	alpha=pdotp; beta=2*vdotp ; gamma=vdotv;
	U=pdotn*pdotn; V=2*vdotn*pdotn; W=vdotn*vdotn;

	double Alpha=alpha+R*R;

	coef[0]=gamma*W;
	coef[1]=gamma*V+W*beta;
	coef[2]=Alpha*W+beta*V+gamma*U-c*c*gamma;
	coef[3]=Alpha*V+beta*U-c*c*beta;
	coef[4]=Alpha*U-c*c*alpha;

	//find the closest root
	rpoly(coef,4,roots1,roots2);

	double t=-1;
	double dist=-1,cand;
	for(int i=0;i<4;i++){
		if(fabs(roots2[i])<1){
			x[0]=p[0]+roots1[i]*v[0];
			x[1]=p[1]+roots1[i]*v[1];
			x[2]=p[2]+roots1[i]*v[2];
			expand(x,R,X);
			cand=fabs((X[0]-o[0])*n[0]+(X[1]-o[1])*n[1]+(X[2]-o[2])*n[2]);		
			if( cand<dist || dist<0){
				dist=cand;
				t=roots1[i];
			}
		}
	}

	if(t!=-1){
		x[0]=p[0]+t*v[0];
		x[1]=p[1]+t*v[1];
		x[2]=p[2]+t*v[2];
		expand(x,R,X);
		double tex[2];
		this->cheap_hack(X,tex);
		if(!(tex[0]>=0 && tex[0]<=1 && tex[1]>=0 && tex[1]<=1)){return -1;}
	}
	//change tangent vector
	expand_derivative(p,v,R,Vo);
	//move in normal direction to avoid spotting
	//X[0]=X[0]-0.01*n[0];
	//X[1]=X[1]-0.01*n[1];
	//X[2]=X[2]-0.01*n[2];
	return t;	
}


void Parallelogram::cheap_hack(double x[], double u[]){
//minimize ||u*x1+v*x2+o-x||^2 by setting (u,v) gradient = 0.
//get 2x2 system

/**matrix<double> A(2,2),X(2,1),b(2,1);

A(0,0)=x1[0]*x1[0]+x1[1]*x1[1]+x1[2]*x1[2];
A(0,1)=x1[0]*x2[0]+x1[1]*x2[1]+x1[2]*x2[2];
A(1,0)=A(0,1);
A(1,1)=x2[0]*x2[0]+x2[1]*x2[1]+x2[2]*x2[2];
b(0,0)=x1[0]*(x[0]-o[0])+x1[1]*(x[1]-o[1])+x1[2]*(x[2]-o[2]);
b(1,0)=x2[0]*(x[0]-o[0])+x2[1]*(x[1]-o[1])+x2[2]*(x[2]-o[2]);**/
double a,b,c,d;
a = x1[0]*x1[0]+x1[1]*x1[1]+x1[2]*x1[2];
b = x1[0]*x2[0]+x1[1]*x2[1]+x1[2]*x2[2];
c = b;
d = x2[0]*x2[0]+x2[1]*x2[1]+x2[2]*x2[2];
double b1,b2;
b1 = x1[0]*(x[0]-o[0])+x1[1]*(x[1]-o[1])+x1[2]*(x[2]-o[2]);
b2 = x2[0]*(x[0]-o[0])+x2[1]*(x[1]-o[1])+x2[2]*(x[2]-o[2]);

//X=[u,v]
//X=A.Solve(b);
//u[0]=X(0,0);
//u[1]=X(1,0);
double det = (a*d-b*c);
u[0] = (d*b1-b*b2)/det;
u[1] = (a*b1-c*b2)/det;



}


void Parallelogram::get_Texture_Coords(double x[], double u[]){
//minimize ||u*x1+v*x2+o-x||^2 by setting (u,v) gradient = 0.
//get 2x2 system
/**matrix<double> A(2,2),X(2,1),b(2,1);

A(0,0)=x1[0]*x1[0]+x1[1]*x1[1]+x1[2]*x1[2];
A(0,1)=x1[0]*x2[0]+x1[1]*x2[1]+x1[2]*x2[2];
A(1,0)=A(0,1);
A(1,1)=x2[0]*x2[0]+x2[1]*x2[1]+x2[2]*x2[2];
b(0,0)=x1[0]*(x[0]-o[0])+x1[1]*(x[1]-o[1])+x1[2]*(x[2]-o[2]);
b(1,0)=x2[0]*(x[0]-o[0])+x2[1]*(x[1]-o[1])+x2[2]*(x[2]-o[2]);

//X=[u,v]
X=A.Solve(b);
u[0]=X(0,0);
u[1]=X(1,0);**/

double a,b,c,d;
a = x1[0]*x1[0]+x1[1]*x1[1]+x1[2]*x1[2];
b = x1[0]*x2[0]+x1[1]*x2[1]+x1[2]*x2[2];
c = b;
d = x2[0]*x2[0]+x2[1]*x2[1]+x2[2]*x2[2];
double b1,b2;
b1 = x1[0]*(x[0]-o[0])+x1[1]*(x[1]-o[1])+x1[2]*(x[2]-o[2]);
b2 = x2[0]*(x[0]-o[0])+x2[1]*(x[1]-o[1])+x2[2]*(x[2]-o[2]);

double det = (a*d-b*c);
u[0] = (d*b1-b*b2)/det;
u[1] = (-c*b1+a*b2)/det;


if(u[0]<0){u[0]=0;}
if(u[1]<0){u[1]=0;}
if(u[0]>1){u[0]=1;}
if(u[1]>1){u[1]=1;}

}

void Parallelogram::get_Normal(double x[], double N[]){
	N[0]=n[0];N[1]=n[1];N[2]=n[2];
}


Upright_Cylinder::Upright_Cylinder(double zb,double zt,double x0, double y0, double r,int red,int green, int blue)
	: BasicShape(red,green,blue)
{
   R=r;
   z_b = zb;
   z_t = zt;
   x_0 = x0;
   y_0 = y0;
   
   is_wireframe=false;  
   
}

bool Upright_Cylinder::wireframe_condition(double p[]){
	return true;
}

double Upright_Cylinder::get_Intersection(double p[],double v[],double x[]){

double a = v[0]*v[0]+v[1]*v[1];
double b = 2*v[0]*(p[0]-x_0)+2*v[1]*(p[1]-y_0);
double c = (p[0]-x_0)*(p[0]-x_0)+(p[1]-y_0)*(p[1]-y_0)-R*R;
double t1=-1;
double t2=-1;
double desc = b*b-4*a*c;
if(fabs(a)>0 && desc>=0){
	double sqrt_desc = sqrt(desc);
	t1 = (-b+sqrt_desc)/(2*a);
	t2 = (-b-sqrt_desc)/(2*a);
	double z1 = p[2]+t1*v[2];
	double z2 = p[2]+t2*v[2];
	if(z1>z_t || z1<z_b){
		t1=-1;
	}
	if(z2>z_t || z2<z_b){
		t2=-1;
	}
}
double t12;
if(t1!=-1 && t2!=-1){
	t12 = fmin(t1,t2);
}
else if(t1==-1 && t2!=-1){
	t12=t2;
}
else if(t1!=-1 && t2==-1){
	t12=t1;
}
else{
	t12=-1;
}

double t3=-1;
double t4=-1;
if(v[2]!=0){
	t3 = (z_t-p[2])/v[2];
	if((p[0]+t3*v[0]-x_0)*(p[0]+t3*v[0]-x_0)+(p[1]+t3*v[1]-y_0)*(p[1]+t3*v[1]-y_0)>R*R){
		t3=-1;
	}
	t4 = (z_b-p[2])/v[2];
	if((p[0]+t4*v[0]-x_0)*(p[0]+t4*v[0]-x_0)+(p[1]+t4*v[1]-y_0)*(p[1]+t4*v[1]-y_0)>R*R){
		t4=-1;
	}
}
double t34;
if(t3!=-1 && t4!=-1){
	t34 = fmin(t3,t4);
}
else if(t3==-1 && t4!=-1){
	t34=t4;
}
else if(t3!=-1 && t4==-1){
	t34=t3;
}
else{
	t34=-1;
}
double t;
if(t12!=-1 && t34!=-1){
	t = fmin(t12,t34);
}
else if(t12==-1 && t34!=-1){
	t=t34;
}
else if(t12!=-1 && t34==-1){
	t=t12;
}
else{
	t=-1;
}

return t;

}

double Upright_Cylinder::get_Intersection2(double P[],double V[],double R,double X[]){

	// not implemented
	return -1;
   
}

void Upright_Cylinder::get_Texture_Coords(double x[], double u[]){
	
	// not implemented
	u[0]=0;
	u[1]=0;
}

void Upright_Cylinder::get_Normal(double x[], double n[]){
	double d_center =  fabs(sqrt( (x[0]-x_0)*(x[0]-x_0)+(x[1]-y_0)*(x[1]-y_0) )-R);
	double d_top = fabs(x[2]-z_t);
	double d_bottom = fabs(x[2]-z_b);
	if(min(d_top,d_bottom)<d_center){
		if(d_top<d_bottom){n[0]=0;n[1]=0;n[2]=1;}
		else{n[0]=0;n[1]=0;n[2]=-1;}
	}
	else{
		n[0]=(x[0]-x_0)/R;
		n[1]=(x[1]-y_0)/R;
		n[2]=0;
	}
}









Cylinder::Cylinder(double P1[3],double P2[3],double R,int red,int green, int blue)
	: BasicShape(red,green,blue)
{
   r=R;
   p1[0]=P1[0];p1[1]=P1[1];p1[2]=P1[2];
   p2[0]=P2[0];p2[1]=P2[1];p2[2]=P2[2];
   is_wireframe=false;
   
   double L = sqrt( (p1[0]-p2[0])*(p1[0]-p2[0])+(p1[1]-p2[1])*(p1[1]-p2[1])+(p1[2]-p2[2])*(p1[2]-p2[2]) );
   H = L;
   e1[0]=(p2[0]-p1[0])/L;
   e1[1]=(p2[1]-p1[1])/L;
   e1[2]=(p2[2]-p1[2])/L;
   if(!(e1[0]==0 && e1[1]==0)){
	   e2[0]=-e1[1];
	   e2[1]=e1[0];
	   e2[2]=0;
   }
   else if(!(e1[0]==0 && e1[2]==0)){
	   e2[0]=-e1[2];
	   e2[2]=e1[0];
	   e2[1]=0;
   }
   L = sqrt( e2[0]*e2[0]+e2[1]*e2[1]+e2[2]*e2[2] );
   e2[0]=e2[0]/L;
   e2[1]=e2[1]/L;
   e2[2]=e2[2]/L;
   
   e3[0]=e1[1]*e2[2]-e2[1]*e1[2];
   e3[1]=e1[2]*e2[0]-e1[0]*e2[2];
   e3[2]=e1[0]*e2[1]-e2[0]*e1[1];
   cout<<e1[0]<<" "<<e1[1]<<" "<<e1[2]<<endl;
   cout<<e2[0]<<" "<<e2[1]<<" "<<e2[2]<<endl;
   cout<<e3[0]<<" "<<e3[1]<<" "<<e3[2]<<endl;
   
   
}

bool Cylinder::wireframe_condition(double p[]){
	return true;
}

double Cylinder::get_Intersection(double p[],double v[],double x[]){

	double v_dot_e2=v[0]*e2[0]+v[1]*e2[1]+v[2]*e2[2];
	double v_dot_e3=v[0]*e3[0]+v[1]*e3[1]+v[2]*e3[2];
	double dp_dot_e2 = (p[0]-p1[0])*e2[0]+(p[1]-p1[1])*e2[1]+(p[2]-p1[2])*e2[2];
	double dp_dot_e3 = (p[0]-p1[0])*e3[0]+(p[1]-p1[1])*e3[1]+(p[2]-p1[2])*e3[2];

	double a=v_dot_e2*v_dot_e2+v_dot_e3*v_dot_e3;
   double b=2*(dp_dot_e2*v_dot_e2+dp_dot_e3*v_dot_e3);
   double c=dp_dot_e2*dp_dot_e2+dp_dot_e3*dp_dot_e3-r*r;
   double desc=b*b-4*a*c;
   if(desc<=0){
	return -1;
   }
   desc=sqrt(desc);
   double t1,t2;;
   if(abs(a)<0.1){
	   if(abs(-b-desc)<0.1){t1=-1;}else{t1=2*c/(-b-desc);}
	   if(abs(-b+desc)<0.1){t2=-1;}else{t2=2*c/(-b+desc);}
	}
   else{
	t1=(-b+desc)/(2*a);
	t2=(-b-desc)/(2*a);
   }
   double t_min;
   if(t1>0 && t1<t2){ x[0]=p[0]+t1*v[0]; x[1]=p[1]+t1*v[1]; x[2]=p[2]+t1*v[2]; t_min = t1;}
   else if(t2>0 && t2<=t1){x[0]=p[0]+t2*v[0]; x[1]=p[1]+t2*v[1]; x[2]=p[2]+t2*v[2]; t_min = t2;}
   else{t_min=-1;}
   
   double h = (x[0]-p1[0])*e1[0]+(x[1]-p1[1])*e1[1]+(x[2]-p1[2])*e1[2];
   if(h>=0 && h<=H){return t_min;}
   else{return -1;}
   return t_min;
   
}

double Cylinder::get_Intersection2(double P[],double V[],double R,double X[]){

	// not implemented
	return -1;
   
}

void Cylinder::get_Texture_Coords(double x[], double u[]){
	
	// not implemented
	u[0]=0;
	u[1]=0;
}

void Cylinder::get_Normal(double x[], double n[]){
	double dx_dot_e2 = (x[0]-p1[0])*e2[0]+(x[1]-p1[1])*e2[1]+(x[2]-p1[2])*e2[2];
	double dx_dot_e3 = (x[0]-p1[0])*e3[0]+(x[1]-p1[1])*e3[1]+(x[2]-p1[2])*e3[2];
	n[0]=dx_dot_e2*e2[0]+dx_dot_e3*e3[0];
	n[1]=dx_dot_e2*e2[1]+dx_dot_e3*e3[1];
	n[2]=dx_dot_e2*e2[2]+dx_dot_e3*e3[2];
	double L = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	n[0]=n[0]/L;n[1]=n[1]/L;n[1]=n[1]/L;
}




Cube::Cube(double c1, double c2, double c3, double L, double W, double H, int r, int g, int b)
	: CompoundShape(r,g,b)
{
	c[0]=c1;c[1]=c2;c[2]=c3;
	dim[0]=L;dim[1]=W;dim[2]=H;

	double o[3],x1[3],x2[3];
	//construct sides

	o[0]=c[0]+0.5*dim[0];o[1]=c[1]-0.5*dim[1];o[2]=c[2]-0.5*dim[2];
	x1[0]=0;x1[1]=dim[1];x1[2]=0;
	x2[0]=0;x2[1]=0;x2[2]=dim[2];
	Parallelogram *P1=new Parallelogram(o,x1,x2,r,g,b);
	parts.push_back(P1);

	o[0]=c[0]+0.5*dim[0];o[1]=c[1]+0.5*dim[1];o[2]=c[2]-0.5*dim[2];
	x1[0]=-dim[0];x1[1]=0;x1[2]=0;
	x2[0]=0;x2[1]=0;x2[2]=dim[2];
	Parallelogram *P2=new Parallelogram(o,x1,x2,r,g,b);
	parts.push_back(P2);

	o[0]=c[0]-0.5*dim[0];o[1]=c[1]+0.5*dim[1];o[2]=c[2]-0.5*dim[2];
	x1[0]=0;x1[1]=-dim[1];x1[2]=0;
	x2[0]=0;x2[1]=0;x2[2]=dim[2];
	Parallelogram *P3=new Parallelogram(o,x1,x2,r,g,b);
	parts.push_back(P3);

	o[0]=c[0]-0.5*dim[0];o[1]=c[1]-0.5*dim[1];o[2]=c[2]-0.5*dim[2];
	x1[0]=dim[0];x1[1]=0;x1[2]=0;
	x2[0]=0;x2[1]=0;x2[2]=dim[2];
	Parallelogram *P4=new Parallelogram(o,x1,x2,r,g,b);
	parts.push_back(P4);

	o[0]=c[0]+0.5*dim[0];o[1]=c[1]+0.5*dim[1];o[2]=c[2]+0.5*dim[2];
	x1[0]=-dim[0];x1[1]=0;x1[2]=0;
	x2[0]=0;x2[1]=-dim[1];x2[2]=0;
	Parallelogram *P5=new Parallelogram(o,x1,x2,r,g,b);
	parts.push_back(P5);

	o[0]=c[0]-0.5*dim[0];o[1]=c[1]-0.5*dim[1];o[2]=c[2]-0.5*dim[2];
	x1[0]=0;x1[1]=dim[1];x1[2]=0;
	x2[0]=dim[0];x2[1]=0;x2[2]=0;
	Parallelogram *P6=new Parallelogram(o,x1,x2,r,g,b);
	parts.push_back(P6);
	
	is_wireframe = false;
	
}

bool Cube::wireframe_condition(double p[]){
	return false;
}

/**ConvexHull::ConvexHull(std::vector<double> X,std::vector<double> Y,std::vector<double> Z,int r, int g, int b)
	: CompoundShape(r,g,b)
{
	int num_points = X.size();
	is_wireframe = false;
	
	using F = double;
	constexpr std::size_t dim = 3;
    using Points = std::vector<std::array<F, dim>>;
	Points points(num_points); 
	for(int i=0;i<num_points;i++){
		points.at(i)={X[i],Y[i],Z[i]};
	}
	
	const auto eps = std::numeric_limits<F>::epsilon();
    quick_hull<typename Points::const_iterator> qh{dim, eps};
    qh.add_points(std::cbegin(points), std::cend(points));
    auto initial_simplex = qh.get_affine_basis();
    if (initial_simplex.size() < dim + 1) {
        std::cout<<"degenerate input set"<<std::endl; // degenerated input set
    }
    qh.create_initial_simplex(std::cbegin(initial_simplex), std::prev(std::cend(initial_simplex)));
    qh.create_convex_hull();
    if (!qh.check()) {
        std::cout<<"not convex"<<std::endl; // resulted structure is not convex (generally due to precision errors)
    }

    qh.facets_; // use as result
	int count = 0;
	for (auto const & facet_ : qh.facets_) {
		std::cout<<"facet number "<<count<<std::endl;
		double p1[3];
		double p2[3];
		double p3[3];
		double n[3];
        auto const & vertices_ = facet_.vertices_;
		count = 0;
        for (auto const & vertex_ : vertices_) {
            for (F const & coordinate_ : *vertex_) {
				std::cout<<"count "<<count<<"count%3 "<<count%3<<std::endl;
				if( 0<=count && count<3){p1[count%3]=coordinate_;count++;std::cout<<" in case 1"<<std::endl;}
				else if ( 3<=count && count<6){p2[count%3]=coordinate_;count++;std::cout<<" in case 2"<<std::endl;}
				else{p3[count%3]=coordinate_;count++;std::cout<<" in case 3"<<std::endl;}
				
            }
		}
		std::cout<<p1[0]<<' '<<p1[1]<<' '<<p1[2]<<std::endl;
		std::cout<<p2[0]<<' '<<p2[1]<<' '<<p2[2]<<std::endl;
		std::cout<<p3[0]<<' '<<p3[1]<<' '<<p3[2]<<std::endl;
		count = 0;
		auto const  normal_ = facet_.normal_;
		for (F const & coordinate_ : normal_) {
			n[count]=coordinate_;
			count++;
		}
		std::cout<<n[0]<<' '<<n[1]<<' '<<n[2]<<std::endl;
		Triangle *T = new Triangle(p1,p2,p3,n,r,g,b);
		parts.push_back(T);
		std::cout<<" parts.size() "<<parts.size()<<std::endl;
	}
}
bool ConvexHull::wireframe_condition(double p[]){
	return false;
}**/


PlusSign3DMesh::PlusSign3DMesh(double c1,double c2,double c3,double Width,int r,int g, int b)
	: CompoundShape(r,g,b)
{
   width=Width;
   c[0]=c1;
   c[1]=c2;
   c[2]=c3;
   double R_s = 0.02;
   Sphere *S = new Sphere(c1+1.5*width,c2+1.5*width,c3+1.5*width,R_s,255,255,255);
   parts.push_back(S);
   is_wireframe=false;
   double R = 0.01;
   for(int i=0;i<2;i++){
	   for(int j=0;j<2;j++){
			double p1[3]={0.0,(1+i)*width+c[1],(1+j)*width+c[2]};
			double p2[3]={width+c[0],(1+i)*width+c[1],(1+j)*width+c[2]};
			Cylinder *C1 = new Cylinder(p1,p2,R,r,g,b);
			Sphere *S1a = new Sphere(p1[0],p1[1],p1[2],R,r,g,b);
			Sphere *S1b = new Sphere(p2[0],p2[1],p2[2],R,r,g,b);
			parts.push_back(C1);
			double p3[3]={1.0,(1+i)*width+c[1],(1+j)*width+c[2]};
			double p4[3]={2*width+c[0],(1+i)*width+c[1],(1+j)*width+c[2]};
			Cylinder *C2 = new Cylinder(p3,p4,R,r,g,b);
			Sphere *S2a = new Sphere(p3[0],p3[1],p3[2],R,r,g,b);
			Sphere *S2b = new Sphere(p4[0],p4[1],p4[2],R,r,g,b);
			parts.push_back(C2);
			parts.push_back(S1a);
			parts.push_back(S1b);
			parts.push_back(S2a);
			parts.push_back(S2b);
			double p5[3]={0.0,(1+i)*width+c[1],(1+j)*width+c[2]};
			if(i==0){
				double p6[3]={0.0,(2+i)*width+c[1],(1+j)*width+c[2]};
				Cylinder *C3 = new Cylinder(p5,p6,R,r,g,b);
				parts.push_back(C3);
			}
			if(j==0){
				double p6[3]={0.0,(1+i)*width+c[1],(2+j)*width+c[2]};
				Cylinder *C3 = new Cylinder(p5,p6,R,r,g,b);
				parts.push_back(C3);
			}
			double p7[3]={1.0,(1+i)*width+c[1],(1+j)*width+c[2]};
			if(i==0){
				double p8[3]={1.0,(2+i)*width+c[1],(1+j)*width+c[2]};
				Cylinder *C3 = new Cylinder(p7,p8,R,r,g,b);
				parts.push_back(C3);
			}
			if(j==0){
				double p8[3]={1.0,(1+i)*width+c[1],(2+j)*width+c[2]};
				Cylinder *C3 = new Cylinder(p7,p8,R,r,g,b);
				parts.push_back(C3);
			}
	   }
   }
   for(int i=0;i<2;i++){
	   for(int j=0;j<2;j++){
			double p1[3]={(1+i)*width+c[0],0.0,(1+j)*width+c[2]};
			double p2[3]={(1+i)*width+c[0],width+c[1],(1+j)*width+c[2]};
			Cylinder *C1 = new Cylinder(p1,p2,R,r,g,b);
			Sphere *S1a = new Sphere(p1[0],p1[1],p1[2],R,r,g,b);
			Sphere *S1b = new Sphere(p2[0],p2[1],p2[2],R,r,g,b);
			parts.push_back(C1);
			double p3[3]={(1+i)*width+c[0],1.0,(1+j)*width+c[2]};
			double p4[3]={(1+i)*width+c[0],2*width+c[1],(1+j)*width+c[2]};
			Cylinder *C2 = new Cylinder(p3,p4,R,r,g,b);
			Sphere *S2a = new Sphere(p3[0],p3[1],p3[2],R,r,g,b);
			Sphere *S2b = new Sphere(p4[0],p4[1],p4[2],R,r,g,b);
			parts.push_back(C2);
			parts.push_back(S1a);
			parts.push_back(S1b);
			parts.push_back(S2a);
			parts.push_back(S2b);
			double p5[3]={(1+i)*width+c[0],0.0,(1+j)*width+c[2]};
			if(i==0){
				double p6[3]={(2+i)*width+c[0],0.0,(1+j)*width+c[2]};
				Cylinder *C3 = new Cylinder(p5,p6,R,r,g,b);
				parts.push_back(C3);
			}
			if(j==0){
				double p6[3]={(1+i)*width+c[0],0.0,(2+j)*width+c[2]};
				Cylinder *C3 = new Cylinder(p5,p6,R,r,g,b);
				parts.push_back(C3);
			}
			double p7[3]={(1+i)*width+c[0],1.0,(1+j)*width+c[2]};
			if(i==0){
				double p8[3]={(2+i)*width+c[0],1.0,(1+j)*width+c[2]};
				Cylinder *C3 = new Cylinder(p7,p8,R,r,g,b);
				parts.push_back(C3);
			}
			if(j==0){
				double p8[3]={(1+i)*width+c[0],1.0,(2+j)*width+c[2]};
				Cylinder *C3 = new Cylinder(p7,p8,R,r,g,b);
				parts.push_back(C3);
			}
	   }
   }
   for(int i=0;i<2;i++){
	   for(int j=0;j<2;j++){
			double p1[3]={(1+i)*width+c[0],(1+j)*width+c[1],0.0};
			double p2[3]={(1+i)*width+c[0],(1+j)*width+c[1],width+c[2]};
			Cylinder *C1 = new Cylinder(p1,p2,R,r,g,b);
			Sphere *S1a = new Sphere(p1[0],p1[1],p1[2],R,r,g,b);
			Sphere *S1b = new Sphere(p2[0],p2[1],p2[2],R,r,g,b);
			parts.push_back(C1);
			double p3[3]={(1+i)*width+c[0],(1+j)*width+c[1],1.0};
			double p4[3]={(1+i)*width+c[0],(1+j)*width+c[1],2*width+c[2]};
			Cylinder *C2 = new Cylinder(p3,p4,R,r,g,b);
			Sphere *S2a = new Sphere(p3[0],p3[1],p3[2],R,r,g,b);
			Sphere *S2b = new Sphere(p4[0],p4[1],p4[2],R,r,g,b);
			parts.push_back(C2);
			parts.push_back(S1a);
			parts.push_back(S1b);
			parts.push_back(S2a);
			parts.push_back(S2b);
			double p5[3]={(1+i)*width+c[0],(1+j)*width+c[1],0.0};
			if(i==0){
				double p6[3]={(2+i)*width+c[0],(1+j)*width+c[1],0.0};
				Cylinder *C3 = new Cylinder(p5,p6,R,r,g,b);
				parts.push_back(C3);
			}
			if(j==0){
				double p6[3]={(1+i)*width+c[0],(2+j)*width+c[1],0.0};
				Cylinder *C3 = new Cylinder(p5,p6,R,r,g,b);
				parts.push_back(C3);
			}
			double p7[3]={(1+i)*width+c[0],(1+j)*width+c[1],1.0};
			if(i==0){
				double p8[3]={(2+i)*width+c[0],(1+j)*width+c[1],1.0};
				Cylinder *C3 = new Cylinder(p7,p8,R,r,g,b);
				parts.push_back(C3);
			}
			if(j==0){
				double p8[3]={(1+i)*width+c[0],(2+j)*width+c[1],1.0};
				Cylinder *C3 = new Cylinder(p7,p8,R,r,g,b);
				parts.push_back(C3);
			}
	   }
   }
   

}

bool PlusSign3DMesh::wireframe_condition(double p[]){
	return true;
}
