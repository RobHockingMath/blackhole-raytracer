#include <math.h>
#include <iostream>
#include "Distance_fractal.h"
#include "Sphere.h"
#include "invisibility.h"

#define PI 3.14159

using namespace std;

Distance_fractal::Distance_fractal(double c1,double c2,double c3,double E,int red,int green, int blue)
	: BasicShape(red,green,blue)
{
   c[0]=c1;c[1]=c2;c[2]=c3;
   epsilon=E;
   bounding_radius=10*sqrt(5);
   //bounding_sphere=new Sphere(0,0,0,sqrt(5.0),255,0,0);
}

double Distance_fractal::get_Intersection(double p[],double V[],double x[])
{
	double v[3];
	double d=sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
	for(int i=0;i<3;i++){v[i]=V[i]/d;}

	//transform to coord sys where set center is origin
	double P[3]={p[0]-c[0],p[1]-c[1],p[2]-c[2]};
	double dist=0;
	/**if(P[0]*P[0]+P[1]*P[1]+P[2]*P[2]>bounding_radius*bounding_radius){
	//define bounding sphere
	Sphere* S=new Sphere(0,0,0,bounding_radius,255,0,0);
	dist=S->get_Intersection(P,v,x);
	free(S);
	}
	else{x[0]=P[0];x[1]=P[1];x[2]=P[2];}**/
	x[0]=P[0];x[1]=P[1];x[2]=P[2];
if(dist!=-1){
	double delta=Distance(x);
	dist+=delta;
	while(delta>epsilon && x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<bounding_radius*bounding_radius+0.1){
		for(int i=0;i<3;i++){
			x[i]=P[i]+dist*v[i];
		}
		delta=Distance(x);
		dist+=delta;
		//cout<<"delta "<<delta<<endl;
	}
	x[0]+=c[0];x[1]+=c[1];x[2]+=c[2];
if(delta<=epsilon){return dist;}//{return dist/d;}
}
return -1;
}

double Distance_fractal::get_Intersection2(double P[],double v[],double R,double x[]){

//I hate this functin

return -1;

}

void Distance_fractal::get_Texture_Basis(double x[],double e1[],double e2[]){
	double n[3];
	get_Normal(x,n);
	if(!(n[0]==0 && n[1]==0)){
	e1[0]=-n[1];e1[1]=n[0];e1[2]=0;}
	else{e1[0]=1;e1[1]=0;e1[2]=0;}
	double len=sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);
	for(int i=0;i<3;i++){e1[i]/=len;}
	e2[0]=n[1]*e1[2]-n[2]*e1[1];
	e2[1]=n[2]*e1[0]-n[0]*e1[2];
	e2[2]=n[0]*e1[1]-n[1]*e1[0];
	
}

//use gradient of distance function
void Distance_fractal::get_Normal(double x[], double n[]){

	double x0[3]={x[0]-c[0],x[1]-c[1],x[2]-c[2]};
	double x1[3]={x[0]-c[0]+epsilon,x[1]-c[1],x[2]-c[2]};
	double x2[3]={x[0]-c[0],x[1]-c[1]+epsilon,x[2]-c[2]};
	double x3[3]={x[0]-c[0],x[1]-c[1],x[2]-c[2]+epsilon};
	double x4[3]={x[0]-c[0]-epsilon,x[1]-c[1],x[2]-c[2]};
	double x5[3]={x[0]-c[0],x[1]-c[1]-epsilon,x[2]-c[2]};
	double x6[3]={x[0]-c[0],x[1]-c[1],x[2]-c[2]-epsilon};
	double x7[3]={x[0]-c[0]-epsilon,x[1]-c[1]-epsilon,x[2]-c[2]};
	double x8[3]={x[0]-c[0]-epsilon,x[1]-c[1]+epsilon,x[2]-c[2]};
	double x9[3]={x[0]-c[0]+epsilon,x[1]-c[1]-epsilon,x[2]-c[2]};
	double x10[3]={x[0]-c[0]+epsilon,x[1]-c[1]+epsilon,x[2]-c[2]};
	double x11[3]={x[0]-c[0]-epsilon,x[1]-c[1],x[2]-c[2]-epsilon};
	double x12[3]={x[0]-c[0]-epsilon,x[1]-c[1],x[2]-c[2]+epsilon};
	double x13[3]={x[0]-c[0]+epsilon,x[1]-c[1],x[2]-c[2]-epsilon};
	double x14[3]={x[0]-c[0]+epsilon,x[1]-c[1],x[2]-c[2]+epsilon};
	double x15[3]={x[0]-c[0],x[1]-c[1]-epsilon,x[2]-c[2]-epsilon};
	double x16[3]={x[0]-c[0],x[1]-c[1]-epsilon,x[2]-c[2]+epsilon};
	double x17[3]={x[0]-c[0],x[1]-c[1]+epsilon,x[2]-c[2]-epsilon};
	double x18[3]={x[0]-c[0],x[1]-c[1]+epsilon,x[2]-c[2]+epsilon};
	
	double x19[3]={x[0]-c[0]-epsilon,x[1]-c[1]-epsilon,x[2]-c[2]+epsilon};
	double x20[3]={x[0]-c[0]-epsilon,x[1]-c[1]+epsilon,x[2]-c[2]+epsilon};
	double x21[3]={x[0]-c[0]+epsilon,x[1]-c[1]-epsilon,x[2]-c[2]+epsilon};
	double x22[3]={x[0]-c[0]+epsilon,x[1]-c[1]+epsilon,x[2]-c[2]+epsilon};
	
	double x23[3]={x[0]-c[0]-epsilon,x[1]-c[1]-epsilon,x[2]-c[2]-epsilon};
	double x24[3]={x[0]-c[0]-epsilon,x[1]-c[1]+epsilon,x[2]-c[2]-epsilon};
	double x25[3]={x[0]-c[0]+epsilon,x[1]-c[1]-epsilon,x[2]-c[2]-epsilon};
	double x26[3]={x[0]-c[0]+epsilon,x[1]-c[1]+epsilon,x[2]-c[2]-epsilon};

	n[0] = 4*Distance(x1)-4*Distance(x4)+2*Distance(x10)-2*Distance(x7)+2*Distance(x9)-2*Distance(x8)+2*Distance(x14)-2*Distance(x11)+2*Distance(x13)-2*Distance(x12)-Distance(x19)-Distance(x20)+Distance(x21)+Distance(x22)-Distance(x23)-Distance(x24)+Distance(x25)+Distance(x26);
	n[1] = 4*Distance(x2)-4*Distance(x5)+2*Distance(x10)-2*Distance(x7)+2*Distance(x8)-2*Distance(x9)+2*Distance(x18)-2*Distance(x15)+2*Distance(x17)-2*Distance(x16)-Distance(x19)+Distance(x20)-Distance(x21)+Distance(x22)-Distance(x23)+Distance(x24)-Distance(x25)+Distance(x26);
	n[2] = 4*Distance(x3)-4*Distance(x6)+2*Distance(x12)-2*Distance(x13)+2*Distance(x14)-2*Distance(x11)+2*Distance(x16)-2*Distance(x17)+2*Distance(x18)-2*Distance(x15)+Distance(x19)+Distance(x20)+Distance(x21)+Distance(x22)-Distance(x23)-Distance(x24)-Distance(x25)-Distance(x26);
	double len=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	for(int i=0;i<3;i++){
		n[i]/=(len);
	}
}
