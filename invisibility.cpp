#include "invisibility.h"

void expand(double x[],double R,double X[]){
	double xdotx=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
	double s=sqrt(1+R*R/xdotx);
	for(int i=0;i<3;i++){
		X[i]=s*x[i];
	}
}

void expand_derivative(double x[],double xp[],double R,double Xp[]){
	double xdotx=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
	double xdotxp=x[0]*xp[0]+x[1]*xp[1]+x[2]*xp[2];
	double s=sqrt(1+R*R/xdotx);
	for(int i=0;i<3;i++){
		Xp[i]=s*xp[i]-R*R*xdotxp*x[i]/(s*xdotx*xdotx);
	}
}

void contract(double X[],double R,double x[]){
	double XdotX=X[0]*X[0]+X[1]*X[1]+X[2]*X[2];
	double s=sqrt(1-R*R/XdotX);
	for(int i=0;i<3;i++){
		x[i]=s*X[i];
	}
}

void contract_derivative(double X[],double Xp[],double R,double xp[]){
	double XdotX=X[0]*X[0]+X[1]*X[1]+X[2]*X[2];
	double XdotXp=X[0]*Xp[0]+X[1]*Xp[1]+X[2]*Xp[2];
	double s=sqrt(1-R*R/XdotX);
	for(int i=0;i<3;i++){
		xp[i]=s*Xp[i]+R*R*XdotXp*X[i]/(s*XdotX*XdotX);
	}
}

void initial_value(double P[],double V[],double R,double p[],double v[]){
	contract(P,R,p);
	contract_derivative(P,V,R,v);
}

void boundary_value(double P1[],double P2[],double R,double p[],double v[]){
	contract(P1,R,p);
	contract(P2,R,v);
	v[0]=v[0]-p[0];
	v[1]=v[1]-p[1];
	v[2]=v[2]-p[2];
}

//on exit, V holds the tangent vector at P1 of the unique
//light ray from P1 to P2.
void get_direction(double P1[],double P2[],double R,double V[]){
	double p[3],v[3];
	boundary_value(P1,P2,R,p,v);
	expand_derivative(p,v,R,V);
}

//At entry, V is assumed to be the tangent
//to the light ray at P.  On exit, it is the tangent
//distance t down the curve in the +V direction.
void transport_tangent(double P[],double t,double R,double V[]){
	double p[3],v[3];
	contract(P,R,p);
	contract_derivative(P,V,R,v);
	for(int i=0;i<3;i++){p[i]+=t*v[i];}
	expand_derivative(p,v,R,V);	
}

void flow(double P[],double V[],double s, double R){
	double p[3],v[3],norm;
	contract(P,R,p);
	contract_derivative(P,V,R,v);

	double dt=0.01;
	double dist=0;
	while(dist<s){
		for(int i=0;i<3;i++){p[i]+=dt*v[i];}
		expand_derivative(p,v,R,V);
		norm=sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);
		dist+=norm*dt;
	}
	expand(p,R,P);
}
