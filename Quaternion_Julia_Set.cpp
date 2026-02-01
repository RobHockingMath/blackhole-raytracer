#include <math.h>
#include <iostream>
#include "Quaternion_Julia_Set.h"
#include "Sphere.h"
#include "invisibility.h"

#define PI 3.14159

using namespace std;

Quaternion_Julia_Set::Quaternion_Julia_Set(double c1,double c2,double c3,double R,double I,double theta,double E,int red,int green, int blue)
	: Distance_fractal(c1,c2,c3,E,red,green,blue)
{
   C[0]=R;C[1]=I;C[2]=0;C[3]=0;
   Theta=theta;
   e2[0]=cos(Theta);e2[1]=-sin(Theta);e2[2]=0;e2[3]=0;
   double e1[]={cos(Theta),sin(Theta),0,0};
   Mult(e1,C);
   is_wireframe=false;
}

bool Quaternion_Julia_Set::wireframe_condition(double p[]){
	return false;
}

//q2 is overwritten with q1*q2.  q1 is unchanged
void Quaternion_Julia_Set::Mult(double q1[],double q2[]){
	double temp[4];
	temp[0]=q1[0]*q2[0]-q1[1]*q2[1]-q1[2]*q2[2]-q1[3]*q2[3];
	temp[1]=q1[1]*q2[0]+q1[0]*q2[1]+q1[2]*q2[3]-q1[3]*q2[2];
	temp[2]=q1[0]*q2[2]+q1[2]*q2[0]+q1[3]*q2[1]-q1[1]*q2[3];
	temp[3]=q1[0]*q2[3]+q1[3]*q2[0]+q1[1]*q2[2]-q1[2]*q2[1];
	for(int i=0;i<4;i++){q2[i]=temp[i];}
}

//two sequences q_n and q'_n are relavent for the distance estimator
//this updates them both.
void Quaternion_Julia_Set::Next(double q[],double qp[]){
	//formula are q'_n+1=2*q_n*q'_n,
	//  q_n+1=q_n^2+C
	
	Mult(q,qp);
	Mult(q,q);
	Mult(e2,q);
	Mult(e2,qp);
	
	
	for(int i=0;i<4;i++){
		q[i]+=C[i];
		qp[i]*=2;
	}
	
	/**Mult(q,qp);
	Mult(q,q);
	for(int i=0;i<4;i++){
		q[i]+=C[i];
		qp[i]*=2;
	}**/
}

//returns a lower bound for the distance from point p to the set
double Quaternion_Julia_Set::Distance(double p[]){

	//use magic formula
	//d(q)=lim n->inf |qn| log |qn| /(2 |q'n|)

	double q[4]={p[0],p[1],p[2],0};
	double qp[4]={1,0,0,0};
	for(int i=0;i<30;i++){
		Next(q,qp);
		if(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]>5){break;}
	}
	
	double len_q=sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
	double len_qp=sqrt(qp[0]*qp[0]+qp[1]*qp[1]+qp[2]*qp[2]+qp[3]*qp[3]);

	//cout<<"here "<<len_qp<<" "<<len_q<<endl;
	if(len_q!=0 && len_qp!=0){
		return len_q*log(len_q)/(2*len_qp);}
	else{return 10;}
	

}
