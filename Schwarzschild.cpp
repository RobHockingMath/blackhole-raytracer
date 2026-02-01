#include "Schwarzschild.h"

Schwarzschild::Schwarzschild(double pos_[3],double tol_, double mass, int dim)
	: Spacetime(tol_)
{
	M=mass;
	n=dim;
	pos[0]=pos_[0];
	pos[1]=pos_[1];
	pos[2]=pos_[2];
}

int Schwarzschild::inside_universe(double *p){
	double r=sqrt((p[0]-pos[0])*(p[0]-pos[0])+(p[1]-pos[1])*(p[1]-pos[1])+(p[2]-pos[2])*(p[2]-pos[2]));
	//event horizon at r2 = M^2
	if(r>1.00*M){
		return 1;
	}
	else{
		return 0;
	}
}

double Schwarzschild::get_h0(double *p){
	double r=sqrt((p[0]-pos[0])*(p[0]-pos[0])+(p[1]-pos[1])*(p[1]-pos[1])+(p[2]-pos[2])*(p[2]-pos[2]));
	return min(0.1,fabs(1-r/M)); // take a big step.
}

// We think of p as being a point and x a vector.  
// This format is useful as it can be used directly in adaptive Runge Kutta routines.

double Schwarzschild::error_norm(double *p,double *x){
	//return sqrt((x[0]-pos[0])*(x[0]-pos[0])+(x[1]-pos[1])*(x[1]-pos[1])+(x[2]-pos[2])*(x[2]-pos[2]));
	return sqrt((x[0])*(x[0])+(x[1])*(x[1])+(x[2])*(x[2]));
	//double r=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
	//return (fabs(x[0])+fabs(x[1])+fabs(x[2]))/fabs(1-r/M);
	
}

bool Schwarzschild::accel_1st_order(double*& yp,double* y, double E){
	
	// y and yp are giant vectors of the form y=[p,v] where p is position and v is velocity.
	
	double r=sqrt((y[0]-pos[0])*(y[0]-pos[0])+(y[1]-pos[1])*(y[1]-pos[1])+(y[2]-pos[2])*(y[2]-pos[2]));
	//double U = 1+M/r;
	//double U = 1+M/pow(r,n-3);
	double r_yz = sqrt((y[1]-pos[1])*(y[1]-pos[1])+(y[2]-pos[2])*(y[2]-pos[2]));
	//double dUdx = -M/(r*r*r)*y[0];
	//double dUdy = -M/(r*r*r)*y[1];
	//double dUdz = -M/(r*r*r)*y[2];
	//double dUdx = -(n-3)*(M/pow(r,n-1))*y[0];
	//double dUdy = -(n-3)*(M/pow(r,n-1))*y[1];
	//double dUdz = -(n-3)*(M/pow(r,n-1))*y[2];
	
	//double dU_dot_v = dUdx*y[3]+dUdy*y[4]+dUdz*y[5];
	
	double speed2 = y[3]*y[3]+y[4]*y[4]+y[5]*y[5];
	
	// position derivative is velocity.
	yp[0]=y[3]; 
	yp[1]=y[4]; 
	yp[2]=y[5]; 
	
	double r0 = 2.0;
	
	//double MM = fmax(M*tanh(r0-r_yz),0.01);
	double MM = M*tanh(r0-r);
	//double MM = M;
	if(MM<=0 && y[0]>=2){
		MM=0;
	}
	//if(MM<0 && y[0]>=2){
	//	MM=0;
	//}
	//double MM = M;
	
	// acceleration is zero.
	if(n==4){
		double U = 1+MM/r;
		double dUdx = -MM*(y[0]-pos[0])/(r*r*r);
		double dUdy = -MM*(y[1]-pos[1])/(r*r*r);
		double dUdz = -MM*(y[2]-pos[2])/(r*r*r);
		double dU_dot_v = dUdx*y[3]+dUdy*y[4]+dUdz*y[5];
		yp[3]=((E*E+speed2)*dUdx-2*dU_dot_v*y[3])/U;
		yp[4]=((E*E+speed2)*dUdy-2*dU_dot_v*y[4])/U;
		yp[5]=((E*E+speed2)*dUdz-2*dU_dot_v*y[5])/U;
		//if(r<0.0001){
		//	yp[3]=0;
		//	yp[4]=0;
		//	yp[5]=0;
		//}
	}
	else{
		double U = 1+M/pow(r,n-3);
		double dUdx = -(n-3)*(M/pow(r,n-1))*(y[0]-pos[0]);
		double dUdy = -(n-3)*(M/pow(r,n-1))*(y[1]-pos[1]);
		double dUdz = -(n-3)*(M/pow(r,n-1))*(y[2]-pos[2]);
		double dU_dot_v = dUdx*y[3]+dUdy*y[4]+dUdz*y[5];
		double U_to_alpha=pow(U,2.0/(n-3)-1);
		double U_to_alpha_m1 = U_to_alpha/U;
		yp[3]=((E*E+(U_to_alpha_m1/(n-3))*speed2)*dUdx-(2.0/(n-3))*U_to_alpha_m1*dU_dot_v*y[3])/U_to_alpha;
		yp[4]=((E*E+(U_to_alpha_m1/(n-3))*speed2)*dUdy-(2.0/(n-3))*U_to_alpha_m1*dU_dot_v*y[4])/U_to_alpha;
		yp[5]=((E*E+(U_to_alpha_m1/(n-3))*speed2)*dUdz-(2.0/(n-3))*U_to_alpha_m1*dU_dot_v*y[5])/U_to_alpha;
	}
	
	return inside_universe(y);
	
}

double Schwarzschild::compute_energy(double *p){
	double r=sqrt((p[0]-pos[0])*(p[0]-pos[0])+(p[1]-pos[1])*(p[1]-pos[1])+(p[2]-pos[2])*(p[2]-pos[2]));
	//double U = 1+M/r;
	double U = 1+M/pow(r,n-3);
	return 1.0/U;
}

void Schwarzschild::scale_velocity(double v[],double *p){
	double r=sqrt((p[0]-pos[0])*(p[0]-pos[0])+(p[1]-pos[1])*(p[1]-pos[1])+(p[2]-pos[2])*(p[2]-pos[2]));
	//double U = 1.0+M/r;
	double U = 1+M/pow(r,n-3);
	v[0]=v[0]/pow(U,1.0/(n-3));
	v[1]=v[1]/pow(U,1.0/(n-3));
	v[2]=v[2]/pow(U,1.0/(n-3));
}




