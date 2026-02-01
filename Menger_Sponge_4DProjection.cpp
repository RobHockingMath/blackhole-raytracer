#include <math.h>
#include <iostream>
#include "Menger_Sponge_4DProjection.h"
#include "Sphere.h"
#include "invisibility.h"
#include <Eigen/Dense>

#define PI 3.14159
#define SQRT3 1.7321

using namespace std;

Menger_Sponge_4DProjection::Menger_Sponge_4DProjection(double c1,double c2,double c3,double scale_,int depth_,double E,double offset_,double n1, double n2, double n3, double n4,int red,int green, int blue)
	: Distance_fractal(c1,c2,c3,E,red,green,blue)
{
    depth = depth_;
    scale = scale_;
    offset = offset_;
   
    Eigen::MatrixXd A(1, 4);
    A(0,0)=n1;A(0,1)=n2;A(0,2)=n3;A(0,3)=n4;
   
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::MatrixXd nullSpaceBasis = svd.matrixV().rightCols(A.cols() - svd.rank());

    // Perform Gram-Schmidt orthogonalization
    Eigen::MatrixXd orthonormalBasis = nullSpaceBasis;
    for (int i = 0; i < orthonormalBasis.cols(); ++i) {
        for (int j = 0; j < i; ++j) {
            orthonormalBasis.col(i) -= orthonormalBasis.col(j).dot(nullSpaceBasis.col(i)) * orthonormalBasis.col(j);
        }
        orthonormalBasis.col(i).normalize();
    }
	
	e1[0]=orthonormalBasis(0,0);
	e1[1]=orthonormalBasis(1,0);
	e1[2]=orthonormalBasis(2,0);
	e1[3]=orthonormalBasis(3,0);
	
	e2[0]=orthonormalBasis(0,1);
	e2[1]=orthonormalBasis(1,1);
	e2[2]=orthonormalBasis(2,1);
	e2[3]=orthonormalBasis(3,1);
	
	e3[0]=orthonormalBasis(0,2);
	e3[1]=orthonormalBasis(1,2);
	e3[2]=orthonormalBasis(2,2);
	e3[3]=orthonormalBasis(3,2);
	
	double len_n = sqrt(n1*n1+n2*n2+n3*n3+n4*n4);
	
	N[0]=n1/len_n;
	N[1]=n2/len_n;
	N[2]=n3/len_n;
	N[3]=n4/len_n;

    // Print the orthonormal basis
    std::cout << "Orthonormal basis for the null space:\n" << orthonormalBasis << std::endl;

}


/**void Menger_Sponge_4DSlice::rot(double p[]){
	double x = p[0];
	double y = p[1];
	double z = p[2];
	
	double theta = PI/3;
	double c = cos(theta);
	double s = sin(theta);
	
	double xx = c*x-s*z;
	z = s*xx+c*z;
	x = xx;
	p[0]=x;
	p[1]=y;
	p[2]=z;
	
}**/


bool Menger_Sponge_4DProjection::wireframe_condition(double p[]){

	return false;

}

//returns a lower bound for the distance from point p to the set
double Menger_Sponge_4DProjection::Distance(double p[]){

	//this->rot(p);

	double X=p[0], Y=p[1], Z=p[2];
	
	X=X/scale;Y=Y/scale;Z=Z/scale; //center it by changing position and scale	
	
	double x = 0.5+X*e1[0]+Y*e2[0]+Z*e3[0];
	double y = 0.5+X*e1[1]+Y*e2[1]+Z*e3[1];
	double z = 0.5+X*e1[2]+Y*e2[2]+Z*e3[2];
	double w = 0.5+X*e1[3]+Y*e2[3]+Z*e3[3];
	
	/**double x = (X+Y-Z)/2.0;
	double y = (-X-Y-Z)/2.0;
	double z = (X-Y+Z)/2.0;
	double w = offset+(-X+Y+Z)/2.0;**/
	
	double offset_p = offset/2+2;
	
	//x = (1+x)/2;y = (1+y)/2; z = (1+z)/2;w=(1+w)/2;

	// find the closest point on the box.

	double x_c = min(max(x,0.0),1.0);
	double y_c = min(max(y,0.0),1.0);
	double z_c = min(max(z,0.0),1.0);
	double w_c = min(max(w,0.0),1.0);
	
	double delta[4]={x_c-x,y_c-y,z_c-z,w_c-w};
	double delta_dot_n = delta[0]*N[0]+delta[1]*N[1]+delta[2]*N[2]+delta[3]*N[3]+delta[4]*N[4];
	double dist = (delta[0]-delta_dot_n*N[0])*(delta[0]-delta_dot_n*N[0])+(delta[1]-delta_dot_n*N[1])*(delta[1]-delta_dot_n*N[1])+(delta[2]-delta_dot_n*N[2])*(delta[2]-delta_dot_n*N[2])+(delta[3]-delta_dot_n*N[3])*(delta[3]-delta_dot_n*N[3])+(delta[4]-delta_dot_n*N[4])*(delta[4]-delta_dot_n*N[4]);
	
	dist = sqrt(dist);
	
	double d=dist;

	double pp=1.0;
	for (int i=1; i<=0; ++i) {
		double xa = fmod(3.0*x*pp,3.0);
		double ya = fmod(3.0*y*pp,3.0);
		double za = fmod(3.0*z*pp,3.0);
		//double wa = fmod(3.0*w*pp,3.0);
		double wa = fmod(3.0*w*pp,3.0);
		pp*=3.0;

		//we can also translate/rotate (xa,ya,za) without affecting the DE estimate

		//double xx=0.5-abs(xa-1.5), yy=0.5-abs(ya-1.5), zz=0.5-abs(za-1.5), ww=0.5-abs(wa-1.5);
		//d1=min(max(xx,zz),min(max(xx,yy),max(yy,zz))) / pp; //w ommitted
		
		// first omit w.
		
		x_c = min(max(xa,1.0),2.0);
		y_c = min(max(ya,1.0),2.0);
		z_c = min(max(za,1.0),2.0);
		w_c = min(max(wa,1.0),2.0);
		
		delta[0]=x_c-xa;delta[1]=y_c-ya;delta[2]=z_c-za;delta[3]=0;
		delta_dot_n = delta[0]*N[0]+delta[1]*N[1]+delta[2]*N[2]+delta[3]*N[3]+delta[4]*N[4];
		double dist_w = (delta[0]-delta_dot_n*N[0])*(delta[0]-delta_dot_n*N[0])+(delta[1]-delta_dot_n*N[1])*(delta[1]-delta_dot_n*N[1])+(delta[2]-delta_dot_n*N[2])*(delta[2]-delta_dot_n*N[2])+(delta[3]-delta_dot_n*N[3])*(delta[3]-delta_dot_n*N[3])+(delta[4]-delta_dot_n*N[4])*(delta[4]-delta_dot_n*N[4]);		
		dist_w=dist_w/pp;
		
		// now omit z
		
		delta[0]=x_c-xa;delta[1]=y_c-ya;delta[2]=0;delta[3]=w_c-wa;
		delta_dot_n = delta[0]*N[0]+delta[1]*N[1]+delta[2]*N[2]+delta[3]*N[3]+delta[4]*N[4];
		double dist_z = (delta[0]-delta_dot_n*N[0])*(delta[0]-delta_dot_n*N[0])+(delta[1]-delta_dot_n*N[1])*(delta[1]-delta_dot_n*N[1])+(delta[2]-delta_dot_n*N[2])*(delta[2]-delta_dot_n*N[2])+(delta[3]-delta_dot_n*N[3])*(delta[3]-delta_dot_n*N[3])+(delta[4]-delta_dot_n*N[4])*(delta[4]-delta_dot_n*N[4]);		
		dist_z=dist_z/pp;
		
		// now omit y
		
		delta[0]=x_c-xa;delta[1]=0;delta[2]=z_c-za;delta[3]=w_c-wa;
		delta_dot_n = delta[0]*N[0]+delta[1]*N[1]+delta[2]*N[2]+delta[3]*N[3]+delta[4]*N[4];
		double dist_y = (delta[0]-delta_dot_n*N[0])*(delta[0]-delta_dot_n*N[0])+(delta[1]-delta_dot_n*N[1])*(delta[1]-delta_dot_n*N[1])+(delta[2]-delta_dot_n*N[2])*(delta[2]-delta_dot_n*N[2])+(delta[3]-delta_dot_n*N[3])*(delta[3]-delta_dot_n*N[3])+(delta[4]-delta_dot_n*N[4])*(delta[4]-delta_dot_n*N[4]);		
		dist_y=dist_y/pp;
		
		// now omit x
		
		delta[0]=0;delta[1]=y_c-ya;delta[2]=z_c-za;delta[3]=w_c-wa;
		delta_dot_n = delta[0]*N[0]+delta[1]*N[1]+delta[2]*N[2]+delta[3]*N[3]+delta[4]*N[4];
		double dist_x = (delta[0]-delta_dot_n*N[0])*(delta[0]-delta_dot_n*N[0])+(delta[1]-delta_dot_n*N[1])*(delta[1]-delta_dot_n*N[1])+(delta[2]-delta_dot_n*N[2])*(delta[2]-delta_dot_n*N[2])+(delta[3]-delta_dot_n*N[3])*(delta[3]-delta_dot_n*N[3])+(delta[4]-delta_dot_n*N[4])*(delta[4]-delta_dot_n*N[4]);		
		dist_x=dist_x/pp;
		
		
		
		
		
		d=max(dist,max(dist_w,max(dist_z,max(dist_y,dist_x))));
	}
	//return d*2.0; //the distance estimate. The *2 is because of the scaling we did at the beginning of the function
	return d*scale;
	

}