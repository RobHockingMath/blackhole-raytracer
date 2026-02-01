#include <math.h>
#include <iostream>
#include "Menger_Sponge_4DSlice.h"
#include "Sphere.h"
#include "invisibility.h"

#define PI 3.14159
#define SQRT3 1.7321

using namespace std;

Menger_Sponge_4DSlice::Menger_Sponge_4DSlice(double c1,double c2,double c3,double scale_,int depth_,int type_,double E,double offset_,int red,int green, int blue)
	: Distance_fractal(c1,c2,c3,E,red,green,blue)
{
   depth = depth_;
   scale = scale_;
   offset = offset_;
   type = type_;
}


void Menger_Sponge_4DSlice::rot(double p[]){
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
	
}


bool Menger_Sponge_4DSlice::wireframe_condition(double p[]){

	double X=p[0]-(this->c[0]), Y=p[1]-(this->c[1]), Z=p[2]-(this->c[2]);
	
	//double X=p[0], Y=p[1], Z=p[2];
	
	X=X/scale;Y=Y/scale;Z=Z/scale; //center it by changing position and scale	
	
	//double x = 0.25*offset+(X+Y-Z)/2.0;
	//double y = 0.25*offset+(-X-Y-Z)/2.0;
	//double z = 0.25*offset+(X-Y+Z)/2.0;
	//double w = 0.25*offset+(-X+Y+Z)/2.0;
	double x = 0.25*offset+(sqrt(8/9.0)*X+0*Y+(-1/3.0)*Z)*sqrt(3/4.0);
	double y = 0.25*offset+(-sqrt(2/9.0)*X+sqrt(2/3.0)*Y+(-1/3.0)*Z)*sqrt(3/4.0);
	double z = 0.25*offset+(-sqrt(2/9.0)*X-sqrt(2/3.0)*Y+(-1/3.0)*Z)*sqrt(3/4.0);
	double w = 0.25*offset+(0*X+0*Y+1*Z)*sqrt(3/4.0);
	/**double x = (X+Y-Z)/2.0;
	double y = (-X-Y-Z)/2.0;
	double z = (X-Y+Z)/2.0;
	double w = offset+(-X+Y+Z)/2.0;**/
	
	//x = (1+x)/2;y = (1+y)/2; z = (1+z)/2;w=(1+w)/2;
	double dx = min(abs(1+x),abs(1-x));
	double dy = min(abs(1+y),abs(1-y));
	double dz = min(abs(1+z),abs(1-z));
	double dw = min(abs(1+w),abs(1-w));
	double dx2 = min(abs(3+x),abs(3-x));
	double dy2 = min(abs(3+y),abs(3-y));
	double dz2 = min(abs(3+z),abs(3-z));
	double dw2 = min(abs(3+w),abs(3-w));
	dx = min(dx,dx2);
	dy = min(dy,dy2);
	dz = min(dz,dz2);
	dw = min(dw,dw2);
	
	double delta = 0.02*2.0/3;
	
	if((dx<delta && dy<delta) || (dx<delta && dz<delta) || (dx<delta && dw<delta) || (dy<delta && dz<delta) || (dy<delta && dw<delta) || (dz<delta && dw<delta)){
		return true;
	}
	else{
		return false;
	}

}



//returns a lower bound for the distance from point p to the set
double Menger_Sponge_4DSlice::Distance(double p[]){

	//this->rot(p);

	double X=p[0], Y=p[1], Z=p[2];
	
	X=X/scale;Y=Y/scale;Z=Z/scale; //center it by changing position and scale	
	
	//double x = 0.25*offset+(X+Y-Z)/2.0;
	//double y = 0.25*offset+(-X-Y-Z)/2.0;
	//double z = 0.25*offset+(X-Y+Z)/2.0;
	//double w = 0.25*offset+(-X+Y+Z)/2.0;
	
	double x = 0.25*offset+(sqrt(8/9.0)*X+0*Y+(-1/3.0)*Z)*sqrt(3/4.0);
	double y = 0.25*offset+(-sqrt(2/9.0)*X+sqrt(2/3.0)*Y+(-1/3.0)*Z)*sqrt(3/4.0);
	double z = 0.25*offset+(-sqrt(2/9.0)*X-sqrt(2/3.0)*Y+(-1/3.0)*Z)*sqrt(3/4.0);
	double w = 0.25*offset+(0*X+0*Y+1*Z)*sqrt(3/4.0);
	
	/**double x = (X+Y-Z)/2.0;
	double y = (-X-Y-Z)/2.0;
	double z = (X-Y+Z)/2.0;
	double w = offset+(-X+Y+Z)/2.0;**/
	
	double offset_p = offset/2+2;
	
	x = (1+x)/2;y = (1+y)/2; z = (1+z)/2;w=(1+w)/2;
	
	double xx=abs(x-0.5)-0.5, yy=abs(y-0.5)-0.5, zz=abs(z-0.5)-0.5, ww=abs(w-0.5)-0.5;//ww=offset_p-x-y-z;//ww=abs(offset_p-x-y-z-0.5)-0.5;//ww=abs(w-0.5)-0.5;
	
	double d1;
	d1=max(xx,max(yy,max(zz,ww))); //distance to the box
	//double dxyz = max(xx,max(yy,zz)); // w omitted
	//double dxyw = max(xx,max(yy,ww)); // z omitted
	//double dxzw = max(xx,max(zz,ww)); // y omitted
	//double dyzw = max(yy,max(zz,ww)); // x omitted
	//d1 = min(dxyz,min(dxyw,min(dxzw,dyzw)));
	//d1=min(max(xx,zz),min(max(xx,yy),max(yy,zz))); //w ommitted
	//d1 = max(d1,min(max(xx,ww),min(max(xx,yy),max(yy,ww)))); // z ommitted
	//d1 = max(d1,min(max(xx,zz),min(max(xx,ww),max(ww,zz)))); // y ommitted 
	//d1 = max(d1,min(max(ww,zz),min(max(ww,yy),max(yy,zz)))); // x ommitted
	
	//double xxx=abs(x-0.5)-1.5, yyy=abs(y-0.5)-1.5, zzz=abs(z-0.5)-1.5, www=abs(w-0.5)-1.5;
	//d1 = min(d1,min(xxx,min(yyy,min(zzz,www))));
	//xxx=abs(x-0.5)-0.5; yyy=abs(y-1.5)-0.5; zzz=abs(z-0.5)-0.5; www=abs(w-0.5)-0.5;
	//d1 = min(d1,min(xxx,min(yyy,min(zzz,www))));
	//xxx=abs(x-0.5)-0.5; yyy=abs(y-0.5)-0.5; zzz=abs(z-1.5)-0.5; www=abs(w-0.5)-0.5;
	//d1 = min(d1,min(xxx,min(yyy,min(zzz,www))));
	//xxx=abs(x-0.5)-0.5; yyy=abs(y-0.5)-0.5; zzz=abs(z-0.5)-0.5; www=abs(w-1.5)-0.5;
	//d1 = min(d1,min(xxx,min(yyy,min(zzz,www))));



	/**double d=d1; //current computed distance
	double pp=1.0;
	for (int i=1; i<=depth; ++i) {
		
		//fmod fucks shit up for negative numbers.
		while(x<0){x=1+x;}
		while(y<0){y=1+y;}
		while(z<0){z=1+z;}
		while(w<0){w=1+w;}
		
		double xa = fmod(3.0*x*pp,3.0);
		double ya = fmod(3.0*y*pp,3.0);
		double za = fmod(3.0*z*pp,3.0);
		//double wa = fmod(3.0*w*pp,3.0);
		double wa = fmod(3.0*w*pp,3.0);
		pp*=3.0;

		//we can also translate/rotate (xa,ya,za) without affecting the DE estimate

		double xx=0.5-abs(xa-1.5), yy=0.5-abs(ya-1.5), zz=0.5-abs(za-1.5), ww=0.5-abs(wa-1.5);
		//double xx=0.5-abs(xa-0.5), yy=0.5-abs(ya-0.5), zz=0.5-abs(za-0.5), ww=0.5-abs(wa-0.5);
		d1=min(max(xx,zz),min(max(xx,yy),max(yy,zz))) / pp; //w ommitted
		d1 = max(d1,min(max(xx,ww),min(max(xx,yy),max(yy,ww))) / pp); // z ommitted
		d1 = max(d1,min(max(xx,zz),min(max(xx,ww),max(ww,zz))) / pp); // y ommitted 
		d1 = max(d1,min(max(ww,zz),min(max(ww,yy),max(yy,zz))) / pp); // x ommitted


		d=max(d,d1); //intersection
	}**/
	double d=d1; //current computed distance
	double pp=1.0;
	
	for (int i=1; i<=depth; ++i) {
		
		
		
		//fmod fucks shit up for negative numbers.
		double x1=x,y1=y,z1=z,w1=w;
		
		/**if(i>1){
		while(x1<0){x1=1+x1;x_neg=true;}
		while(y1<0){y1=1+y1;y_neg=true;}
		while(z1<0){z1=1+z1;z_neg=true;}
		while(w1<0){w1=1+w1;w_neg=true;}
		}**/
		
		
		// x % k = x - k*floor(x/k)
		
		double xa =3.0*x1*pp-floor(3.0*x1*pp/3.0)*3.0;
		double ya =3.0*y1*pp-floor(3.0*y1*pp/3.0)*3.0;
		double za =3.0*z1*pp-floor(3.0*z1*pp/3.0)*3.0;
		double wa =3.0*w1*pp-floor(3.0*w1*pp/3.0)*3.0;
		
		//double xa = fmod(3.0*x1*pp,3.0);
		//double ya = fmod(3.0*y1*pp,3.0);
		//double za = fmod(3.0*z1*pp,3.0);
		//double wa = fmod(3.0*w*pp,3.0);
		//double wa = fmod(3.0*w1*pp,3.0);
		pp*=3.0;

		//we can also translate/rotate (xa,ya,za) without affecting the DE estimate

		double xx=0.5-abs(xa-1.5), yy=0.5-abs(ya-1.5), zz=0.5-abs(za-1.5), ww=0.5-abs(wa-1.5);
		
		/**if(i==1){
			if(x1>=1){xx=-100;}
			if(y1>=1){yy=-100;}
			if(z1>=1){zz=-100;}
			if(w1>=1){ww=-100;}
		}**/
		
		//d1=min(max(xx,zz),min(max(xx,yy),max(yy,zz))) / pp; //distance inside the 3 axis-aligned square tubes
		
		//double xx=0.5-abs(xa-0.5), yy=0.5-abs(ya-0.5), zz=0.5-abs(za-0.5), ww=0.5-abs(wa-0.5);
		d1=0;
		//if(!(w_neg && i==2) && !(x_lt_m1 || y_lt_m1 || z_lt_m1 || w_lt_m1)){
		d1 = max(d1,min(max(xx,zz),min(max(xx,yy),max(yy,zz))) / pp);//} //w ommitted
		//if(!(z_neg && i==2) && !(x_lt_m1 || y_lt_m1 || z_lt_m1 || w_lt_m1)){
		d1 = max(d1,min(max(xx,ww),min(max(xx,yy),max(yy,ww))) / pp);//} // z ommitted
		//if(!(y_neg && i==2) && !(x_lt_m1 || y_lt_m1 || z_lt_m1 || w_lt_m1)){
		d1 = max(d1,min(max(xx,zz),min(max(xx,ww),max(ww,zz))) / pp);//} // y ommitted 
		//if(!(x_neg && i==2) && !(x_lt_m1 || y_lt_m1 || z_lt_m1 || w_lt_m1)){
		d1 = max(d1,min(max(ww,zz),min(max(ww,yy),max(yy,zz))) / pp);//} // x ommitted


		d=max(d,d1); //intersection
	}
	//return d*2.0; //the distance estimate. The *2 is because of the scaling we did at the beginning of the function
	
	//d = max(d,max(xxx,max(yyy,max(zzz,www))));
	
	return d*scale;
	

}

/**void Menger_Sponge_4DSlice::get_Normal(double x[], double n[]){
	
	double h = epsilon;
	double v[3]={n[0],n[1],n[2]};
	double x1[3]={x[0]+h,x[1],x[2]};double x1p[3]={0,0,0};
	double x2[3]={x[0],x[1]+h,x[2]};double x2p[3]={0,0,0};
	double x3[3]={x[0],x[1],x[2]+h};double x3p[3]={0,0,0};
	
	this->get_Intersection(x1,n,x1p);
	this->get_Intersection(x2,n,x2p);
	this->get_Intersection(x3,n,x3p);
	
	double x01[3]={x1p[0]-x[0],x1p[1]-x[1],x1p[2]-x[2]};
	double x02[3]={x2p[0]-x[0],x2p[1]-x[1],x2p[2]-x[2]};
	n[0]=x01[1]*x02[2]-x02[1]*x01[2];
	n[1]=x01[2]*x02[0]-x02[2]*x01[0];
	n[2]=x01[0]*x02[1]-x02[1]*x01[0];
	double len = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	n[0]=n[0]/len;
	n[1]=n[1]/len;
	n[2]=n[2]/len;
	double n_dot_v = n[0]*v[0]+n[1]*v[1]+n[2]*v[2];
	if(n_dot_v<0){
		n[0]=-n[0];
		n[1]=-n[1];
		n[2]=-n[2];
	}
	
	
	/**double x1[3]={x[0],x[1],x[2]};
	double x2[3]={x[0],x[1],x[2]};
	double xp[3]={x1[0]+epsilon,x1[1]+epsilon,x1[2]+epsilon};
	double xm[3]={x1[0]-epsilon,x1[1]-epsilon,x1[2]-epsilon};
	
	
	double n1[3]={1.0/SQRT3 , 1.0/SQRT3 , -1.0/SQRT3};
    double n2[3]={-1.0/SQRT3, -1.0/SQRT3 ,-1.0/SQRT3};
    double n3[3]={1.0/SQRT3, -1.0/SQRT3 , 1.0/SQRT3};
    double n4[3]={-1.0/SQRT3,1.0/SQRT3,1.0/SQRT3};


	double min_dist = this->Distance(x);

	double n_c[3]={0,0,0};
	double x_p_n[3]={0,0,0};
	double x_m_n[3]={0,0,0};
	double dist_p,dist_m;

	for(int i=0;i<4;i++){
		if(i==0){
			n_c[0]=n1[0];n_c[1]=n1[1];n_c[2]=n1[2];
		}
		if(i==1){
			n_c[0]=n2[0];n_c[1]=n2[1];n_c[2]=n2[2];
		}
		if(i==2){
			n_c[0]=n3[0];n_c[1]=n3[1];n_c[2]=n3[2];
		}
		if(i==3){
			n_c[0]=n4[0];n_c[1]=n4[1];n_c[2]=n4[2];
		}
		x_p_n[0]=x[0]+0.5*epsilon*n_c[0];
		x_p_n[1]=x[1]+0.5*epsilon*n_c[1];
		x_p_n[2]=x[2]+0.5*epsilon*n_c[2];
		x_m_n[0]=x[0]-0.5*epsilon*n_c[0];
		x_m_n[1]=x[1]-0.5*epsilon*n_c[1];
		x_m_n[2]=x[2]-0.5*epsilon*n_c[2];
		dist_p = this->Distance(x_p_n);
		//if(this->is_Outside(x_p_n)==-1){dist_p=0;}
		dist_m = this->Distance(x_m_n);
		//if(this->is_Outside(x_m_n)==-1){dist_m=0;}
		if(dist_p<min_dist ){min_dist=dist_p;n[0]=n_c[0];n[1]=n_c[1];n[2]=n_c[2];}
		if(dist_m<min_dist ){min_dist=dist_m;n[0]=-n_c[0];n[1]=-n_c[1];n[2]=-n_c[2];}
		n[0]=-n[0];
		n[1]=-n[1];
		n[2]=-n[2];
		
	}**/

	/**for(int i=0;i<3;i++){
		x1[i]+=epsilon;
		x2[i]-=epsilon;
		n[i]=this->Distance(x1)-this->Distance(x2);
		x1[i]-=epsilon;
		x2[i]+=epsilon;
	}
	double len=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	for(int i=0;i<3;i++){
		n[i]/=(len);
	}
	double x_p_n[3]={x[0]+2*epsilon*n[0],x[1]+2*epsilon*n[1],x[2]+2*epsilon*n[2]};
	if(this->is_Outside(x_p_n)==1){
		//n[0]=-n[0];
		//n[1]=-n[1];
		//n[2]=-n[2];
	}
	else{
		n[0]=-n[0];
		n[1]=-n[1];
		n[2]=-n[2];
	}
	return;**/
	
//}


double Menger_Sponge_4DSlice::is_Outside(double p[]){
	
	double X=p[0], Y=p[1], Z=p[2];
	
	X=X/scale;Y=Y/scale;Z=Z/scale; //center it by changing position and scale	
	
	//double x = 0.25*offset+(X+Y-Z)/2.0;
	//double y = 0.25*offset+(-X-Y-Z)/2.0;
	//double z = 0.25*offset+(X-Y+Z)/2.0;
	//double w = 0.25*offset+(-X+Y+Z)/2.0;
	
	double x = 0.25*offset+(sqrt(8/9.0)*X+0*Y+(-1/3.0)*Z)*sqrt(3/4.0);
	double y = 0.25*offset+(-sqrt(2/9.0)*X+sqrt(2/3.0)*Y+(-1/3.0)*Z)*sqrt(3/4.0);
	double z = 0.25*offset+(-sqrt(2/9.0)*X-sqrt(2/3.0)*Y+(-1/3.0)*Z)*sqrt(3/4.0);
	double w = 0.25*offset+(0*X+0*Y+1*Z)*sqrt(3/4.0);
	
	
	x = (1+x)/2;y = (1+y)/2; z = (1+z)/2; w = (1+w)/2;
	
	if(x<0 || x>1 || y<0 || y> 1 || z<0 || z>1 || w<0 || w> 1){
		return 1;
	}
	
	int base = 3;
	int middle = floor(base/2.0);
	
	double px,py,pz,pw;
	int ax,ay,az,aw;
	
	bool inside = true;
	for(int i=1;i<=depth;i++){
		int middle_sum = 0;
		if(i==1){
			px=base*x;
			py=base*y;
			pz=base*z;
			pw=base*w;
			ax=floor(px);
			ay=floor(py);
			az=floor(pz);
			aw=floor(pw);
		}
		else{
			px = base*(px-ax);
			ax = floor(px);
			py = base*(py-ay);
			ay = floor(py);
			pz = base*(pz-az);
			az = floor(pz);
			pw = base*(pw-aw);
			aw = floor(pw);
		}
		// cavity Menger Sponge
		//if((ax==1 && ay==1 && az==1)){
		//	inside=false;
		//}
			// standard Menger Sponge
		if(ax==middle){
			middle_sum=middle_sum+1;
		}
		if(ay==middle){
			middle_sum=middle_sum+1;
		}
		if(az==middle){
			middle_sum=middle_sum+1;
		}
		if(aw==middle){
			middle_sum=middle_sum+1;
		}
		if(middle_sum>1){
			inside=false;
		}
	}
	if(inside){return -1;}
	else{return 1;}
	
}

//use sphere texture coords
void Menger_Sponge_4DSlice::get_Texture_Coords(double x[], double u[]){
	double n[3]={0,0,0};
	
	this->get_Normal(x, n);
	double X[3];X[0]=x[0]-c[0];X[1]=x[1]-c[1];X[2]=x[2]-c[2];
	/**double r=sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
	u[0]=(atan2(X[1],X[0])+PI)/(2*PI)+r;
	u[1]=acos(X[2]/r)/PI+r;
	//if(u[0]<0){u[0]=0;}
	//if(u[1]<0){u[1]=0;}
	//if(u[0]>1){u[0]=1;}
	//if(u[1]>1){u[1]=1;}
	u[0]=u[0]-floor(u[0]/1)*1;
	u[1]=u[1]-floor(u[1]/1)*1;**/
	
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
	
	double s = scale*scale;
	
	if(d_min == d12 || d_min == d12m){
		u[0]=(X[0]*e12[0]+X[1]*e12[1]+X[2]*e12[2]-(-s*2))/(s*4);
		u[1]=(X[0]*u12[0]+X[1]*u12[1]+X[2]*u12[2]-(-s*2))/(s*4);
	}
	else if(d_min == d23 || d_min == d23m){
		u[0]=(X[0]*e23[0]+X[1]*e23[1]+X[2]*e23[2]-(-s*2))/(s*4);
		u[1]=(X[0]*u23[0]+X[1]*u23[1]+X[2]*u23[2]-(-s*2))/(s*4);
	}
	else if(d_min == d31 || d_min == d31m){
		u[0]=(X[0]*e31[0]+X[1]*e31[1]+X[2]*e31[2]-(-s*2))/(s*4);
		u[1]=(X[0]*u31[0]+X[1]*u31[1]+X[2]*u31[2]-(-s*2))/(s*4);
	}
	else if(d_min == d4 || d_min == d4m){
		u[0]=(X[0]-(-s*2))/(s*4);
		u[1]=(X[1]-(-s*2))/(s*4);
	}
	u[0]=u[0]-floor(u[0]/1)*1;
	u[1]=u[1]-floor(u[1]/1)*1;
	
	
	/**double d1 = (n[0]-1)*(n[0]-1)+n[1]*n[1]+n[2]*n[2];
	double d2 = (n[0]+1)*(n[0]+1)+n[1]*n[1]+n[2]*n[2];
	double d3 = n[0]*n[0]+(n[1]-1)*(n[1]-1)+n[2]*n[2];
	double d4 = n[0]*n[0]+(n[1]+1)*(n[1]+1)+n[2]*n[2];
	double d5 = n[0]*n[0]+n[1]*n[1]+(n[2]-1)*(n[2]-1);
	double d6 = n[0]*n[0]+n[1]*n[1]+(n[2]+1)*(n[2]+1);
	if(min(d1,d2)<min(min(d3,d4),min(d5,d6))){
		u[0]=(X[1]-(-2))/4.0;
		u[1]=(X[2]-(-2))/4.0;
	}
	else if(min(d3,d4)<min(min(d5,d6),min(d1,d2))){
		u[0]=(X[0]-(-2))/4.0;
		u[1]=(X[2]-(-2))/4.0;		
	}
	else if(min(d5,d6)<min(min(d1,d2),min(d3,d4))){
		u[0]=(X[0]-(-2))/4.0;
		u[1]=(X[1]-(-2))/4.0;			
	}
	u[0]=u[0]-floor(u[0]/1)*1;
	u[1]=u[1]-floor(u[1]/1)*1;**/
	/**if(u[0]<0){u[0]=0;}
	if(u[1]<0){u[1]=0;}
	if(u[0]>1){u[0]=1;}
	if(u[1]>1){u[1]=1;}*/
	
}

void Menger_Sponge_4DSlice::get_Bump_Perturbation(double x[], vector<Texture*> textures, double dn[]){
	double n[3]={0,0,0};
	
	this->get_Normal(x, n);
	double X[3];X[0]=x[0]-c[0];X[1]=x[1]-c[1];X[2]=x[2]-c[2];
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
	
	if(d_min == d12 || d_min == d12m){
		u[0]=(X[0]*e12[0]+X[1]*e12[1]+X[2]*e12[2]-(-2))/4;
		u[1]=(X[0]*u12[0]+X[1]*u12[1]+X[2]*u12[2]-(-2))/4;
	}
	else if(d_min == d23 || d_min == d23m){
		u[0]=(X[0]*e23[0]+X[1]*e23[1]+X[2]*e23[2]-(-2))/4;
		u[1]=(X[0]*u23[0]+X[1]*u23[1]+X[2]*u23[2]-(-2))/4;
	}
	else if(d_min == d31 || d_min == d31m){
		u[0]=(X[0]*e31[0]+X[1]*e31[1]+X[2]*e31[2]-(-2))/4;
		u[1]=(X[0]*u31[0]+X[1]*u31[1]+X[2]*u31[2]-(-2))/4;
	}
	else if(d_min == d4 || d_min == d4m){
		u[0]=(X[0]-(-2))/4.0;
		u[1]=(X[1]-(-2))/4.0;
	}
	u[0]=u[0]-floor(u[0]/1)*1;
	u[1]=u[1]-floor(u[1]/1)*1;
	
	int col[3];
	
	if(this->bump_ID>-1 && this->bump_ID<textures.size()){
		int w=textures[bump_ID]->get_width()-1;
		int h=textures[bump_ID]->get_height()-1;
		int i=int(u[0]*w);int j=int((1-u[1])*h);
		//textures[bump_ID]->get_at(i,j,col);
		textures[bump_ID]->bilinear_get_at(u, col);
	}
	else{
		col[0]=color[0];
		col[1]=color[1];
		col[2]=color[2];
	}
	if(d_min == d12m){
		dn[0]=e12[0]*double((col[0]-127)/127.0)+u12[0]*double((col[1]-127)/127.0);
		dn[1]=e12[1]*double((col[0]-127)/127.0)+u12[1]*double((col[1]-127)/127.0);
		dn[2]=e12[2]*double((col[0]-127)/127.0)+u12[2]*double((col[1]-127)/127.0);
	}
	else if(d_min == d12){
		dn[0]=-e12[0]*double((col[0]-127)/127.0)-u12[0]*double((col[1]-127)/127.0);
		dn[1]=-e12[1]*double((col[0]-127)/127.0)-u12[1]*double((col[1]-127)/127.0);
		dn[2]=-e12[2]*double((col[0]-127)/127.0)-u12[2]*double((col[1]-127)/127.0);
	}
	else if(d_min == d23m){
		dn[0]=e23[0]*double((col[0]-127)/127.0)+u23[0]*double((col[1]-127)/127.0);
		dn[1]=e23[1]*double((col[0]-127)/127.0)+u23[1]*double((col[1]-127)/127.0);
		dn[2]=e23[2]*double((col[0]-127)/127.0)+u23[2]*double((col[1]-127)/127.0);
	}
	else if(d_min == d23){
		dn[0]=-e23[0]*double((col[0]-127)/127.0)-u23[0]*double((col[1]-127)/127.0);
		dn[1]=-e23[1]*double((col[0]-127)/127.0)-u23[1]*double((col[1]-127)/127.0);
		dn[2]=-e23[2]*double((col[0]-127)/127.0)-u23[2]*double((col[1]-127)/127.0);
	}
	else if(d_min == d31m){
		dn[0]=e31[0]*double((col[0]-127)/127.0)+u31[0]*double((col[1]-127)/127.0);
		dn[1]=e31[1]*double((col[0]-127)/127.0)+u31[1]*double((col[1]-127)/127.0);
		dn[2]=e31[2]*double((col[0]-127)/127.0)+u31[2]*double((col[1]-127)/127.0);
	}
	else if(d_min == d31){
		dn[0]=-e31[0]*double((col[0]-127)/127.0)-u31[0]*double((col[1]-127)/127.0);
		dn[1]=-e31[1]*double((col[0]-127)/127.0)-u31[1]*double((col[1]-127)/127.0);
		dn[2]=-e31[2]*double((col[0]-127)/127.0)-u31[2]*double((col[1]-127)/127.0);
	}
	else if(d_min == d4m){
		dn[0]=double((col[0]-127)/127.0);
		dn[1]=double((col[1]-127)/127.0);
		dn[2]=0;
	}
	else if(d_min == d4){
		dn[0]=-double((col[0]-127)/127.0);
		dn[1]=-double((col[1]-127)/127.0);
		dn[2]=0;
	}

	
/**	double d1 = (n[0]-1)*(n[0]-1)+n[1]*n[1]+n[2]*n[2];
	double d2 = (n[0]+1)*(n[0]+1)+n[1]*n[1]+n[2]*n[2];
	double d3 = n[0]*n[0]+(n[1]-1)*(n[1]-1)+n[2]*n[2];
	double d4 = n[0]*n[0]+(n[1]+1)*(n[1]+1)+n[2]*n[2];
	double d5 = n[0]*n[0]+n[1]*n[1]+(n[2]-1)*(n[2]-1);
	double d6 = n[0]*n[0]+n[1]*n[1]+(n[2]+1)*(n[2]+1);
	if(min(d1,d2)<min(min(d3,d4),min(d5,d6))){
		u[0]=(X[1]-(-2))/4.0;
		u[1]=(X[2]-(-2))/4.0;
	}
	else if(min(d3,d4)<min(min(d5,d6),min(d1,d2))){
		u[0]=(X[0]-(-2))/4.0;
		u[1]=(X[2]-(-2))/4.0;		
	}
	else if(min(d5,d6)<min(min(d1,d2),min(d3,d4))){
		u[0]=(X[0]-(-2))/4.0;
		u[1]=(X[1]-(-2))/4.0;			
	}
	u[0]=u[0]-floor(u[0]/1)*1;
	u[1]=u[1]-floor(u[1]/1)*1;
	
	int col[3];
	
	if(this->bump_ID>-1 && this->bump_ID<textures.size()){
		int w=textures[bump_ID]->get_width()-1;
		int h=textures[bump_ID]->get_height()-1;
		int i=int(u[0]*w);int j=int((1-u[1])*h);
		//textures[bump_ID]->get_at(i,j,col);
		textures[bump_ID]->bilinear_get_at(u, col);
	}
	else{
		col[0]=color[0];
		col[1]=color[1];
		col[2]=color[2];
	}
	if(min(d1,d2)<min(min(d3,d4),min(d5,d6))){
		dn[0]=0;
		dn[1]=(col[0]-127)/127.0;
		dn[2]=(col[1]-127)/127.0;
	}
	else if(min(d3,d4)<min(min(d5,d6),min(d1,d2))){
		dn[0]=(col[0]-127)/127.0;
		dn[1]=0;
		dn[2]=(col[1]-127)/127.0;		
	}
	else if(min(d5,d6)<min(min(d1,d2),min(d3,d4))){
		dn[0]=(col[0]-127)/127.0;
		dn[1]=(col[1]-127)/127.0;
		dn[2]=0;		
	}**/
	
}
