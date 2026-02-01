#include <math.h>
#include <iostream>
#include "Stellated_Menger_Sponge_4DSlice.h"
#include "Sphere.h"
#include "invisibility.h"

#define PI 3.14159
#define SQRT3 1.7321

using namespace std;

Stellated_Menger_Sponge_4DSlice::Stellated_Menger_Sponge_4DSlice(double c1,double c2,double c3,double scale_,int depth_,int type_,double E,double offset_,int red,int green, int blue)
	: Distance_fractal(c1,c2,c3,E,red,green,blue)
{
   depth = depth_;
   scale = scale_;
   offset = offset_;
   type = type_;
}


void Stellated_Menger_Sponge_4DSlice::rot(double p[]){
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


bool Stellated_Menger_Sponge_4DSlice::wireframe_condition(double p[]){

	double X=p[0]-(this->c[0]), Y=p[1]-(this->c[1]), Z=p[2]-(this->c[2]);
	
	//double X=p[0], Y=p[1], Z=p[2];
	
	X=X/scale;Y=Y/scale;Z=Z/scale; //center it by changing position and scale	
	
	double x = 0.25*offset+(X+Y-Z)/2.0;
	double y = 0.25*offset+(-X-Y-Z)/2.0;
	double z = 0.25*offset+(X-Y+Z)/2.0;
	double w = 0.25*offset+(-X+Y+Z)/2.0;
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
double Stellated_Menger_Sponge_4DSlice::Distance(double p[]){

	//this->rot(p);

	double X=p[0], Y=p[1], Z=p[2];
	
	X=X/scale;Y=Y/scale;Z=Z/scale; //center it by changing position and scale	
	
	double x = 0.25*offset+(X+Y-Z)/2.0;
	double y = 0.25*offset+(-X-Y-Z)/2.0;
	double z = 0.25*offset+(X-Y+Z)/2.0;
	double w = 0.25*offset+(-X+Y+Z)/2.0;
	
	/**double x = (X+Y-Z)/2.0;
	double y = (-X-Y-Z)/2.0;
	double z = (X-Y+Z)/2.0;
	double w = offset+(-X+Y+Z)/2.0;**/
	
	double offset_p = offset/2+2;
	
	x = (1+x)/2;y = (1+y)/2; z = (1+z)/2;w=(1+w)/2;

	bool x_lt_m1 = false;
	bool y_lt_m1 = false;
	bool z_lt_m1 = false;
	bool w_lt_m1 = false;
	if(x<-1 || x>2){x_lt_m1=true;}
	if(y<-1 || y>2){y_lt_m1=true;}
	if(z<-1 || z>2){z_lt_m1=true;}
	if(w<-1 || w>2){w_lt_m1=true;}
	

	double xx=abs(x-0.5)-0.5, yy=abs(y-0.5)-0.5, zz=abs(z-0.5)-0.5, ww=abs(w-0.5)-0.5;//ww=offset_p-x-y-z;//ww=abs(offset_p-x-y-z-0.5)-0.5;//ww=abs(w-0.5)-0.5;
	
	double d1;
	//d1=max(xx,max(yy,max(zz,ww))); //distance to the box
	//double dxyz = max(xx,max(yy,zz)); // w omitted
	//double dxyw = max(xx,max(yy,ww)); // z omitted
	//double dxzw = max(xx,max(zz,ww)); // y omitted
	//double dyzw = max(yy,max(zz,ww)); // x omitted
	//d1 = min(dxyz,min(dxyw,min(dxzw,dyzw)));
	d1=min(max(xx,zz),min(max(xx,yy),max(yy,zz))); //w ommitted
	d1 = max(d1,min(max(xx,ww),min(max(xx,yy),max(yy,ww)))); // z ommitted
	d1 = max(d1,min(max(xx,zz),min(max(xx,ww),max(ww,zz)))); // y ommitted 
	d1 = max(d1,min(max(ww,zz),min(max(ww,yy),max(yy,zz)))); // x ommitted
	
	double xxx=abs(x-0.5)-1.5, yyy=abs(y-0.5)-1.5, zzz=abs(z-0.5)-1.5, www=abs(w-0.5)-1.5;
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
	bool x_neg=false;
	bool y_neg=false;
	bool z_neg=false;
	bool w_neg=false;
	
	for (int i=1; i<=depth; ++i) {
		
		
		
		//fmod fucks shit up for negative numbers.
		double x1=x,y1=y,z1=z,w1=w;
		
		if(i>1){
		while(x1<0){x1=1+x1;x_neg=true;}
		while(y1<0){y1=1+y1;y_neg=true;}
		while(z1<0){z1=1+z1;z_neg=true;}
		while(w1<0){w1=1+w1;w_neg=true;}
		}
		
		double xa = fmod(3.0*x1*pp,3.0);
		double ya = fmod(3.0*y1*pp,3.0);
		double za = fmod(3.0*z1*pp,3.0);
		//double wa = fmod(3.0*w*pp,3.0);
		double wa = fmod(3.0*w1*pp,3.0);
		pp*=3.0;

		//we can also translate/rotate (xa,ya,za) without affecting the DE estimate

		double xx=0.5-abs(xa-1.5), yy=0.5-abs(ya-1.5), zz=0.5-abs(za-1.5), ww=0.5-abs(wa-1.5);
		
		if(i==1){
			if(x1>=1){xx=-100;}
			if(y1>=1){yy=-100;}
			if(z1>=1){zz=-100;}
			if(w1>=1){ww=-100;}
		}
		
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
	
	d = max(d,max(xxx,max(yyy,max(zzz,www))));
	
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


double Stellated_Menger_Sponge_4DSlice::is_Outside(double p[]){
	
	double X=p[0], Y=p[1], Z=p[2];
	
	X=X/scale;Y=Y/scale;Z=Z/scale; //center it by changing position and scale	
	
	double x = 0.25*offset+(X+Y-Z)/2.0;
	double y = 0.25*offset+(-X-Y-Z)/2.0;
	double z = 0.25*offset+(X-Y+Z)/2.0;
	double w = 0.25*offset+(-X+Y+Z)/2.0;
	
	
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
