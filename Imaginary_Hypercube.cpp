#include "Imaginary_Hypercube.h"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>

using namespace std;

Imaginary_Hypercube::Imaginary_Hypercube(double N[],bool is_T_,int r,int g,int b)
	: BasicShape(r,g,b)
{
	
	is_T = is_T_;
	n[0]=N[0];n[1]=N[1];n[2]=N[2];n[3]=N[3];
	double len = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]+n[3]*n[3]);
	n[0]/=len;n[1]/=len;n[2]/=len;n[3]/=len;
	
	// let's build an orthonormal basis for the hyperplane orthogonal to n
	
	Eigen::MatrixXd A(1, 4);
    A(0,0)=n[0];A(0,1)=n[1];A(0,2)=n[2];A(0,3)=n[3];
   
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
	
	
	std::ifstream inputFile("V_T.txt",std::ifstream::in);
	if (inputFile.is_open()) {
		std::cout<<"File read successful!"<<std::endl;
		for (int i = 0; i < 17; ++i){
			double projection;			
			for (int j = 0; j < 4; ++j){
				inputFile >> Vertices_T[i][j];
			}
			projection = Vertices_T[i][0]*n[0]+Vertices_T[i][1]*n[1]+Vertices_T[i][2]*n[2]+Vertices_T[i][3]*n[3];
			X_T3.push_back( (Vertices_T[i][0]-projection*n[0])*e1[0]+(Vertices_T[i][1]-projection*n[1])*e1[1]+(Vertices_T[i][2]-projection*n[2])*e1[2]+(Vertices_T[i][3]-projection*n[3])*e1[3] );
			Y_T3.push_back( (Vertices_T[i][0]-projection*n[0])*e2[0]+(Vertices_T[i][1]-projection*n[1])*e2[1]+(Vertices_T[i][2]-projection*n[2])*e2[2]+(Vertices_T[i][3]-projection*n[3])*e2[3] );
			Z_T3.push_back( (Vertices_T[i][0]-projection*n[0])*e3[0]+(Vertices_T[i][1]-projection*n[1])*e3[1]+(Vertices_T[i][2]-projection*n[2])*e3[2]+(Vertices_T[i][3]-projection*n[3])*e3[3] );
			
		}
	}
	else{
		std::cerr << "Unable to open file." << std::endl;
	}
	inputFile.close();
	
	std::ifstream inputFile2("T_T.txt",std::ifstream::in);
	if (inputFile2.is_open()) {
		std::cout<<"File read successful!"<<std::endl;
		for (int i = 0; i < 27; ++i){
			double projection;
			for (int j = 0; j < 4; ++j){
				inputFile2 >> Translations_T[i][j];
			}
			projection = Translations_T[i][0]*n[0]+Translations_T[i][1]*n[1]+Translations_T[i][2]*n[2]+Translations_T[i][3]*n[3];
			Translations_T3[i][0] = (Translations_T[i][0]-projection*n[0])*e1[0]+(Translations_T[i][1]-projection*n[1])*e1[1]+(Translations_T[i][2]-projection*n[2])*e1[2]+(Translations_T[i][3]-projection*n[3])*e1[3];
			Translations_T3[i][1] = (Translations_T[i][0]-projection*n[0])*e2[0]+(Translations_T[i][1]-projection*n[1])*e2[1]+(Translations_T[i][2]-projection*n[2])*e2[2]+(Translations_T[i][3]-projection*n[3])*e2[3];
			Translations_T3[i][2] = (Translations_T[i][0]-projection*n[0])*e3[0]+(Translations_T[i][1]-projection*n[1])*e3[1]+(Translations_T[i][2]-projection*n[2])*e3[2]+(Translations_T[i][3]-projection*n[3])*e3[3];
		}
	}
	else{
		std::cerr << "Unable to open file." << std::endl;
	}
	inputFile2.close();
	
	std::ifstream inputFile3("V_H.txt",std::ifstream::in);
	if (inputFile3.is_open()) {
		std::cout<<"File read successful!"<<std::endl;
		for (int i = 0; i < 14; ++i){
			double projection;
			for (int j = 0; j < 4; ++j){
				inputFile3 >> Vertices_H[i][j];
			}
			projection = Vertices_H[i][0]*n[0]+Vertices_H[i][1]*n[1]+Vertices_H[i][2]*n[2]+Vertices_H[i][3]*n[3];
			X_H3.push_back( (Vertices_H[i][0]-projection*n[0])*e1[0]+(Vertices_H[i][1]-projection*n[1])*e1[1]+(Vertices_H[i][2]-projection*n[2])*e1[2]+(Vertices_H[i][3]-projection*n[3])*e1[3] );
			Y_H3.push_back( (Vertices_H[i][0]-projection*n[0])*e2[0]+(Vertices_H[i][1]-projection*n[1])*e2[1]+(Vertices_H[i][2]-projection*n[2])*e2[2]+(Vertices_H[i][3]-projection*n[3])*e2[3] );
			Z_H3.push_back( (Vertices_H[i][0]-projection*n[0])*e3[0]+(Vertices_H[i][1]-projection*n[1])*e3[1]+(Vertices_H[i][2]-projection*n[2])*e3[2]+(Vertices_H[i][3]-projection*n[3])*e3[3] );
		}
	}
	else{
		std::cerr << "Unable to open file." << std::endl;
	}
	inputFile3.close();
	
	std::ifstream inputFile4("T_H.txt",std::ifstream::in);
	if (inputFile4.is_open()) {
		std::cout<<"File read successful!"<<std::endl;
		for (int i = 0; i < 27; ++i){
			double projection;
			for (int j = 0; j < 4; ++j){
				inputFile4 >> Translations_H[i][j];
			}
			projection = Translations_H[i][0]*n[0]+Translations_H[i][1]*n[1]+Translations_H[i][2]*n[2]+Translations_H[i][3]*n[3];
			Translations_H3[i][0] = (Translations_H[i][0]-projection*n[0])*e1[0]+(Translations_H[i][1]-projection*n[1])*e1[1]+(Translations_H[i][2]-projection*n[2])*e1[2]+(Translations_H[i][3]-projection*n[3])*e1[3];
			Translations_H3[i][1] = (Translations_H[i][0]-projection*n[0])*e2[0]+(Translations_H[i][1]-projection*n[1])*e2[1]+(Translations_H[i][2]-projection*n[2])*e2[2]+(Translations_H[i][3]-projection*n[3])*e2[3];
			Translations_H3[i][2] = (Translations_H[i][0]-projection*n[0])*e3[0]+(Translations_H[i][1]-projection*n[1])*e3[1]+(Translations_H[i][2]-projection*n[2])*e3[2]+(Translations_H[i][3]-projection*n[3])*e3[3];
		}
	}
	else{
		std::cerr << "Unable to open file." << std::endl;
	}
	inputFile4.close();
	for(int i=0;i<X_H3.size();i++){
		std::cout<<"i = "<<i<<"X_H3= "<<X_H3[i]<<" "<<"Y_H3= "<<Y_H3[i]<<" "<<"Z_H3= "<<Z_H3[i]<<std::endl;
	}
	for(int i=0;i<X_T3.size();i++){
		std::cout<<"i = "<<i<<"X_T3= "<<X_T3[i]<<" "<<"Y_T3= "<<Y_T3[i]<<" "<<"Z_T3= "<<Z_T3[i]<<std::endl;
	}
	double s = 3.0;
	if(is_T){
		ConvexHull *T3 = new ConvexHull(X_T3,Y_T3,Z_T3,r,g,b);
		for(int i=0;i<(T3->parts).size();i++){
			Triangle* T = new Triangle(((Triangle*)(T3->parts[i]))->p1,((Triangle*)(T3->parts[i]))->p2,((Triangle*)(T3->parts[i]))->p3,((Triangle*)(T3->parts[i]))->n,r,g,b);
			Triangles.push_back(T);
		}
		free(T3);
	}
	else{
		ConvexHull *H3 = new ConvexHull(X_H3,Y_H3,Z_H3,r,g,b);
		for(int i=0;i<(H3->parts).size();i++){
			Triangle* T = new Triangle(((Triangle*)(H3->parts[i]))->p1,((Triangle*)(H3->parts[i]))->p2,((Triangle*)(H3->parts[i]))->p3,((Triangle*)(H3->parts[i]))->n,r,g,b);
			Triangles.push_back(T);	
		}
		free(H3);
	}
	std::cout<<"Triangles.size()== "<<Triangles.size();
	
	
	
	is_wireframe=false;
	
}

double Imaginary_Hypercube::get_Intersection(double p[],double v[],double x[]){
	return get_Intersection_depth(p,v,x,recursion_depth);
}

double Imaginary_Hypercube::get_Intersection_depth(double p[],double v[],double x[],int depth){
	
	//if(depth == 1){
	//	std::cout<<"depth = "<<depth<<std::endl;
	//}
	//if(depth == 0){return -1;}
	
	double t=-1;
	double t_i;
	double y[4];
	for(int i=0;i<Triangles.size();i++){
		t_i=Triangles[i]->get_Intersection(p,v,x);
		if(t_i>0.01 && (t_i<t || t<0) ){
			t=t_i;
			y[0]=x[0];y[1]=x[1];y[2]=x[2];y[3]=double(i);
		}
	}
	if(t == -1){return t;}
	else if(depth==1){x[0]=y[0];x[1]=y[1];x[2]=y[2];x[3]=y[3];return t;}
	else{
		std::vector<double> t_s;
		std::vector<double> x_s;
		std::vector<double> y_s;
		std::vector<double> z_s;
		std::vector<double> w_s;
		double dx,dy,dz;
		for(int ii=0;ii<27;ii++){
			if(is_T){  
				dx=2*Translations_T3[ii][0];
				dy=2*Translations_T3[ii][1];
				dz=2*Translations_T3[ii][2];
			}
			else{
				dx=2*Translations_H3[ii][0];
				dy=2*Translations_H3[ii][1];
				dz=2*Translations_H3[ii][2];
			}
			//for(int i=0;i<Triangles.size();i++){
			//	Triangles[i]->scaled_translation(dx,dy,dz,s);
			//}
			double pp[3]={s*p[0]-dx,s*p[1]-dy,s*p[2]-dz};
			double vv[3]={s*v[0],s*v[1],s*v[2]};
			
			//t_s.push_back(this->get_Intersection_depth(p,v,x,depth-1));
			t_s.push_back(this->get_Intersection_depth(pp,vv,x,depth-1));
			//x_s.push_back(x[0]);
			//y_s.push_back(x[1]);
			//z_s.push_back(x[2]);
			//w_s.push_back(x[3]);
			x_s.push_back(p[0]+t_s[ii]*v[0]);
			y_s.push_back(p[1]+t_s[ii]*v[1]);
			z_s.push_back(p[2]+t_s[ii]*v[2]);
			w_s.push_back(x[3]);
			//for(int i=0;i<Triangles.size();i++){
			//	Triangles[i]->scaled_translation(-dx/s,-dy/s,-dz/s,1.0/s);
			//}
		}
		double t_min=-1;
		for(int ii=0;ii<27;ii++){
			if(t_s[ii]!=-1){
				//std::cout<<"t_s[ii]=="<<t_s[ii]<<std::endl;
				if( t_min==-1 || (t_min!=-1 && t_s[ii]<t_min) ){
					t_min=t_s[ii];
					x[0]=x_s[ii];
					x[1]=y_s[ii];
					x[2]=z_s[ii];
					x[3]=w_s[ii];
				}
			}
		}
		return t_min;
	}	
}

bool Imaginary_Hypercube::wireframe_condition(double p[]){
	return false;
}

//4th component of x tells which part is relevant
void Imaginary_Hypercube::get_Texture_Coords(double x[], double u[]){
	Triangles[int(x[3])]->get_Texture_Coords(x,u);
}

//4th component of x tells which part is relevant
void Imaginary_Hypercube::get_Normal(double x[], double n[]){
	//std::cout<<" part number "<<int(x[3])<<std::endl;
	Triangles[int(x[3])]->get_Normal(x,n);
}

void Imaginary_Hypercube::set_Texture_IDs(int IDs[]){
	for(int i=0;i<Triangles.size();i++){
		Triangles[i]->texture_ID=IDs[i];
	}
}

void Imaginary_Hypercube::get_Color(double x[], vector<Texture*> textures, int col[]){
	Triangles[int(x[3])]->get_Color(x,textures,col);
}



