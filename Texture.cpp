#include <iostream>
#include <stdio.h>
#include "Texture.h"
#include "math.h"
#include "bmp.h"

using namespace std;

Texture::Texture(char* path){
	this->load(path);
}

Texture::~Texture(){
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){
			delete[] data[i][j];
		}	
		delete[] data[i];
	}
}


void Texture::load(char* path){
	FILE *myfile;
	myfile=fopen(path,"r");
	fscanf(myfile,"%d %d",&w,&h);
	data = new uint8_t**[w];
	for(int i=0;i<w;i++){
		data[i]=new uint8_t*[h];
		for(int j=0;j<h;j++){
			data[i][j]=new uint8_t[3];
		}
	}
	for(int i=0;i<w;i++){
		//cout <<"loading texture "<<path<<" "<<i<<"/"<<w<<endl;
		for(int j=0;j<h;j++){
			int temp0, temp1, temp2;  // Temporary int variables
			fscanf(myfile,"%d %d %d",&temp0,&temp1,&temp2);
			data[i][j][0] = static_cast<uint8_t>(temp0);
			data[i][j][1] = static_cast<uint8_t>(temp1);
			data[i][j][2] = static_cast<uint8_t>(temp2);
		}
	}
	fclose(myfile);
	
	//char *filename_i = new char[filename_str.length() + 1];
	//strcpy(filename_i, filename_str.c_str());
	
	name = path;
}

/**void Texture::save(int thread_id){
	
    string string_extension = "_texture.bmp";
	string temporary_name = name+to_string(thread_id)+string_extension;
	char* save_name = new char[temporary_name.length() + 1];
	strcpy(save_name, temporary_name.c_str());
    
	save_bmp(w,h,save_name,data);	
}**/

void Texture::get_at(int i,int j, int color[3]){
	if(i<0 || j<0 || i>w-1 || j>h-1){
		cout<<"texture read out of bounds: i,j = "<<i<<","<<j<<" but w,h = "<<w<<","<<h<<endl;
	}
	else{
		color[0]=static_cast<int>(data[i][j][0]);
		color[1]=static_cast<int>(data[i][j][1]);
		color[2]=static_cast<int>(data[i][j][2]);
	}
}

void Texture::scale_data(double s_r, double s_g, double s_b){
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){
			data[i][j][0]=static_cast<uint8_t>(floor(s_r*static_cast<int>(data[i][j][0])));
			data[i][j][1]=static_cast<uint8_t>(floor(s_g*static_cast<int>(data[i][j][1])));
			data[i][j][2]=static_cast<uint8_t>(floor(s_b*static_cast<int>(data[i][j][2])));
		}
	}
}

void Texture::bilinear_get_at(double u[], int color[]){
	
	
	
	int ii=floor(u[0]*double(w-1));
	int ip=min(ii+1,w-1);
	int jj=floor(u[1]*double(h-1));
	int jp=min(jj+1,h-1);
	
	int C00[3]={0,0,0};
	int C01[3]={0,0,0};
	int C10[3]={0,0,0};
	int C11[3]={0,0,0};
	
	//if(color[2]==-2){
	//	cout<<"w = "<<w<<" h = "<<h<<endl;
	//	cout<<"ii = "<<ii<<" jj = "<<jj<<endl;
	//}
	
	this->get_at(ii,jj,C00);
	this->get_at(ip,jj,C10);
	this->get_at(ii,jp,C01);
	this->get_at(ip,jp,C11);
	
	double s = u[0]*double(w-1)-ii;
	double t = u[1]*double(h-1)-jj;
	
	//if(color[2]==-2){
	//	cout<<"s = "<<s<<" t = "<<t<<endl;
	//}
	
	color[0]=int((1-s)*(1-t)*C00[0]+s*(1-t)*C10[0]+(1-s)*t*C01[0]+s*t*C11[0]);
	color[1]=int((1-s)*(1-t)*C00[1]+s*(1-t)*C10[1]+(1-s)*t*C01[1]+s*t*C11[1]);
	color[2]=int((1-s)*(1-t)*C00[2]+s*(1-t)*C10[2]+(1-s)*t*C01[2]+s*t*C11[2]);
}

int Texture::get_width(){
	return w;
}

int Texture::get_height(){
	return h;
}
