#include "Shape.h"
#include "invisibility.h"
#include <vector>


class Light{
public:
	Light(double,int[]);
	virtual void get_Contribution(double[],double[],double v[],vector<Texture*> textures,double,double[])=0;
	virtual void get_Specular_Contribution(double r[],double x[],double specular_term[])=0;
	virtual bool is_Visible(double[],double,double[],std::vector<Shape*>)=0;
	virtual double percent_Visible(double x[],double R,vector<Shape*> shapes)=0;
	int illumination_map_ID;
	double pp;
	bool is_flat;
protected:
	int color[3];
	double brightness;
};


class ambientLight : public Light{
public:
	ambientLight(double,int[]);
	void get_Contribution(double[],double[],double v[],vector<Texture*> textures,double,double[]);
	void get_Specular_Contribution(double r[],double x[],double specular_term[]);
	bool is_Visible(double[],double,double[],std::vector<Shape*>);
	double percent_Visible(double x[],double R,vector<Shape*> shapes);
};

class pointLight : public Light{
public:
	pointLight(double[],double,int[]);
	void color_from_illumination_map(double x[],double n[],vector<Texture*> textures, int C[]);
	void get_Contribution(double[],double[],double v[],vector<Texture*> textures,double,double[]);
	void get_Specular_Contribution(double r[],double x[],double specular_term[]);
	bool is_Visible(double[],double,double[],std::vector<Shape*>);
	double percent_Visible(double x[],double R,vector<Shape*> shapes);
private:
	double pos[3];
};

class planeLight : public Light{
public:
	planeLight(double Z, double b, int C[],Texture*);
	void get_Contribution(double[],double[],double v[],vector<Texture*> textures,double,double[]);
	void getColor(double x, double y, int color[3]);
	void get_Specular_Contribution(double r[],double x[],double specular_term[]);
	bool is_Visible(double[],double,double[],std::vector<Shape*>);
	double percent_Visible(double x[],double R,vector<Shape*> shapes);
private:
	double z;
	Texture* texture;
	bool use_texture = false;
};

