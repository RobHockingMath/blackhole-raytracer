#ifndef TEXTURE_H
#define TEXTURE_H

#include <string>
#include <cstring>
#include <cstdint>

class Texture {
public:
	Texture(char*);
	~Texture();
	void load(char*);
	//void save(int);
	void get_at(int,int,int[]);
	void bilinear_get_at(double u[], int color[]);
	void scale_data(double s_r, double s_g, double s_b);
	int get_width();
	int get_height();
private:
	int w,h;
	uint8_t ***data;
	std::string name;
};

#endif
