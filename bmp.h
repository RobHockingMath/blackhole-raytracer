// Utility for saving a bitmap to a file using a callback function.
//

#include <string.h>

typedef struct {
	double red;
	double green;
	double blue;
} COLOR;

typedef struct {
	long filesize;
	char reserved[2];
	long headersize;
	long infoSize;
	long width;
	long height;
	short biPlanes;
	short bits;
	long biCompression;
	long biSizeImage;
	long biXPelsPerMeter;
	long biYPelsPerMeter;
	long biClrUsed;
	long biClrImportant;
} BMPHEAD;

/* saves a bitmap file */
int save_bmp(long width, long depth, char* filename, int ***img);
