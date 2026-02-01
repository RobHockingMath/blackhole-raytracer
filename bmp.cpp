// Utility for saving a bitmap to a file using a callback function.
//

#include "bmp.h"
#include <stdlib.h>
#include <stdio.h>
#include "math.h"
#include <iostream>

using namespace std;

int save_bmp(long width, long height, char* filename, int ***img) {
	printf("huh.");
	COLOR color;
	char *ptr;
	char id[2];
	BMPHEAD bh;

	memset ((char *)&bh,0,sizeof(BMPHEAD)); /* sets entire header to 0 */
	id[0] = 'B';
	id[1] = 'M';
	/* bh.filesize  = calculated size of your file (see below)*/
	/* bh.reserved  = /* two zero bytes */
	bh.headersize  = 54L;  /* (for 24 bit images) */
	bh.infoSize  =  0x28L;  /*(for 24 bit images) */
	bh.width     = width; /* width in pixels of your image */
	bh.height     = height; /* depth in pixels of your image */
	bh.biPlanes  =  1; /* (for 24 bit images) */
	bh.bits      = 24; /* (for 24 bit images) */
	bh.biCompression = 0L; /* (no compression) */

	int bytesPerLine;

	bytesPerLine = bh.width * 3;  /* (for 24 bit images) */
	/* round up to a dword boundary */
	if (bytesPerLine & 0x0003) 
	{
	bytesPerLine |= 0x0003;
	++bytesPerLine;
	}

	bh.filesize = bh.headersize+(long)bytesPerLine*bh.height;
	FILE * bmpfile;
    FILE * fp;
    int i = 0;
    
	bmpfile = fopen(filename, "wb");

	if (bmpfile == NULL)
	{
		printf("Error opening output file\n");
		/* -- close all open files and free any allocated memory -- */
		return 1;
	}

	fwrite(id, 2, sizeof (char), bmpfile);
	fwrite(&bh, 1, sizeof (bh), bmpfile);

	char linebuf[30000];

	//linebuf = (char *) calloc(1, bytesPerLine);

	if (linebuf == NULL)
	{
		printf ("Error allocating memory\n");
		/* -- close all open files and free any allocated memory -- */
		fclose(bmpfile);
		return(1);
	}

	for (int line = bh.height - 1; line >= 0; line --)  {
		/* fill line linebuf with the image data for that line */
		//for (ptr = linebuf; ptr <= linebuf+bytesPerLine-3; ptr+=3) {
		ptr = linebuf;
		for (int col = 0; col < bh.width; col++) {
			
			///callback(line,col,&color);
            color.red=img[col][line][0];
            color.green=img[col][line][1];
            color.blue=img[col][line][2];
      
			// write colors for pixel to bitmap memory.
			*ptr = (int)(color.blue);
			*(ptr+1) = (int)(color.green);
			*(ptr+2) = (int)(color.red);
			ptr += 3;
		}
		fwrite(linebuf, 1, bytesPerLine, bmpfile);
        fflush(bmpfile);
	}

	fclose(bmpfile);
	return 0;
}
