#pragma once
#include "stdafx.h"
#include <string>
#include <iostream>

inline int ifloor(real x) {
	return (x > 0) ? (int)x : (int)x - 1;
}

inline int loopIndex(int pos, int dim) {
	if (pos >= 0) {
		return pos % dim;
	}
	else {
		return pos % dim + dim;
	}
}

inline real fieldMax(real *f, int l, bool checkZero = false) {
	real max = 0;
	bool zero = false, allZero = true;
	for (int i = 0; i < l; i++) {
		if (f[i] != f[i]) {
			std::cout << "INF" << std::endl; 
			exit(0);
		}
		if (f[i] == 0) {
			zero = true;
		}
		else {
			allZero = false;
		}
		if (f[i] * f[i] > max) {
			max = f[i] * f[i];
		}
	}
	if (zero && checkZero) {
		std::cout << "zero value occurred." << std::endl;
	}
	if (allZero && checkZero) {
		std::cout << "all zero field."<< std::endl;
	}
	return max;
}

inline int getIndex(int x, int y, int z, int dimX, int dimY, int dimZ) {
	return x + y * dimX + z * dimX * dimY;
}

inline void getPos(int index, int dimX, int dimY, int dimZ, int &i, int &j, int &k) {
	i = index % dimX;
	j = (index % (dimX * dimY)) / dimX;
	k = index / (dimX * dimY);
}

inline void setBorderValue(real *field, int dimX, int dimY, int dimZ, real value) {
	for (int z = 0; z < dimZ; z++) {
		for (int y = 0; y < dimY; y++) {
			field[getIndex(0, y, z, dimX, dimY, dimZ)] = value;
			field[getIndex(dimX - 1, y, z, dimX, dimY, dimZ)] = value;
		}
	}
	for (int x = 0; x < dimX; x++) {
		for (int z = 0; z < dimZ; z++) {
			field[getIndex(x, 0, z, dimX, dimY, dimZ)] = value;
			field[getIndex(x, dimY - 1, z, dimX, dimY, dimZ)] = value;
		}
	}
	for (int x = 0; x < dimX; x++) {
		for (int y = 0; y < dimY; y++) {
			field[getIndex(x, y, 0, dimX, dimY, dimZ)] = value;
			field[getIndex(x, y, dimZ-1, dimX, dimY, dimZ)] = value;
		}
	}
}

inline void setBorderType(CELL_TYPE *field, int dimX, int dimY, int dimZ, CELL_TYPE value) {
	for (int z = 0; z < dimZ; z++) {
		for (int y = 0; y < dimY; y++) {
			field[getIndex(0, y, z, dimX, dimY, dimZ)] = value;
			field[getIndex(dimX - 1, y, z, dimX, dimY, dimZ)] = value;
		}
	}
	for (int x = 0; x < dimX; x++) {
		for (int z = 0; z < dimZ; z++) {
			field[getIndex(x, 0, z, dimX, dimY, dimZ)] = value;
			field[getIndex(x, dimY - 1, z, dimX, dimY, dimZ)] = value;
		}
	}
	for (int x = 0; x < dimX; x++) {
		for (int y = 0; y < dimY; y++) {
			field[getIndex(x, y, 0, dimX, dimY, dimZ)] = value;
			field[getIndex(x, y, dimZ - 1, dimX, dimY, dimZ)] = value;
		}
	}
}

inline real max(real a, real b) {
	return a > b ? a : b;
}
inline real max(real a, real b, real c, real d) {
	return max(max(a, b), max(c, d));
}
inline real max(real a, real b, real c, real d, real e, real f, real g, real h) {
	return max(max(a, b, c, d), max(e, f, g, h));
}

inline real min(real a, real b) {
	return a < b ? a : b;
}
inline real min(real a, real b, real c, real d) {
	return min(min(a, b), min(c, d));
}
inline real min(real a, real b, real c, real d, real e, real f, real g, real h) {
	return min(min(a, b, c, d), min(e, f, g, h));
}
//
//inline real debug(real *a, int l) {
//	real sum = 0;
//	for (int i = 0; i < l; i++) {
//		sum += a[i] * i;
//	}
//	return sum;
//}
#include "zlib/zlib.h"
#pragma warning (disable:4244)
#ifdef _M_X64
//#ifdef _DEBUG
//#pragma comment (lib, "lpng144/x64/libpng14_64d.lib")
//#else
#pragma comment (lib, "lpng16/libpng16.lib")
#pragma comment (lib, "zlib/zlib1_64.lib")
//#endif
#else
#ifdef _DEBUG
#pragma comment (lib, "lpng144/libpng14d.lib")
#else
#pragma comment (lib, "lpng144/libpng14.lib")
#endif
#endif
#include "lpng16/png.h"

#define _CRT_SECURE_NO_WARNINGS	
#define _CRT_SECURE_NO_WARNINGS_GLOBALS
#pragma warning(disable : 4996)  
static int writePng(const char *fileName, unsigned char **rowsp, int w, int h, bool normalize)
{
	// defaults 
	const int colortype = PNG_COLOR_TYPE_RGBA;
	const int bitdepth = 8;
	png_structp png_ptr = NULL;
	png_infop info_ptr = NULL;
	png_bytep *rows = rowsp;

	FILE *fp = NULL;
	std::string doing = "open for writing";
	if (!(fp = fopen(fileName, "wb"))) goto fail;

	if (!png_ptr) {
		doing = "create png write struct";
		if (!(png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL))) goto fail;
	}
	if (!info_ptr) {
		doing = "create png info struct";
		if (!(info_ptr = png_create_info_struct(png_ptr))) goto fail;
	}

	if (setjmp(png_jmpbuf(png_ptr))) goto fail;
	doing = "init IO";
	png_init_io(png_ptr, fp);
	doing = "write header";
	png_set_IHDR(png_ptr, info_ptr, w, h, bitdepth, colortype, PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
	doing = "write info";
	png_write_info(png_ptr, info_ptr);
	doing = "write image";
	png_write_image(png_ptr, rows);
	doing = "write end";
	png_write_end(png_ptr, NULL);
	doing = "write destroy structs";
	png_destroy_write_struct(&png_ptr, &info_ptr);

	fclose(fp);
	return 0;

fail:
	std::cerr << "writePng: could not " << doing << " !\n";
	if (fp) fclose(fp);
	if (png_ptr || info_ptr) png_destroy_write_struct(&png_ptr, &info_ptr);
	return -1;
}


static void dumpNumberedPNG(int counter, std::string prefix, float* field, int xRes, int yRes)
{
	char buffer[256];
	sprintf(buffer, "%04i", counter);
	std::string number = std::string(buffer);

	unsigned char *pngbuf = new unsigned char[xRes*yRes * 4];
	unsigned char **rows = new unsigned char*[yRes];
	float *pfield = field;
	for (int j = 0; j<yRes; j++) {
		for (int i = 0; i<xRes; i++) {
			float val = *pfield;
			if (val > 0) {
				if (val > 1.) val = 1.;
				pngbuf[(j*xRes + i) * 4 + 0] = (unsigned char)(abs(val)*255.);
				pngbuf[(j*xRes + i) * 4 + 1] = (unsigned char)(0);
				pngbuf[(j*xRes + i) * 4 + 2] = (unsigned char)(0);
			}
			else if (val < 0) {
				if (val < -1) val = -1;
				pngbuf[(j*xRes + i) * 4 + 0] = (unsigned char)(0);
				pngbuf[(j*xRes + i) * 4 + 1] = (unsigned char)(0);
				pngbuf[(j*xRes + i) * 4 + 2] = (unsigned char)(abs(val)*255.);
			}
			else {
				pngbuf[(j*xRes + i) * 4 + 0] = (unsigned char)(0);
				pngbuf[(j*xRes + i) * 4 + 1] = (unsigned char)(255);
				pngbuf[(j*xRes + i) * 4 + 2] = (unsigned char)(0);
			}
			pfield++;
			pngbuf[(j*xRes + i) * 4 + 3] = 255;
		}
		rows[j] = &pngbuf[(yRes - j - 1)*xRes * 4];
	}
	std::string filenamePNG = prefix + number + std::string(".png");
	writePng(filenamePNG.c_str(), rows, xRes, yRes, false);
	delete[]rows;
	delete[]pngbuf;
	printf("Writing %s\n", filenamePNG.c_str());
}
