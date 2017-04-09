#pragma once
#include "meta.h"
#include <fstream>

// Fields are stored with x-axis changing fastest
class Field
{
	const int slabSize;
public:
	const int sizeX;
	const int sizeY;
	const int sizeZ;
	const int totalSize;
	real* content;
	inline int getIndex(int i, int j, int k) const {
		int modular[3] = { i % sizeX, j % sizeY, k % sizeZ };
		i = (modular[0] >= 0) ? modular[0] : modular[0] + sizeX;
		j = (modular[1] >= 0) ? modular[1] : modular[1] + sizeY;
		k = (modular[2] >= 0) ? modular[2] : modular[2] + sizeZ;
		return i + j * sizeX + k * slabSize;
	}
	void writeSlabPreviewToFile(const std::string& filename, int z = -1);
	void copyFrom(const Field*);
	void dumpFieldToFile(std::ofstream& fout);
	void loadFiledFromFile(std::ifstream& fin);
	Field(int sx, int sy, int sz, bool clear = false);
	Field(const Field&);
	Field(const Field*, bool copy = true);
	~Field();
};

