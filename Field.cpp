#include "stdafx.h"
#include "Field.h"
#include <fstream>

real Field::getMax() const
{
	real maxVal = std::numeric_limits<real>::min();
	for (int i = 0; i < totalSize; ++i) {
		// Make sure no Nan exists.
		CHECK_EQ(content[i], content[i]);
		if (content[i] > maxVal) maxVal = content[i];
	}
	return maxVal;
}

void Field::writeSlabPreviewToFile(const std::string& filename, int z /* = -1 */)
{
	CHECK_LT(z, sizeZ);
	if (z < 0) z = sizeZ / 2;
	std::ofstream fout(filename);
	for (int j = 0; j < sizeY; ++j) {
		for (int i = 0; i < sizeX; ++i) {
			fout << content[getIndex(i, j, z)] << '\t';
		}
		fout << '\n';
	}
	fout.close();
	LOG(INFO) << "Slice preview (z = " << z << ") dumped to " << filename;
}

void Field::copyFrom(const Field *f)
{
	CHECK_NOTNULL(f);
	CHECK_EQ(sizeX, f->sizeX);
	CHECK_EQ(sizeY, f->sizeY);
	CHECK_EQ(sizeZ, f->sizeZ);
	memcpy(content, f->content, sizeof(real) * totalSize);
}

void Field::dumpFieldToFile(std::ofstream & fout)
{
	fout.write((char*)this->content, sizeof(real) * totalSize);
}

void Field::loadFiledFromFile(std::ifstream & fin)
{
	fin.read((char*)this->content, sizeof(real) * totalSize);
}

void Field::initTopSlice()
{
	topSlice = new real[sizeX * sizeZ];
}

Field::Field(int sx, int sy, int sz, bool clear /* = false */)
	: sizeX(sx), sizeY(sy), sizeZ(sz), slabSize(sx * sy), totalSize(sx * sy * sz)
{
	content = new real[totalSize];
	if (clear) {
		memset(content, 0, sizeof(real) * totalSize);
	}
	topSlice = 0;
}

Field::Field(const Field &f)
	: sizeX(f.sizeX), sizeY(f.sizeY), sizeZ(f.sizeZ), slabSize(f.sizeX * f.sizeY),
	totalSize(f.sizeX * f.sizeY * f.sizeZ)
{
	content = new real[totalSize];
	memcpy(content, f.content, sizeof(real) * totalSize);
	topSlice = 0;
}

Field::Field(const Field *f, bool copy /* = true */)
	: sizeX(f->sizeX), sizeY(f->sizeY), sizeZ(f->sizeZ), slabSize(f->sizeX * f->sizeY),
	totalSize(f->sizeX * f->sizeY * f->sizeZ)
{
	content = new real[totalSize];
	if (copy) memcpy(content, f->content, sizeof(real) * totalSize);
	topSlice = 0;
}

Field::~Field()
{
	delete[] content;
	if (topSlice) delete[] topSlice;
}
