#pragma once

#include <string>

class Domain
{
public:
	Domain(std::string file);
	~Domain();

	unsigned char* readBMP(std::string& filename);

	void processData(unsigned char* data);
	void findBoundaries(unsigned char* data, int width, int height, int& startWidth,
		int& endWidth, int& startHeight, int& endHeight);
	void calcDomain(unsigned char* data, int startWidth, int startHeight);
	void calcHorizontalZones(unsigned char* data, int startWidth, int startHeight);
	void calcVerticalZones(unsigned char* data, int startWidth, int startHeight);

	int getWidth() { return width_; }
	int getHeight() { return height_; }

	int getHorizontalZone(int j) { return horizontalZone_[j]; }
	int getVerticalZone(int j) { return verticalZone_[j]; }

	int getDomain(int i, int j) { return domain_[i*width_ + j]; }

private:
	int* horizontalZone_;
	int* verticalZone_;
	int* domain_;
	int width_;
	int height_;
	int dataWidth_;
	int dataHeight_;
};

