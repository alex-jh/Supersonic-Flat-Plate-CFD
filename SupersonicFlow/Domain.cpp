#include "Domain.h"

#include <iostream>

Domain::Domain(std::string file) {
	unsigned char* data = readBMP(file);
	processData(data);

	delete data;
}

Domain::~Domain() {
	delete domain_;
	delete horizontalZone_;
	delete verticalZone_;
}

unsigned char* Domain::readBMP(std::string& filename) {
	int i;
	FILE *f;
	errno_t err;
	err = fopen_s(&f, filename.c_str(), "rb");
	unsigned char info[54];
	fread(info, sizeof(unsigned char), 54, f); // read the 54-byte header

	dataWidth_ = *(int*)&info[18];
	dataHeight_ = *(int*)&info[22];

	int size = 3 * dataWidth_ * dataHeight_;
	unsigned char* data = new unsigned char[size]; // allocate 3 bytes per pixel
	fread(data, sizeof(unsigned char), size, f); // read the rest of the data at once
	fclose(f);

	for (i = 0; i < size; i += 3)
	{
		unsigned char tmp = data[i];
		data[i] = data[i + 2];
		data[i + 2] = tmp;
	}

	return data;
}

void Domain::processData(unsigned char* data) {

	int width = dataWidth_;
	int height = dataHeight_;

	int startWidth = 0;
	int endWidth = width - 1;
	int startHeight = 0;
	int endHeight = height - 1;

	findBoundaries(data, width, height, startWidth, endWidth, startHeight, endHeight);

	width_ = endWidth - startWidth + 1;
	height_ = endHeight - startHeight + 1;

	calcDomain(data, startWidth, startHeight);

	calcHorizontalZones(data, startWidth, startHeight);

	calcVerticalZones(data, startWidth, startHeight);

}

void Domain::findBoundaries(unsigned char* data, int width, int height, int& startWidth,
	int& endWidth, int& startHeight, int& endHeight) {

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (data[((dataHeight_ - 1 - i)*dataWidth_ + j) * 3] == 0) {
				startWidth = j;
				break;
			}
		}
		for (int j = width - 1; j >= 0; j--) {
			if (data[((dataHeight_ - 1 - i)*dataWidth_ + j) * 3] == 0) {
				endWidth = j;
				break;
			}
		}
	}

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			if (data[((dataHeight_ - 1 - j)*dataWidth_ + i) * 3] == 0) {
				startHeight = j;
				break;
			}
		}
		for (int j = height - 1; j >= 0; j--) {
			if (data[((dataHeight_ - 1 - j)*dataWidth_ + i) * 3] == 0) {
				endHeight = j;
				break;
			}
		}
	}

}

void Domain::calcDomain(unsigned char* data, int startWidth, int startHeight) {
	domain_ = new int[width_*height_];

	for (int i = 0; i < height_; i++) {
		for (int j = 0; j < width_; j++) {
			if (data[((dataHeight_ - 1 - i - startHeight)*dataWidth_ + j + startWidth) * 3 + 2] == 0) {
				domain_[(i*width_ + j)] = 0;
			}
			else {
				domain_[(i*width_ + j)] = 1;
			}
		}
	}
}

void Domain::calcHorizontalZones(unsigned char* data, int startWidth, int startHeight) {
	horizontalZone_ = new int[width_];

	int i = height_ / 2;
	int idx = 1;

	for (int j = 0; j < width_; j++) {
		if (data[((dataHeight_ - 1 - i - startHeight)*dataWidth_ + j + startWidth) * 3 + 2] == 255
			&& data[((dataHeight_ - 1 - i - startHeight)*dataWidth_ + j + startWidth) * 3 + 1] == 0) {
			idx++;
		}
		horizontalZone_[j] = idx;
	}
}

void Domain::calcVerticalZones(unsigned char* data, int startWidth, int startHeight) {
	verticalZone_ = new int[height_];

	int j = width_ / 2;
	int idx = 1;
	for (int i = 0; i < height_; i++) {
		if (data[((dataHeight_ - 1 - i - startHeight)*dataWidth_ + j + startWidth) * 3 + 2] == 255
			&& data[((dataHeight_ - 1 - i - startHeight)*dataWidth_ + j + startWidth) * 3 + 1] == 0) {
			idx++;
		}
		verticalZone_[i] = idx;
	}
}