#pragma once
#include "math_utils.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

static bool WriteInFile(const Array2D<double>& data, double deltax, double deltay, std::string filename) {
	std::ofstream file;
	file.open(filename + ".txt");
	file << std::fixed;
	file << data.XSize() << " " << data.YSize() << std::endl;
	for (int i = 0; i < data.XSize(); i++) {
		for (int j = 0; j < data.YSize(); j++) {
			file << std::setprecision(8) << i * deltax << " " << j * deltay << " " << data.Get(i, j) << std::endl;
		}
	}
	file.close();

    return true;
}