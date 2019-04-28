#include "Mesh.h"
#include "Domain.h"
#include <fstream>
#include <math.h>
#include <algorithm>

Mesh::Mesh(std::string file, Domain& domain) {
	readData(file);
	createStructures();
	initializeCoordinates(domain);
	setWalls(domain);
	writeDebugFile();
}

Mesh::Mesh() {}

void Mesh::createSimpleRectangularMesh(int x, int y) {
	nodes_.resize(y + 1);
	for (int i = 0; i < y + 1; i++) {
		nodes_[i].resize(x + 1);
		for (int j = 0; j < x + 1; j++) {
			nodes_[i][j].idxX = j;
			nodes_[i][j].idxY = i;
			nodes_[i][j].x = 1.0*j / x;
			nodes_[i][j].y = 1.0*i / y;
		}
	}
}

Mesh::~Mesh() {
}

void Mesh::readData(std::string& file) {
	std::ifstream ifs;
	ifs.open(file, std::ifstream::in);

	std::string s;

	std::getline(ifs, s);

	ifs >> yLength_ >> xLength_;

	std::getline(ifs, s);
	std::getline(ifs, s);

	int n;
	int tmp;

	ifs >> n;

	horizontalSizes_.clear();
	horizontalSizes_.push_back(0);
	for (int i = 0; i < n; i++) {
		ifs >> tmp;
		horizontalSizes_.push_back(tmp);
	}

	std::getline(ifs, s);
	std::getline(ifs, s);

	ifs >> n;

	verticalSizes_.clear();
	verticalSizes_.push_back(0);
	for (int i = 0; i < n; i++) {
		ifs >> tmp;
		verticalSizes_.push_back(tmp);
	}
}

void Mesh::createStructures() {
	height_ = 1;
	for (auto& i : verticalSizes_) {
		height_ += i;
	}
	width_ = 1;
	for (auto& i : horizontalSizes_) {
		width_ += i;
	}

	nodes_.resize(height_, std::vector<Node>(width_));
}

void Mesh::initializeCoordinates(Domain& domain) {
	initializeCoordinatesX(domain);
	initializeCoordinatesY(domain);

	for (int i = 0; i < height_; i++) {
		for (int j = 0; j < width_; j++) {
			nodes_[i][j].idxX = i;
			nodes_[i][j].idxY = j;
		}
	}
}

void Mesh::initializeCoordinatesX(Domain& domain) {
	for (int i = 1; i < horizontalSizes_.size(); i++) {
		horizontalSizes_[i] += horizontalSizes_[i - 1];
	}

	double step = 1.0 / domain.getWidth() * xLength_;

	for (int k = 0; k < height_; k++) {
		int idx = 0;
		int st = -1;
		int end = -1;
		double x = 0.0;
		int pxl = 0;
		for (int i = 0; i < width_; i++) {
			if (i == horizontalSizes_[idx] && i != width_-1) {
				idx++;
				st = -1;
				end = domain.getWidth();
				for (int j = 0; j < domain.getWidth(); j++) {
					if (domain.getHorizontalZone(j) == idx) {
						if (st == -1) {
							st = j;
						}
					}
					else {
						if (st != -1 && end == domain.getWidth()) {
							end = j;
						}
					}
				}
			}

			double interval = 1.0*(end - st) / domain.getWidth();

			nodes_[k][i].x = (1.0*st / domain.getWidth() + interval * (i - horizontalSizes_[idx - 1]) /
				(horizontalSizes_[idx] - horizontalSizes_[idx - 1])) * xLength_;

			while (nodes_[k][i].x >= x + step) {
				x += step;
				pxl = std::min(pxl+1, domain.getWidth() - 1);
			}

			nodes_[k][i].domainX = pxl;
		}
	}
}

void Mesh::initializeCoordinatesY(Domain& domain) {
	for (int i = 1; i < verticalSizes_.size(); i++) {
		verticalSizes_[i] += verticalSizes_[i - 1];
	}

	double step = 1.0 / domain.getHeight() * yLength_;

	for (int k = 0; k < width_; k++) {
		int idx = 0;
		int st = -1;
		int end = -1;
		double y = 0.0;
		int pxl = 0;
		for (int i = 0; i < height_; i++) {
			if (i == verticalSizes_[idx] && i != height_ - 1) {
				idx++;
				st = -1;
				end = domain.getHeight();
				for (int j = 0; j < domain.getHeight(); j++) {
					if (domain.getVerticalZone(j) == idx) {
						if (st == -1) {
							st = j;
						}
					}
					else {
						if (st != -1 && end == domain.getHeight()) {
							end = j;
						}
					}
				}
			}

			double interval = 1.0*(end - st) / domain.getHeight();

			nodes_[i][k].y = (1.0*st / domain.getHeight() + interval * (i - verticalSizes_[idx - 1]) /
				(verticalSizes_[idx] - verticalSizes_[idx - 1])) * yLength_;

			while (nodes_[i][k].y >= y + step) {
				y += step;
				pxl = std::min(pxl+1, domain.getHeight() - 1);
			}
			nodes_[i][k].domainY = pxl;
		}
	}
}

void Mesh::writeDebugFile() {
	std::ofstream ofs;
	ofs.open("mesh_debug.txt", std::ofstream::out);

	ofs << nodes_[0].size() << std::endl;
	for (int i = 0; i < nodes_[0].size(); i++) {
		ofs << nodes_[0][i].x << " ";
	}
	ofs << std::endl;

	ofs << nodes_.size() << std::endl;
	for (int i = 0; i < nodes_.size(); i++) {
		ofs << nodes_[i][0].y << " ";
	}
	ofs << std::endl;
}

std::vector<std::vector<Node*> > Mesh::findLeftBoundary(Domain& domain, int size) {

	std::vector<Node*> boundary_points;

	for (int i = 0; i < height_; i++) {
		for (int j = 1; j < width_; j++) {
			if (domain.getDomain(nodes_[i][j].domainY, nodes_[i][j].domainX) != 0
				&& (domain.getDomain(nodes_[i][j].domainY, nodes_[i][std::max(0, j - size)].domainX) == 0
					|| (domain.getDomain(nodes_[std::min(height_ - 1, i + size)][j].domainY, nodes_[i][std::max(0, j - size)].domainX) == 0
						&& domain.getDomain(nodes_[std::min(height_ - 1, i + size)][j].domainY, nodes_[i][j].domainX) != 0)
					|| (domain.getDomain(nodes_[std::max(0, i - size)][j].domainY, nodes_[i][std::max(0, j - size)].domainX) == 0
						&& domain.getDomain(nodes_[std::max(0, i - size)][j].domainY, nodes_[i][j].domainX) != 0))) {
				boundary_points.push_back(&nodes_[i][j]);
			}
		}
	}

	std::vector<std::vector<Node*> > intervals_boundary;
	std::vector<std::vector<Node*> > tmp;
	tmp.push_back(boundary_points);

	for (int k = 0; k < 2; k++) {

		intervals_boundary.clear();
		for (auto& v : tmp) {

			if (k == 0) {
				std::sort(v.begin(), v.end(),
					[](const Node* a, const Node* b) -> bool
				{
					return a->idxY < b->idxY;
				});
			}
			else {
				std::sort(v.begin(), v.end(),
					[](const Node* a, const Node* b) -> bool
				{
					return a->idxX < b->idxX;
				});
			}

			for (int i = 0; i < v.size(); i++) {
				if (k == 0) {
					if (i == 0 || abs(v[i]->idxY -
						v[i-1]->idxY) > size) {
						intervals_boundary.push_back(std::vector<Node*>());
					}
				}
				else {
					if (i == 0 || abs(v[i]->idxX -
						v[i - 1]->idxX) > size) {
						intervals_boundary.push_back(std::vector<Node*>());
					}
				}
				intervals_boundary.back().push_back(v[i]);
			}
		}
		tmp = intervals_boundary;
	}

	return intervals_boundary;
}

std::vector<std::vector<Node*> > Mesh::findRightBoundary(Domain& domain, int size) {

	std::vector<Node*> boundary_points;

	for (int i = 0; i < height_; i++) {
		for (int j = 0; j < width_ - 1; j++) {
			if (domain.getDomain(nodes_[i][j].domainY, nodes_[i][j].domainX) != 0
				&& domain.getDomain(nodes_[i][j].domainY, nodes_[i][std::min(width_ - 1, j + size)].domainX) == 0) {
				boundary_points.push_back(&nodes_[i][j]);
			}
		}
	}

	std::sort(boundary_points.begin(), boundary_points.end(),
		[](const Node* a, const Node* b) -> bool
	{
		return a->y < b->y;
	});

	std::vector<std::vector<Node*> > intervals_boundary;

	for (int i = 0; i < boundary_points.size(); i++) {
		if (i == 0 || abs(boundary_points[i]->idxY -
			intervals_boundary.back()[0]->idxY) > size) {
			intervals_boundary.push_back(std::vector<Node*>());
		}
		intervals_boundary.back().push_back(boundary_points[i]);
	}

	return intervals_boundary;
}

std::vector<std::vector<Node*> > Mesh::findUpperBoundary(Domain& domain, int size) {

	std::vector<Node*> boundary_points;

	for (int i = 0; i < height_ - 1; i++) {
		for (int j = 0; j < width_; j++) {
			if (domain.getDomain(nodes_[i][j].domainY, nodes_[i][j].domainX) != 0
					&& (domain.getDomain(nodes_[std::min(height_ - 1, i + size)][j].domainY, nodes_[i][j].domainX) == 0
					|| (domain.getDomain(nodes_[std::min(height_ - 1, i + size)][j].domainY, nodes_[i][std::min(width_ - 1, j + size)].domainX) == 0
						&& domain.getDomain(nodes_[i][j].domainY, nodes_[i][std::min(width_ - 1, j + size)].domainX) != 0)
					|| (domain.getDomain(nodes_[std::min(height_ - 1, i + size)][j].domainY, nodes_[i][std::max(0, j - size)].domainX) == 0
						&& domain.getDomain(nodes_[i][j].domainY, nodes_[i][std::max(0, j - size)].domainX) != 0))) {
				boundary_points.push_back(&nodes_[i][j]);
			}
		}
	}

	std::sort(boundary_points.begin(), boundary_points.end(),
		[](const Node* a, const Node* b) -> bool
	{
		return a->y < b->y;
	});

	std::vector<std::vector<Node*> > intervals_boundary;

	for (int i = 0; i < boundary_points.size(); i++) {
		if (i == 0 || abs(boundary_points[i]->idxX -
			intervals_boundary.back()[0]->idxX) > size) {
			intervals_boundary.push_back(std::vector<Node*>());
		}
		intervals_boundary.back().push_back(boundary_points[i]);
	}

	return intervals_boundary;
}

std::vector<std::vector<Node*> > Mesh::findDownBoundary(Domain& domain, int size) {

	std::vector<Node*> boundary_points;

	for (int i = 1; i < height_; i++) {
		for (int j = 0; j < width_; j++) {
			if (domain.getDomain(nodes_[i][j].domainY, nodes_[i][j].domainX) != 0
					&& (domain.getDomain(nodes_[std::max(0, i - size)][j].domainY, nodes_[i][j].domainX) == 0
					|| (domain.getDomain(nodes_[std::max(0, i - size)][j].domainY, nodes_[i][std::max(0, j - size)].domainX) == 0
						&& domain.getDomain(nodes_[i][j].domainY, nodes_[i][std::max(0, j - size)].domainX) != 0)
					|| (domain.getDomain(nodes_[std::max(0, i - size)][j].domainY, nodes_[i][std::min(width_-1, j + size)].domainX == 0
						&& domain.getDomain(nodes_[i][j].domainY, nodes_[i][std::min(width_ - 1, j + size)].domainX) != 0)))) {
				boundary_points.push_back(&nodes_[i][j]);
			}
		}
	}

	std::sort(boundary_points.begin(), boundary_points.end(),
		[](const Node* a, const Node* b) -> bool
	{
		return a->y < b->y;
	});

	std::vector<std::vector<Node*> > intervals_boundary;

	for (int i = 0; i < boundary_points.size(); i++) {
		if (i == 0 || abs(boundary_points[i]->idxX -
			intervals_boundary.back()[0]->idxX) > size) {
			intervals_boundary.push_back(std::vector<Node*>());
		}
		intervals_boundary.back().push_back(boundary_points[i]);
	}

	return intervals_boundary;
}

void Mesh::setWalls(Domain& domain) {
	for (int i = 0; i < height_; i++) {
		for (int j = 0; j < width_; j++) {
			nodes_[i][j].isWall = (domain.getDomain(nodes_[i][j].domainY, nodes_[i][j].domainX) == 0);
		}
	}
}