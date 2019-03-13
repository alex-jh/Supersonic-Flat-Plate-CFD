#pragma once
#include <string>
#include <vector>

class Domain;

class Node {
public:
	double x;
	double y;
	int idxX;
	int idxY;
	int domainX;
	int domainY;
	bool isWall;
};

class Mesh {
public:
	Mesh(std::string file, Domain& domain);
	Mesh();
	~Mesh();

	void readData(std::string& file);
	void createStructures();
	void initializeCoordinates(Domain& domain);
	void initializeCoordinatesX(Domain& domain);
	void initializeCoordinatesY(Domain& domain);

	void createSimpleRectangularMesh(int x, int y);

	std::vector<std::vector<Node*> > findLeftBoundary(Domain& domain, int size=1);
	std::vector<std::vector<Node*> > findRightBoundary(Domain& domain, int size = 1);
	std::vector<std::vector<Node*> > findUpperBoundary(Domain& domain, int size = 1);
	std::vector<std::vector<Node*> > findDownBoundary(Domain& domain, int size = 1);

	void writeDebugFile();

	int getXNodes() { return (int) nodes_[0].size(); }
	int getYNodes() { return (int) nodes_.size(); }

	std::vector<Node>& operator[](int i) { return nodes_[i]; }

	bool isWall(int y, int x) { return nodes_[y][x].isWall; }

	void setWalls(Domain& domain);

private:
	int height_;
	int width_;
	double xLength_;
	double yLength_;
	std::vector<int> horizontalSizes_;
	std::vector<int> verticalSizes_;

	std::vector<std::vector<Node> > nodes_;
};

