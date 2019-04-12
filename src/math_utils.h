#pragma once
#include <vector>

template <class T>
class Array2D
{
public:
	Array2D() {}
	Array2D(int imax, int jmax);
	~Array2D() {};

	void operator=(const Array2D<T>& other) {
		if (other.XSize() > XSize()) {
			data_.resize(other.XSize());
		}
		if (other.YSize() > YSize()) {
			for (int i = 0; i < XSize(); i++) {
				data_[i].resize(other.YSize());
			}
		}
		for (int i = 0; i < XSize(); i++) {
			for (int j = 0; j < YSize(); j++) {
				data_[i][j] = other.Get(i, j);
			}
		}
	}

	int XSize() const { return (int)data_.size(); }
	int YSize() const  { return data_.size() > 0 ? (int)data_[0].size() : 0; }

	T& Get(int x, int y) { return data_[x][y]; }
	const T& Get(int x, int y) const { return data_[x][y]; }

	void Set(int x, int y, const T& val) { data_[x][y] = val; }

private:
	std::vector<std::vector<T> > data_;
};

