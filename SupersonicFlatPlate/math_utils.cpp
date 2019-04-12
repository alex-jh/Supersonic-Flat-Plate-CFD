#include "math_utils.h"

template Array2D<double>;
template Array2D<int>;

template <class T>
Array2D<T>::Array2D(int imax, int jmax) {
	data_.resize(imax);
	for (int i = 0; i < imax; i++) {
		data_[i].resize(jmax);
	}
}


