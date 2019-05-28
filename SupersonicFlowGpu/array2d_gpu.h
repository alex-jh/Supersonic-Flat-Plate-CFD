#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>

template <class T>
class Array2DGpu
{
public:
	Array2DGpu() {}
	__device__ Array2DGpu(int imax, int jmax);
	__device__ ~Array2DGpu();

	__device__ void operator=(const Array2DGpu<T>& other) {
		if (other.XSize() != XSize() || other.YSize() != YSize()) {
			printf("failed copying array contents!");
		}
		//memcpy(data_, other.ptr(), XSize()*YSize());
	}

	__device__ int XSize() const { return sizex_; }
	__device__ int YSize() const { return sizey_; }

	__device__ T& Get(int x, int y) { return *(data_ + x*sizey_ + y); }

	__device__ const T& Get(int x, int y) const { return *(data_ + x * sizey_ + y); }
	__device__ void Set(int x, int y, const T& val) { *(data_ + x * sizey_ + y) = val; }

	__device__ const T& Get(int x) const { return *(data_ + x); }
	__device__ void Set(int x, const T& val) { *(data_ + x) = val; }

	__host__ void copyArray2D(int n, Array2DGpu<T>* other, int numThreads, int numBlocks);

	__host__ void memsetToArray(int start, int end, double val, int numThreads, int numBlocks);

public:
	__host__ __device__ T* ptr() const { return data_; }

	__host__ __device__ void* data_adress() const { return (void*)&data_; }
private:
	T* data_;
	int sizex_;
	int sizey_;
};

template <class T>
__global__ void initArray2D(void* buf, int imax, int jmax)
{
	if (threadIdx.x == 0 && blockIdx.x == 0) {
		new (buf) Array2DGpu<T>(imax, jmax);
	}
}

template <class T>
__global__ void destroyArray2D(Array2DGpu<T>* buf)
{
	if (threadIdx.x == 0 && blockIdx.x == 0) {
		buf->~Array2DGpu<T>();
	}
}

template <class T>
__global__ void copyArray2D(const Array2DGpu<T>* buf, T* a)
{
	a[blockIdx.x*blockDim.x + threadIdx.x] = buf->Get(blockIdx.x, threadIdx.x);
}

template <class T>
__device__ Array2DGpu<T>::Array2DGpu(int imax, int jmax) {
	data_ = new T[imax*jmax];
	memset(data_, T(), imax*jmax);
	sizex_ = imax;
	sizey_ = jmax;
}

template <class T>
__device__ Array2DGpu<T>::~Array2DGpu() {
	delete data_;
}

template <class T>
__global__ void copyArray2DGpu(int n, const Array2DGpu<T>* a, Array2DGpu<T>* b)
{
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;

	while (idx < n) {

		b->Set(idx, a->Get(idx));

		idx += gridSize;
	}
}

template <class T>
__host__ void Array2DGpu<T>::copyArray2D(int n, Array2DGpu<T>* other, int numThreads, int numBlocks) {
	copyArray2DGpu<T> << < numBlocks, numThreads, 1 >> > (n, other, this);
}
