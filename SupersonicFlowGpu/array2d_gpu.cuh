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

	__device__ T& Get(int x, int y) { return *(data_ + x * sizey_ + y); }

	__device__ const T& Get(int x, int y) const { return *(data_ + x * sizey_ + y); }
	__device__ void Set(int x, int y, const T& val) { *(data_ + x * sizey_ + y) = val; }

	__device__ const T& Get(int x) const { return *(data_ + x); }
	__device__ void Set(int x, const T& val) { *(data_ + x) = val; }

	void copyArray2D(int n, const Array2DGpu<T>* other, int numThreads, int numBlocks);

	void memsetToArray(int start, int end, double val, int numThreads, int numBlocks);

	void memsetToArrayLinear(int start, int end, double val_s, double val_e, int numThreads, int numBlocks);

	void partialSum(int start, int end, int numThreads, int numBlocks);

public:
	__host__ __device__ T* ptr() const { return data_; }

	__host__ __device__ void* data_adress() const { return (void*)&data_; }
private:
	T* data_;
	int sizex_;
	int sizey_;
};

template <class T>
class ArrayWrapper4D
{
public:

	__device__ __host__ Array2DGpu<T>* operator[](int i) { return arr[i]; }

	Array2DGpu<T>* arr[4];
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
__global__ void copyArray2D(const Array2DGpu<T>* buf, T* a, int n)
{
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;

	while (idx < n) {

		a[idx] = buf->Get(idx);

		idx += gridSize;
	}
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
void copyArray2DFromArrayGpu(int n, const Array2DGpu<T>* a, Array2DGpu<T>* b,
	int numThreads, int numBlocks);

template <class T>
__global__ void copyArray2DGpu(int n, const Array2DGpu<T>* a, Array2DGpu<T>* b);

template <class T>
void Array2DGpu<T>::copyArray2D(int n, const Array2DGpu<T>* other,
	int numThreads, int numBlocks) {
	copyArray2DFromArrayGpu<T>(n, other, this, numThreads, numBlocks);
}

template <class T>
__global__ void memsetToArrayGpu(int start, int end, double val, Array2DGpu<T>* data);

template <class T>
void memsetToArrayGpu(int start, int end, double val, int numThreads, int numBlocks, Array2DGpu<T>* data);

template <class T>
void Array2DGpu<T>::memsetToArray(int start, int end, double val, int numThreads, int numBlocks) {
	memsetToArrayGpu<T>(start, end, val, numThreads, numBlocks, this);
}

template <class T>
__global__ void memsetToArrayGpuLinear(int start, int end, double val_s, double val_e, Array2DGpu<T>* data);

template <class T>
void memsetToArrayGpuLinear(int start, int end, double val_s, double val_e, int numThreads, int numBlocks, Array2DGpu<T>* data);

template <class T>
void Array2DGpu<T>::memsetToArrayLinear(int start, int end, double val_s, double val_e, int numThreads, int numBlocks) {
	memsetToArrayGpuLinear<T>(start, end, val_s, val_e, numThreads, numBlocks, this);
}

template <class T>
__global__ void partialSumGPU(int start, int end, Array2DGpu<T>* data);

template <class T>
void partialSumGPU(int start, int end, int numThreads, int numBlocks, Array2DGpu<T>* data);

template <class T>
void Array2DGpu<T>::partialSum(int start, int end, int numThreads, int numBlocks) {
	partialSumGPU<T>(start, end, numThreads, numBlocks, this);
}