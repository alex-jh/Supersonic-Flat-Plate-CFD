#include "array2d_gpu.cuh"

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

template
void copyArray2DFromArrayGpu(int n, const Array2DGpu<double>* a, Array2DGpu<double>* b,
	int numThreads, int numBlocks);

template <class T>
void copyArray2DFromArrayGpu(int n, const Array2DGpu<T>* a, Array2DGpu<T>* b,
	int numThreads, int numBlocks) {
	copyArray2DGpu<T> << < numBlocks, numThreads, 1 >> > (n, a, b);
}

template <class T>
__global__ void memsetToArrayGpu(int start, int end, double val, Array2DGpu<T>* data) {
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;

	while (idx < end) {

		if (idx >= start) {
			data->Set(idx, val);
		}

		idx += gridSize;
	}
}

template
void memsetToArrayGpu(int start, int end, double val, int numThreads, int numBlocks, Array2DGpu<double>* data);

template <class T>
void memsetToArrayGpu(int start, int end, double val, int numThreads, int numBlocks, Array2DGpu<T>* data) {
	memsetToArrayGpu<T> << < numBlocks, numThreads, 1 >> > (start, end, val, data);
}

template <class T>
__global__ void memsetToArrayGpuLinear(int start, int end, double val_s, double val_e, Array2DGpu<T>* data) {
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;

	while (idx < end) {

		if (idx >= start) {
			double val = end > start + 1 ? val_s + (val_e - val_s) / (end - start - 1) * (idx - start) : val_s;
			data->Set(idx, val);
		}

		idx += gridSize;
	}
}

template
void memsetToArrayGpuLinear(int start, int end, double val_s, double val_e, int numThreads, int numBlocks, Array2DGpu<double>* data);

template <class T>
void memsetToArrayGpuLinear(int start, int end, double val_s, double val_e, int numThreads, int numBlocks, Array2DGpu<T>* data) {
	memsetToArrayGpuLinear<T> << < numBlocks, numThreads, 1 >> > (start, end, val_s, val_e, data);
}

template <class T>
__global__ void partialSumGPU(int start, int end, Array2DGpu<T>* data) {
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx == 0) {
		data->Set(start + end, data->Get(start));
		for (int i = start + 1; i < end; i++) {
			data->Set(i + end, data->Get(i) + data->Get(i + end - 1));
		}
	}
}

template
void partialSumGPU(int start, int end, int numThreads, int numBlocks, Array2DGpu<double>* data);

template <class T>
void partialSumGPU(int start, int end, int numThreads, int numBlocks, Array2DGpu<T>* data) {
	partialSumGPU<T> << < numBlocks, numThreads, 1 >> > (start, end, data);
}
