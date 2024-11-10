#include "Kernel.cuh"

__host__ Vector* CUDACalculateAcceleration(Particle* particles, const int n)
{
	cudaError_t cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		throw 1;
	}

	Particle* particlesDevice;
	const size_t sizeParticlesInBytes = n * sizeof(Particle);

	Vector* accelerationDevice;
	const size_t sizeAccelerationInBytes = n * sizeof(Vector);

	dim3 dimBlock(1024);
	dim3 dimGrid(n / 1024 + 1);

	cudaStatus = cudaMalloc((void**)&particlesDevice, sizeParticlesInBytes);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "cudaMalloc failed!");
		throw 1;
	}

	cudaStatus = cudaMemcpy(particlesDevice, particles, sizeParticlesInBytes, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "cudaMemcpy failed!");
		throw 1;
	}

	cudaStatus = cudaMalloc((void**)&accelerationDevice, sizeAccelerationInBytes);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "cudaMalloc failed!");
		throw 1;
	}

	calculateForce<<<dimGrid, dimBlock>>>(particlesDevice, n, accelerationDevice);
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
		throw - 1;
	}

	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		throw 1;
	}

	Vector* acceleration = new Vector[n];

	cudaStatus = cudaMemcpy(acceleration, accelerationDevice, sizeAccelerationInBytes, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "cudaMemcpy failed!");
		throw 1;
	}

	cudaFree(particlesDevice);
	cudaFree(accelerationDevice);

	return acceleration;
}

__global__ void calculateForce(Particle* particles, const size_t n, Vector* acceleration)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < n)
	{
		Vector force;
		force.x = 0.0;
		force.y = 0.0;
		force.z = 0.0;

		for (int j = 0; j < n; ++j)
		{
			if (i != j)
			{
				double distanceX = particles[j].position.x - particles[i].position.x;
				double distanceY = particles[j].position.y - particles[i].position.y;
				double distanceZ = particles[j].position.z - particles[i].position.z;

				double vector = sqrt(distanceX * distanceX + distanceY * distanceY + distanceZ * distanceZ);
				double denominator = vector * vector * vector;

				force.x = force.x + distanceX * particles[i].mass * particles[j].mass / denominator;
				force.y = force.y + distanceY * particles[i].mass * particles[j].mass / denominator;
				force.z = force.z + distanceZ * particles[i].mass * particles[j].mass / denominator;
			}
		}

		acceleration[i].x = force.x / particles[i].mass;
		acceleration[i].y = force.y / particles[i].mass;
		acceleration[i].z = force.z / particles[i].mass;
	}
}