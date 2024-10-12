#include "Particle.h"
#include "Vector.h"
#include <fstream>
#include <iostream>
#include <string>
#include <filesystem>
#include <chrono>
#include <stdio.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

/*Solution of N-body problem with Particle-Particle (direct sum) method.
This program calculates in 3D space,
with nondimensialization (G = 1),
equations of motion are solved with Euler's method.*/


Particle* InitializeNBodySystem(const std::string path, int& n);
void SetInitialParameters(
	const std::string path,
	double& timeStep,
	double& cuttingRadius,
	double& limitOfLoop,
	int& outputFrequency);

__global__ void calculateForce(Particle* particles, const size_t n)
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

		particles[i].acceleration.x = force.x / particles[i].mass;
		particles[i].acceleration.y = force.y / particles[i].mass;
		particles[i].acceleration.z = force.z / particles[i].mass;
	}
}


int main()
{
	cudaError_t cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		return 1;
	}

	int n;
	double timeStep;
	double cuttingRadius;
	double limitOfLoop;
	int outputFrequency;

	SetInitialParameters(
		"initialParameters.txt",
		timeStep,
		cuttingRadius,
		limitOfLoop,
		outputFrequency);

	std::fstream logFile;
	logFile.open("log.txt");

	std::cout << "Reading input file\n";

	auto startReading = std::chrono::high_resolution_clock::now();
	Particle* particles = InitializeNBodySystem("Particles.txt", n);
	auto finishReading = std::chrono::high_resolution_clock::now();
	auto msReading = std::chrono::duration_cast<std::chrono::milliseconds>(finishReading - startReading);

	std::cout << "File is read in " << msReading.count() << " milliseconds\n\n\n";
	logFile << "File is read in " << msReading.count() << " milliseconds\n\n\n";

	Particle* particlesDevice;
	const size_t sizeBytes = n * sizeof(Particle);

	dim3 dimBlock(1024);
	dim3 dimGrid(n / 1024 + 1);


	std::filesystem::path path = L"coordinates";
	if (std::filesystem::exists(path))
	{
		std::filesystem::remove_all(path);
	}

	if (!std::filesystem::create_directory(path))
	{
		std::cout << "Error making a directory\n";
		return 1;
	}

	double time = 0.0;
	long count = 0;
	while (time < limitOfLoop)
	{
		std::cout << "\n\n============================================================\n";
		std::cout << count << " iteration\n";
		std::cout << "============================================================\n";
		logFile << "\n\n============================================================\n";
		logFile << count << " iteration\n";
		logFile << "============================================================\n";


		cudaStatus = cudaMalloc((void**)&particlesDevice, sizeBytes);
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaMalloc failed!");
			return 1;
		}

		cudaStatus = cudaMemcpy(particlesDevice, particles, sizeBytes, cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaMemcpy failed!");
			return 1;
		}


		if (count % outputFrequency == 0)
		{
			std::ofstream fileCoordinates;
			std::string countStr = std::to_string(count);
			fileCoordinates.open("coordinates\\" + countStr + ".csv");
			fileCoordinates << "x;y;z\n";

			std::cout << "\nWriting to file all positions\n";
			auto startWriting = std::chrono::high_resolution_clock::now();
			for (int i = 0; i < n; ++i)
			{
				fileCoordinates << particles[i].position << std::endl;
			}
			auto finishWriting = std::chrono::high_resolution_clock::now();
			auto msWriting = std::chrono::duration_cast<std::chrono::milliseconds>(finishWriting - startWriting);
			std::cout << "Writing to file has finished in " << msWriting.count() << " milliseconds\n\n";
			logFile << "Writing to file has finished in " << msWriting.count() << " milliseconds\n\n";

			fileCoordinates.close();
		}
		


		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		std::cout << "\nCalculating forces on GPU\n";
		cudaEventRecord(start, 0);

		calculateForce<<<dimGrid, dimBlock>>> (particlesDevice, n);
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return -1;
		}

		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		float elapsedTime;
		cudaEventElapsedTime(&elapsedTime, start, stop);
		std::cout << "Calculating forces on GPU finished in " << elapsedTime << " milliseconds\n\n";
		logFile << "Calculating forces on GPU finished in " << elapsedTime << " milliseconds\n\n";

		cudaEventDestroy(start);
		cudaEventDestroy(stop);



		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
			return 1;
		}

		cudaStatus = cudaMemcpy(particles, particlesDevice, sizeBytes, cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess)
		{
			fprintf(stderr, "cudaMemcpy failed!");
			return 1;
		}

		std::cout << "\nCalculating velocities and positions\n";
		auto startEuler = std::chrono::high_resolution_clock::now();
		for (int i = 0; i < n; ++i)
		{
			particles[i].velocity = particles[i].velocity + particles[i].acceleration * timeStep;

			particles[i].position = particles[i].position + particles[i].velocity * timeStep;
		}
		auto finishEuler = std::chrono::high_resolution_clock::now();
		auto msEuler = std::chrono::duration_cast<std::chrono::milliseconds>(finishEuler - startEuler);
		std::cout << "Calculating velocities and positions has finished in " << msEuler.count() << " milliseconds\n\n\n";
		logFile << "Calculating velocities and positions has finished in " << msEuler.count() << " milliseconds\n\n\n";

		time += timeStep;
		++count;
	}

	delete[] particles;
	logFile.close();
	std::system("pause");
	return 0;
}

Particle* InitializeNBodySystem(const std::string path, int& n)
{
	std::ifstream fileParticles;
	fileParticles.open(path);

	char tempString[256];
	fileParticles.getline(tempString, 256, ':');

	fileParticles >> n;
	Particle* particles = new Particle[n];

	fileParticles.get();
	fileParticles.get();

	fileParticles.getline(tempString, 256);

	for (int i = 0; i < n; ++i)
	{
		fileParticles >> particles[i].mass;
		fileParticles.get();
		fileParticles >> particles[i].velocity.x >> particles[i].velocity.y >> particles[i].velocity.z;
		fileParticles.get();
		fileParticles >> particles[i].position.x >> particles[i].position.y >> particles[i].position.z;
	}

	fileParticles.close();
	return particles;
}


void SetInitialParameters(
	const std::string path,
	double& timeStep,
	double& cuttingRadius,
	double& limitOfLoop,
	int& outputFrequency)
{
	std::ifstream initialFile;
	initialFile.open(path);

	char tempString[32];
	initialFile.getline(tempString, 32, ';');

	initialFile.getline(tempString, 32, '=');
	initialFile >> timeStep;

	initialFile.getline(tempString, 32, '=');
	initialFile >> cuttingRadius;

	initialFile.getline(tempString, 32, '=');
	initialFile >> limitOfLoop;

	initialFile.getline(tempString, 32, '=');
	initialFile >> outputFrequency;

	initialFile.close();
}
