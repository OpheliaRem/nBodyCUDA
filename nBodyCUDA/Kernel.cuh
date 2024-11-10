#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "Particle.h"

__host__ Vector* CUDACalculateAcceleration(Particle* particles, const int n);

__global__ void calculateForce(Particle* particles, const size_t n, Vector* acceleration);