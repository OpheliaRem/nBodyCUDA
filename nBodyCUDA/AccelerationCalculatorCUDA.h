#pragma once
#include "AccelerationCalculator.h"
#include "Kernel.cuh"

class AccelerationCalculatorCUDA : public AccelerationCalculator
{
	Vector* calculateAcceleration(Particle* particles, const int n) const override {
		return CUDACalculateAcceleration(particles, n);
	}
};