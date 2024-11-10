#pragma once
#include "Particle.h"

class AccelerationCalculator
{
public:
	virtual Vector* calculateAcceleration(Particle* particles, const int n) const = 0;
};