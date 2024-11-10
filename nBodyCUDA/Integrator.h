#pragma once
#include "Particle.h"
#include "DTOParameters.h"
#include "AccelerationCalculator.h"

class Integrator
{
public:
	virtual void integrate(
		Particle* particles,
		const int n,
		const DTOParameters& parameters,
		const AccelerationCalculator& calculator
	) const = 0;
};