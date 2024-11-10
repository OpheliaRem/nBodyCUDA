#pragma once
#include "Integrator.h"

class EulerIntegrator : public Integrator
{
public:
	void integrate(
		Particle* particles,
		const int n,
		const DTOParameters& parameters,
		const AccelerationCalculator& calculator
	) const override
	{
		Vector* acceleration = calculator.calculateAcceleration(particles, n);

		for (int i = 0; i < n; ++i)
		{
			particles[i].position = particles[i].position + particles[i].velocity * parameters.timeStep;
			particles[i].velocity = particles[i].velocity + acceleration[i] * parameters.timeStep;
		}
	}
};