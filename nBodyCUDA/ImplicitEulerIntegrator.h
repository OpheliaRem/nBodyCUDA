#pragma once
#include "Integrator.h"

class ImplicitEulerIntegrator : public Integrator
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

		Particle* temporaryParticles = new Particle[n];

		for (int i = 0; i < n; ++i)
		{
			temporaryParticles[i] = particles[i];
		}

		for (int i = 0; i < n; ++i)
		{
			temporaryParticles[i].position = particles[i].position + temporaryParticles[i].velocity * parameters.timeStep;
			temporaryParticles[i].velocity = particles[i].velocity + acceleration[i] * parameters.timeStep;
		}

		Vector* intermediateAcceleration = calculator.calculateAcceleration(temporaryParticles, n);

		for (int i = 0; i < n; ++i)
		{
			particles[i].position = particles[i].position + (particles[i].velocity + temporaryParticles[i].velocity) * parameters.timeStep / 2;
			particles[i].velocity = particles[i].velocity + (acceleration[i] + intermediateAcceleration[i]) * parameters.timeStep / 2;
		}

		delete[] temporaryParticles;
	}
};