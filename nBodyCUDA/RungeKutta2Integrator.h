#pragma once
#include "Integrator.h"

class RungeKutta2Integrator : public Integrator
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

		Vector* k1 = new Vector[n];
		Vector* k2 = new Vector[n];
		Vector* l1 = new Vector[n];
		Vector* l2 = new Vector[n];

		Particle* temporaryParticles = new Particle[n];

		for (int i = 0; i < n; ++i)
		{
			temporaryParticles[i] = particles[i];
		}

		for (int i = 0; i < n; ++i)
		{
			l1[i] = particles[i].velocity * parameters.timeStep;
			k1[i] = acceleration[i] * parameters.timeStep;
			temporaryParticles[i].position = particles[i].position + l1[i] / 2.0;
			temporaryParticles[i].velocity = particles[i].velocity + k1[i] / 2.0;
		}

		delete[] acceleration;
		acceleration = calculator.calculateAcceleration(temporaryParticles, n);

		for (int i = 0; i < n; ++i)
		{
			l2[i] = temporaryParticles[i].velocity * parameters.timeStep;
			k2[i] = acceleration[i] * parameters.timeStep;
		}

		for (int i = 0; i < n; ++i)
		{
			particles[i].position = particles[i].position + l2[i];
			particles[i].velocity = particles[i].velocity + k2[i];
		}

		delete[] acceleration;
		delete[] k1;
		delete[] k2;
		delete[] l1;
		delete[] l2;
		delete[] temporaryParticles;
	}
};