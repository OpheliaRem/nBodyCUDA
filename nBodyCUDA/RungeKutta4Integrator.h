#pragma once
#include "Integrator.h"

class RungeKutta4Integrator : public Integrator
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

		Vector* k1 = new Vector[n];
		Vector* k2 = new Vector[n];
		Vector* k3 = new Vector[n];
		Vector* k4 = new Vector[n];
		Vector* l1 = new Vector[n];
		Vector* l2 = new Vector[n];
		Vector* l3 = new Vector[n];
		Vector* l4 = new Vector[n];

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
			temporaryParticles[i].position = particles[i].position + l2[i] / 2.0;
			temporaryParticles[i].velocity = particles[i].velocity + k2[i] / 2.0;
		}

		delete[] acceleration;
		acceleration = calculator.calculateAcceleration(temporaryParticles, n);

		for (int i = 0; i < n; ++i)
		{
			l3[i] = temporaryParticles[i].velocity * parameters.timeStep;
			k3[i] = acceleration[i] * parameters.timeStep;
			temporaryParticles[i].position = particles[i].position + l3[i];
			temporaryParticles[i].velocity = particles[i].velocity + k3[i];
		}

		delete[] acceleration;
		acceleration = calculator.calculateAcceleration(temporaryParticles, n);

		for (int i = 0; i < n; ++i)
		{
			l4[i] = temporaryParticles[i].velocity * parameters.timeStep;
			k4[i] = acceleration[i] * parameters.timeStep;
		}

		for (int i = 0; i < n; ++i)
		{
			particles[i].position = particles[i].position + (l1[i] + l2[i] * 2 + l3[i] * 2 + l4[i]) / 6.0;
			particles[i].velocity = particles[i].velocity + (k1[i] + k2[i] * 2 + k3[i] * 2 + k4[i]) / 6.0;
		}

		delete[] k1;
		delete[] k2;
		delete[] k3;
		delete[] k4;
		delete[] l1;
		delete[] l2;
		delete[] l3;
		delete[] l4;

		delete[] temporaryParticles;

		delete[] acceleration;
	}
};