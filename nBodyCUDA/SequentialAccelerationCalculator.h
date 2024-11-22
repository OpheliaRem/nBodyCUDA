#pragma once
#include "AccelerationCalculator.h"
#include "SimpleMath.h"

class SequentialAccelerationCalculator : public AccelerationCalculator
{
	Vector* calculateAcceleration(Particle* particles, const int n) const override {
		Vector* acceleration = new Vector[n];

		for (int i = 0; i < n; ++i) {
			Vector force;
			force.setZeroVector();

			for (int j = 0; j < n; ++j) {
				if (i == j) {
					continue;
				}

				Vector distance = particles[j].position - particles[i].position;
				Vector toAdd = -distance * particles[i].mass * particles[j].mass / SimpleMath::cube(distance.abs());
				force = force + toAdd;
			}

			acceleration[i] = force / particles[i].mass;
		}

		return acceleration;
	}
};