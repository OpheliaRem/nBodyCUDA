#pragma once
#include "PhysicalCharacteristic.h"
#include "SimpleMath.h"

class Energy : public PhysicalCharacteristic
{
	std::ofstream file;

public:
	Energy(std::string fileName = "Energy.txt") :
		file{ std::ofstream() }
	{
		file.open(fileName);
	}

	~Energy()
	{
		file.close();
	}

	void checkLawOfConservation(Particle* particles, const int n) override
	{
		double cineticEnergy = 0.0;

		for (int i = 0; i < n; ++i)
		{
			cineticEnergy += particles[i].mass * SimpleMath::square(particles[i].velocity.abs()) / 2.0;
		}

		double potentialEnergy = 0.0;

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				if (i != j)
				{
					Vector distance = particles[j].position - particles[i].position;
					potentialEnergy += particles[i].mass * particles[j].mass / distance.abs();
				}
			}
		}

		file << cineticEnergy + potentialEnergy << "\n";
	}
};