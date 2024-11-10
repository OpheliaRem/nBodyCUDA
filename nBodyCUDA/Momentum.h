#pragma once
#include "PhysicalCharacteristic.h"

class Momentum : public PhysicalCharacteristic
{
	std::ofstream file;

public:
	Momentum(const std::string fileName) : file{ std::ofstream() }
	{
		file.open(fileName);
	}

	Momentum() : Momentum("Momentum.txt") {}

	~Momentum()
	{
		file.close();
	}

	void checkLawOfConservation(Particle* particles, const int n) override
	{
		Vector momentum;
		momentum.setZeroVector();

		for (int i = 0; i < n; ++i)
		{
			momentum = momentum + particles[i].velocity * particles[i].mass;
		}

		file << momentum << "\n";
	}
};