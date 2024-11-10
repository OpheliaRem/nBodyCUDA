#pragma once
#include "Vector.h"
#include "Particle.h"

class MomentumData
{
	Vector momentum;
	std::ofstream file;

	void calculate(Particle* particles, const int n)
	{
		Vector sum;
		sum.setZeroVector();

		for (int i = 0; i < n; ++i)
		{
			sum = sum + particles[i].velocity * particles[i].mass;
		}

		momentum = sum;
	}
public:
	MomentumData(const std::string fileName) : momentum{Vector()}, file{std::ofstream()}
	{
		file.open(fileName);
	}

	MomentumData() : MomentumData("Momentum.txt") {}

	~MomentumData()
	{
		file.close();
	}

	void checkLawOfConservationOfMomentum(Particle* particles, const int n)
	{
		calculate(particles, n);

		file << momentum << "\n";
	}
};
