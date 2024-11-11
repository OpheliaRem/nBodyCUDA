#pragma once
#include "FileWriter.h"
#include "Particle.h"
#include <string>

class FileCoordinatesWriter : public FileWriter
{
	long count;
	Particle* particles;
	int n;

public:
	FileCoordinatesWriter(
		long count,
		Particle* particles,
		int n
	) : 
		count{count},
		particles{particles},
		n{n}
	{}

	void write() override
	{
		std::ofstream fileCoordinates;
		std::string countStr = std::to_string(count);
		fileCoordinates.open("coordinates\\" + countStr + ".csv");
		fileCoordinates << "x;y;z\n";

		for (int i = 0; i < n; ++i)
		{
			fileCoordinates << particles[i].position << std::endl;
		}

		fileCoordinates.close();
	}
};