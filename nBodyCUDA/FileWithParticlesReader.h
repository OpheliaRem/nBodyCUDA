#pragma once
#include "FileReader.h"
#include "Particle.h"
#include <iostream>
#include <fstream>

class FileWithParticlesReader : public FileReader
{
	const std::string path;
	Particle*& particles;
	int& n;
public:
	FileWithParticlesReader(
		const std::string path,
		Particle*& particles,
		int& n
	) : path{path}, particles{particles}, n{n} {}

	void readFile() override
	{
		std::ifstream fileParticles;
		fileParticles.open(path);

		char tempString[256];
		fileParticles.getline(tempString, 256, ':');

		fileParticles >> n;
		particles = new Particle[n];

		fileParticles.get();
		fileParticles.get();

		fileParticles.getline(tempString, 256);

		for (int i = 0; i < n; ++i)
		{
			fileParticles >> particles[i].mass;
			fileParticles.get();
			fileParticles >> particles[i].velocity.x >> particles[i].velocity.y >> particles[i].velocity.z;
			fileParticles.get();
			fileParticles >> particles[i].position.x >> particles[i].position.y >> particles[i].position.z;
		}

		fileParticles.close();
	}
};