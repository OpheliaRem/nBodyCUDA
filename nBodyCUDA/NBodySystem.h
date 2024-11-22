#pragma once
#include "Particle.h"
#include "Vector.h"
#include "SolvingDEMethod.h"
#include "DTOParameters.h"
#include "SimpleMath.h"
#include "NBodyInitializer.h"
#include "Integrator.h"
#include "EulerIntegrator.h"
#include "PredictorCorrectorIntegrator.h"
#include "RungeKutta2Integrator.h"
#include "RungeKutta4Integrator.h"
#include "AccelerationCalculator.h"
#include "SequentialAccelerationCalculator.h"
#include "AccelerationCalculatorCUDA.h"
#include "FileWriter.h"
#include "FileCoordinatesWriter.h"
#include "PhysicalCharacteristic.h"
#include "Momentum.h"
#include "Energy.h"
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>

class NBodySystem
{
	Particle* particles;
	int n;

	DTOParameters parameters;

	NBodyInitializer* initializer;

	Integrator* integrator;

	AccelerationCalculator* calculator;

	PhysicalCharacteristic* momentum;

	PhysicalCharacteristic* energy;


	void initializeSystemWithLog(std::ofstream& logFile) {
		std::cout << "Initializing system\n";

		auto start = std::chrono::high_resolution_clock::now();

		initializer->initialize();

		auto finish = std::chrono::high_resolution_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);

		std::cout << "System is initialized in " << ms.count() << " milliseconds\n\n\n";
		logFile << "System is initialized in " << ms.count() << " milliseconds\n\n\n";
	}

	void writeToFileWithLog(const long count, std::ofstream& logFile) {
		if (count % parameters.outputFrequency != 0)
			return;

		std::cout << "\nWriting to file all positions\n";

		FileWriter* writer = new FileCoordinatesWriter(
			count, particles, n
		);

		auto start = std::chrono::high_resolution_clock::now();

		writer->write();

		auto finish = std::chrono::high_resolution_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);

		delete writer;

		std::cout << "Writing to file has finished in " << ms.count() << " milliseconds\n\n";
		logFile << "Writing to file has finished in " << ms.count() << " milliseconds\n\n";
	}

	void checkLawsOfConservation(const long count) {
		long limitOfIterations = (long)(parameters.limitOfLoop / parameters.timeStep) + 1L;

		if (
			parameters.physicalCharacteristicsOutputFrequency >= limitOfIterations ||
			count % parameters.physicalCharacteristicsOutputFrequency != 0
			)
			return;

		momentum->checkLawOfConservation(particles, n);
		energy->checkLawOfConservation(particles, n);
	}

	void calculateVelocitiesAndPositionsWithLog(std::ofstream& logFile) {
		std::cout << "\nCalculating velocities and positions\n";

		auto start = std::chrono::high_resolution_clock::now();

		integrator->integrate(particles, n, parameters, *calculator);

		auto finish = std::chrono::high_resolution_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);

		std::cout << "Calculating velocities and positions has finished in " << ms.count() << " milliseconds\n\n\n";
		logFile << "Calculating velocities and positions has finished in " << ms.count() << " milliseconds\n\n\n";
	}

public:
	NBodySystem() : 
		particles{ nullptr },
		n{ 0 },
		parameters{ DTOParameters() },
		integrator{ nullptr },
		calculator{ nullptr },
		initializer{ new NBodyInitializer(
			"Particles.txt",
			"initialParameters.txt",
			particles,
			n,
			parameters,
			integrator,
			calculator
		)},
		momentum { new Momentum() },
		energy { new Energy() }
	{}

	~NBodySystem()
	{
		delete[] particles;
		delete initializer;
		delete integrator;
		delete calculator;
		delete momentum;
		delete energy;
	}

	void simulate()
	{
		std::ofstream logFile;
		logFile.open("log.txt");

		initializeSystemWithLog(logFile);

		double time = 0.0;
		long limitOfIterations = (long)((parameters.limitOfLoop - time) / parameters.timeStep) + 1L;

		for (long count = 0; count < limitOfIterations; ++count)
		{
			std::cout << "\n\n============================================================\n";
			std::cout << count << " iteration\n";
			std::cout << "============================================================\n";
			logFile << "\n\n============================================================\n";
			logFile << count << " iteration\n";
			logFile << "============================================================\n";

			writeToFileWithLog(count, logFile);

			checkLawsOfConservation(count);

			calculateVelocitiesAndPositionsWithLog(logFile);

			time += parameters.timeStep;
		}

		logFile.close();
	}
};