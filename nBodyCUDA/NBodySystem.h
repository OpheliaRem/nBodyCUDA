#pragma once
#include "Particle.h"
#include "Vector.h"
#include "SolvingDEMethod.h"
#include "DTOParameters.h"
#include "SimpleMath.h"
#include "Integrator.h"
#include "EulerIntegrator.h"
#include "ImplicitEulerIntegrator.h"
#include "RungeKutta4Integrator.h"
#include "AccelerationCalculator.h"
#include "SequentialAccelerationCalculator.h"
#include "AccelerationCalculatorCUDA.h"
#include "PhysicalCharacteristic.h"
#include "Momentum.h"
#include "Energy.h"
#include <fstream>
#include <iostream>
#include <string>
#include <filesystem>
#include <chrono>

class NBodySystem
{
	//fields

	Particle* particles;
	int n;

	DTOParameters parameters;

	Integrator* integrator;

	AccelerationCalculator* calculator;

	PhysicalCharacteristic* momentum;

	PhysicalCharacteristic* energy;

	//methods

	void setInitialParameters(const std::string path)
	{
		std::ifstream initialFile;
		initialFile.open(path);

		char tempString[320];
		initialFile.getline(tempString, 320, ';');

		initialFile.getline(tempString, 320, '=');
		initialFile >> parameters.timeStep;

		initialFile.getline(tempString, 320, '=');
		initialFile >> parameters.limitOfLoop;

		initialFile.getline(tempString, 320, '=');
		initialFile >> parameters.outputFrequency;

		initialFile.getline(tempString, 320, '=');
		initialFile >> parameters.physicalCharacteristicsOutputFrequency;

		initialFile.getline(tempString, 320, '=');
		std::string solvingDEMethod;
		initialFile >> solvingDEMethod;

		if (solvingDEMethod == "Euler;")
		{
			parameters.solvingDEMethod = SolvingDEMethod::Euler;
		}

		else if (solvingDEMethod == "ImplicitEuler;")
		{
			parameters.solvingDEMethod = SolvingDEMethod::ImplicitEuler;
		}

		else if (solvingDEMethod == "RungeKutta4;")
		{
			parameters.solvingDEMethod = SolvingDEMethod::RungeKutta4;
		}

		else
		{
			throw std::string("Unknown method for solving DEs");
		}

		initialFile.getline(tempString, 32, '=');
		std::string calculateWithCUDA;
		initialFile >> calculateWithCUDA;

		parameters.calculateWithCUDA = calculateWithCUDA == "true;";

		initialFile.close();
	}

	void readInitialVelocitiesAndCoordinatates(const std::string path)
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

	void defineCalculationParameters()
	{
		switch (parameters.solvingDEMethod)
		{
		case Euler:
			integrator = new EulerIntegrator();
			break;
		case ImplicitEuler:
			integrator = new ImplicitEulerIntegrator();
			break;
		case RungeKutta4:
			integrator = new RungeKutta4Integrator();
			break;
		default:
			throw std::string("Unknown method for solving DEs");
		}

		//calculator = parameters.calculateWithCUDA ? new AccelerationCalculatorCUDA() : new SequentialAccelerationCalculator();

		if (parameters.calculateWithCUDA)
		{
			calculator = new AccelerationCalculatorCUDA();
		}
		else
		{
			calculator = new SequentialAccelerationCalculator();
		}
	}

public:
	NBodySystem() : 
		particles{ nullptr },
		n{ 0 },
		parameters{ DTOParameters() },
		integrator{ nullptr },
		calculator{ nullptr },
		momentum { new Momentum() },
		energy { new Energy() }
	{}

	void initializeNBodySystem()
	{
		setInitialParameters("initialParameters.txt");

		readInitialVelocitiesAndCoordinatates("Particles.txt");

		defineCalculationParameters();
	}

	void writeInFileAllPositions(const long count)
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

	void initializeSystemWithLog(std::ofstream& logFile)
	{
		std::cout << "Reading input file\n";

		auto start = std::chrono::high_resolution_clock::now();

		initializeNBodySystem();

		auto finish = std::chrono::high_resolution_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);

		std::cout << "File is read in " << ms.count() << " milliseconds\n\n\n";
		logFile << "File is read in " << ms.count() << " milliseconds\n\n\n";
	}

	bool createCoordinatesDirectory()
	{
		std::filesystem::path path = L"coordinates";
		if (std::filesystem::exists(path))
		{
			std::filesystem::remove_all(path);
		}

		return std::filesystem::create_directory(path);
	}

	void writeToFileWithLog(const long count, std::ofstream& logFile)
	{
		if (count % parameters.outputFrequency != 0)
			return;

		std::cout << "\nWriting to file all positions\n";

		auto start = std::chrono::high_resolution_clock::now();

		writeInFileAllPositions(count);

		auto finish = std::chrono::high_resolution_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);

		std::cout << "Writing to file has finished in " << ms.count() << " milliseconds\n\n";
		logFile << "Writing to file has finished in " << ms.count() << " milliseconds\n\n";
	}

	void calculateVelocitiesAndPositionsWithLog(std::ofstream& logFile)
	{
		std::cout << "\nCalculating velocities and positions\n";

		auto start = std::chrono::high_resolution_clock::now();

		integrator->integrate(particles, n, parameters, *calculator);

		auto finish = std::chrono::high_resolution_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);

		std::cout << "Calculating velocities and positions has finished in " << ms.count() << " milliseconds\n\n\n";
		logFile << "Calculating velocities and positions has finished in " << ms.count() << " milliseconds\n\n\n";
	}

	void checkLawsOfConservation(const long count)
	{
		if (count % parameters.physicalCharacteristicsOutputFrequency != 0)
			return;

		momentum->checkLawOfConservation(particles, n);
		energy->checkLawOfConservation(particles, n);
	}

	void simulate()
	{
		std::ofstream logFile;
		logFile.open("log.txt");

		initializeSystemWithLog(logFile);

		if (!createCoordinatesDirectory())
		{
			throw std::string{ "Error making a directory\n" };
		}

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

		delete[] particles;
		delete integrator;
		delete calculator;
		delete momentum;
		delete energy;

		logFile.close();
	}
};