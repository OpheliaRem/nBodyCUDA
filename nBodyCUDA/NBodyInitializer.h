#pragma once
#include "FileReader.h"
#include "FileWithParticlesReader.h"
#include "FileWithInputParametersReader.h"
#include "Integrator.h"
#include "EulerIntegrator.h"
#include "PredictorCorrectorIntegrator.h"
#include "RungeKutta2Integrator.h"
#include "RungeKutta4Integrator.h"
#include "AccelerationCalculator.h"
#include "SequentialAccelerationCalculator.h"
#include "AccelerationCalculatorCUDA.h"
#include <filesystem>

class NBodyInitializer
{
	FileReader* particlesReader;

	FileReader* parametersReader;

	DTOParameters& parameters;

	Integrator*& integrator;

	AccelerationCalculator*& calculator;

	void defineCalculationParameters()
	{
		switch (parameters.solvingDEMethod)
		{
		case Euler:
			integrator = new EulerIntegrator();
			break;
		case PredictorCorrector:
			integrator = new PredictorCorrectorIntegrator();
			break;
		case RungeKutta2:
			integrator = new RungeKutta2Integrator();
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

	bool createCoordinatesDirectory()
	{
		std::filesystem::path path = L"coordinates";
		if (std::filesystem::exists(path))
		{
			std::filesystem::remove_all(path);
		}

		return std::filesystem::create_directory(path);
	}

public:
	NBodyInitializer(
		const std::string particlePath,
		const std::string parametersPath,
		Particle*& particles,
		int& n,
		DTOParameters& parameters,
		Integrator*& integrator,
		AccelerationCalculator*& calculator
	) : 
		particlesReader{ new FileWithParticlesReader(particlePath, particles, n) },
		parametersReader{ new FileWithInputParametersReader(parametersPath, parameters) },
		parameters{parameters},
		integrator{integrator},
		calculator{calculator}
	{}

	void initialize()
	{
		particlesReader->readFile();
		parametersReader->readFile();

		defineCalculationParameters();

		if (!createCoordinatesDirectory())
		{
			throw std::string{ "Error making a directory \"coordinates\"\n" };
		}
	}
};