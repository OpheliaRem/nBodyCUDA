#pragma once
#include "FileReader.h"
#include "DTOParameters.h"
#include <iostream>
#include <fstream>

class FileWithInputParametersReader : public FileReader
{
	const std::string path;
	DTOParameters& parameters;

public:
	FileWithInputParametersReader(
		const std::string path,
		DTOParameters& parameters
	) : path{path}, parameters{parameters} {}

	void readFile() override
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

		else if (solvingDEMethod == "PredictorCorrector;")
		{
			parameters.solvingDEMethod = SolvingDEMethod::PredictorCorrector;
		}

		else if (solvingDEMethod == "RungeKutta2;")
		{
			parameters.solvingDEMethod = SolvingDEMethod::RungeKutta2;
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
};