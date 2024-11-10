#pragma once
#include "SolvingDEMethod.h"

class DTOParameters
{
public:
	double timeStep;
	double limitOfLoop;
	int outputFrequency;
	int physicalCharacteristicsOutputFrequency;
	SolvingDEMethod solvingDEMethod;
	bool calculateWithCUDA;
};