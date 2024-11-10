#pragma once
#include "Particle.h"

class PhysicalCharacteristic
{
public:
	virtual void checkLawOfConservation(Particle* particles, const int n) = 0;
};