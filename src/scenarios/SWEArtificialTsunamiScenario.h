/*
 * SWEArtificialTsunamiScenario.h
 *
 *  Created on: Nov 13, 2014
 *      Author: kyu
 */

#ifndef SWEARTIFICIALTSUNAMISCENARIO_H_
#define SWEARTIFICIALTSUNAMISCENARIO_H_

#include "SWE_Scenario.hh"

class SWE_ArtificialTsunamiScenario: public SWE_Scenario {
public:
	SWE_ArtificialTsunamiScenario();
	float getWaterHeight(float x, float y);

	float getBathymetry(float x, float y);

};

#endif /* SWEARTIFICIALTSUNAMISCENARIO_H_ */
