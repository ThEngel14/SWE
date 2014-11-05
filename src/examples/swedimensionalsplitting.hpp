/*
 * swedimensionalsplitting.h
 *
 *  Created on: Nov 5, 2014
 *      Author: kyu
 */

#ifndef SWEDIMENSIONALSPLITTING_HPP_
#define SWEDIMENSIONALSPLITTING_HPP_

#include "../scenarios/SWE_Scenario.hh"
#include "../blocks/SWE_Block.hh"
#include "../solvers/FWave.hpp"
#include "../scenarios/SWE_simple_scenarios.hh"

class swe_dimensionalsplitting : SWE_Block{
public:
	swe_dimensionalsplitting(int l_nx, int l_ny, float l_dx, float l_dy)
		:SWE_Block(l_nx, l_ny, l_dx, l_dy){
		SWE_RadialDamBreakScenario scenario();

		initScenario(0.0f, 0.0f, scenario, false);



	}

	virtual void updateUnknowns(float dt){

	}
	void computeXSweep(){

	}
	void computeYSweep(){

	}

};

#endif /* SWEDIMENSIONALSPLITTING_HPP_ */
