/*
 * swedimensionalsplitting.hh
 *
 *  Created on: Nov 5, 2014
 *      Author: kyu
 */

#ifndef SWEDIMENSIONALSPLITTING_HH_
#define SWEDIMENSIONALSPLITTING_HH_

#include "tools/help.hh"
#include "scenarios/SWE_Scenario.hh"
#include "blocks/SWE_Block.hh"
#include "solvers/FWave.hpp"
#include "scenarios/SWE_simple_scenarios.hh"

class swe_dimensionalsplitting :public SWE_Block{
private:
		const float zeroTol;
		const float dryTol;

		solver::FWave<float> fwave;
		Float2D hLeft;
	    Float2D hRight;
	    Float2D huLeft;
	    Float2D huRight;
	    Float2D hBelow;
	    Float2D hAbove;
	    Float2D hvBelow;
	    Float2D hvAbove;

public:

	swe_dimensionalsplitting(int l_nx, int l_ny, float l_dx, float l_dy);

	/**
	 * Not needed for this case.
	 */
	void computeNumericalFluxes();

	/**
	 * Compute netUpdates for the given input-values using x-and y-sweep
	 * @param dt time step for the update
	 */
	void updateUnknowns(float dt);

};

#endif /* SWEDIMENSIONALSPLITTING_HH_ */
