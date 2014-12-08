/*
 * CheckpointingTest.h
 *
 *  Created on: 08.12.2014
 *      Author: vmuser
 */

#ifndef CHECKPOINTINGTEST_H_
#define CHECKPOINTINGTEST_H_

#include <cxxtest/TestSuite.h>
#include "../scenarios/SWE_Checkpoint.hh"
#include "../scenarios/SWE_simple_scenarios.hh"
#include <netcdf.h>
#include <cstdlib>

class CheckpointingTest : public CxxTest::TestSuite
{
public:
	void testSameHeightNoSpeed(void) {
		SWE_CheckpointScenario scenario("unittest.nc");
		SWE_RadialDamBreakScenario dambreak;

		float l_dx, l_dy;

		//dimensions
		int xDim = scenario.getxDim();
		int yDim = scenario.getyDim();

		TS_ASSERT_EQUALS(xDim, 5);
		TS_ASSERT_EQUALS(yDim, 10);

		//Boundary Pos
		float left, right, top, bottom;
		left = scenario.getBoundaryPos(BND_LEFT);
		right = scenario.getBoundaryPos(BND_RIGHT);
		top = scenario.getBoundaryPos(BND_TOP);
		bottom = scenario.getBoundaryPos(BND_BOTTOM);

		TS_ASSERT_EQUALS(left, dambreak.getBoundaryPos(BND_LEFT));
		TS_ASSERT_EQUALS(right, dambreak.getBoundaryPos(BND_RIGHT));
		TS_ASSERT_EQUALS(top, dambreak.getBoundaryPos(BND_TOP));
		TS_ASSERT_EQUALS(bottom, dambreak.getBoundaryPos(BND_BOTTOM));

		l_dx = (right-left)/xDim;
		l_dy = (top - bottom)/yDim;

		//Boundary Type
		int s_left = scenario.getBoundaryType(BND_LEFT);
		int d_left;
		switch (dambreak.getBoundaryType(BND_LEFT)) {
		case OUTFLOW: d_left = 1; break;
		case WALL: d_left = 2; break;
		case INFLOW: d_left = 3; break;
		case CONNECT: d_left = 4; break;
		case PASSIVE: d_left = 5; break;
		}

		TS_ASSERT_EQUALS(s_left, d_left);

		//Bathymetry, hu, hv, water
		for(int i = 0; i < xDim; i++) {
			for(int j = 0; j < yDim; j++) {
				float ax, ay;
				ax = l_dx*i;
				ay = l_dy*j;

				TS_ASSERT_EQUALS(scenario.getWaterHeight(ax, ay), dambreak.getWaterHeight(ax, ay));
				TS_ASSERT_EQUALS(scenario.getBathymetry(ax, ay), dambreak.getBathymetry(ax, ay));
				TS_ASSERT_EQUALS(scenario.getVeloc_u(ax, ay), dambreak.getVeloc_u(ax, ay));
				TS_ASSERT_EQUALS(scenario.getVeloc_v(ax, ay), dambreak.getVeloc_v(ax, ay));
			}
		}

		//other
		TS_ASSERT_EQUALS(scenario.endSimulation(), 15);
		TS_ASSERT_EQUALS(scenario.calculatedSteps(), 3);
		TS_ASSERT_DELTA(scenario.continueSimulationAt(), 21, 0.5);

	}
};


#endif /* CHECKPOINTINGTEST_H_ */
