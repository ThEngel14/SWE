/*
 * FWaveTest.h
 *
 *  Created on: 20.10.2014
 *      Author: vmuser
 */

#ifndef DIMENSIONALSPLITTINGTEST_H_
#define DIMENSIONALSPLITTINGTEST_H_

#include <cxxtest/TestSuite.h>
#include "../solvers/FWave.hpp"
#include "../blocks/swe_dimensionalsplitting.cpp"
#include "../blocks/SWE_Block.cpp"

class DimensionalSplittingTest : public CxxTest::TestSuite
{
private:
	/** Zero: 0.0f */
	float zero;

public:
	/** Delta used for comparing expected and actual values */
	static const float delta = 0.0001f;

	/** FWave to test */
	solver::FWave<float> fwave;

	DimensionalSplittingTest() {
		zero = 0.0f;
	}

	void testDimensionalSplitting() {
		swe_dimensionalsplitting dimSplitting(50, 50, 0.5f, 0.5f);
		SWE_RadialDamBreakScenario scenario;

		dimSplitting.initScenario(50,50,scenario);

		int dx = (scenario.getBoundaryPos(BND_RIGHT) - scenario.getBoundaryPos(BND_LEFT) )/50;
		int dy = (scenario.getBoundaryPos(BND_TOP) - scenario.getBoundaryPos(BND_BOTTOM) )/50;

		dimSplitting.computeNumericalFluxes();
		dimSplitting.updateUnknowns(dimSplitting.getMaxTimestep());

		float *h;
		float *hu;
		float *hv;

		h = dimSplitting.getWaterHeight();
		hu = dimSplitting.getDischarge_hu();
		hv = dimSplitting.getDischarge_hv();

		float h_res[50][50];
		float hu_res[50][50];
		float hv_res[50][50];
		float b[50][50];

		for(int i = 0; i < 50; i++) {
			for(int j = 0; j < 50; j++) {
				h_res[i][j] = scenario.getWaterHeight(scenario.getBoundaryPos(BND_LEFT) + i*dx, scenario.getBoundaryPos(BND_BOTTOM)+j*dy);
				b[i][j] = scenario.getBathymetry(scenario.getBoundaryPos(BND_LEFT) + i*dx, scenario.getBoundaryPos(BND_BOTTOM)+j*dy);
				hu_res[i][j] = scenario.getVeloc_u(scenario.getBoundaryPos(BND_LEFT) + i*dx, scenario.getBoundaryPos(BND_BOTTOM)+j*dy);
				hu_res[i][j] *= h_res[i][j];
				hv_res[i][j] = scenario.getVeloc_v(scenario.getBoundaryPos(BND_LEFT) + i*dx, scenario.getBoundaryPos(BND_BOTTOM)+j*dy);
				hv_res[i][j] *= h_res[i][j];
			}
		}

		float maxWaveSpeed;

		float hLeft[52][52];
		float hRight[52][52];
		float huLeft[52][52];
		float huRight[52][52];

		float hBelow[52][52];
		float hAbove[52][52];
		float hvBelow[52][52];
		float hvAbove[52][52];

		float maxSpeed = 0.0f;

		int nx = 50;
		int ny = 50;
		//x-sweep (vertical edges)
		for (int i = 1; i < nx+2; i++) {
			for (int j=1; j < ny+1; j++) {
				float maxEdgeSpeed;
					fwave.computeNetUpdates (
					h[i - 1][j], h[i][j],
					hu[i - 1][j], hu[i][j],
					b[i - 1][j], b[i][j],
					hLeft[i - 1][j - 1], hRight[i - 1][j - 1],
					huLeft[i - 1][j - 1], huRight[i - 1][j - 1],
					maxEdgeSpeed
				);
					maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
			}
		}

		float dt = 0.4 *(dx/maxWaveSpeed);
		for (int i = 1; i < nx+1; i++) {
			for (int j=1; j < ny+2; j++) {
				h[i][j] -= dt / dx * (hRight[i - 1][j - 1] + hLeft[i][j-1]);
				hu[i][j] -= dt / dx * (huRight[i - 1][j-1] + huLeft[i][j-1]);
			}
		}


		//y-sweep (horizontal edges)
		for (int i = 1; i < nx+1; i++) {
			for (int j=1; j < ny+2; ++j) {
				float maxEdgeSpeed;

				fwave.computeNetUpdates (
					h[i][j - 1], h[i][j],
					hv[i][j - 1], hv[i][j],
					b[i][j - 1], b[i][j],
					hBelow[i - 1][j - 1], hAbove[i - 1][j - 1],
					hvBelow[i - 1][j - 1], hvAbove[i - 1][j - 1],
					maxEdgeSpeed
				);

					maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
			}
		}

		dt = 0.4*(dy / maxWaveSpeed);
		for (int i = 1; i < nx+1; i++) {
			for (int j=1; j < ny+2; j++) {
				h[i][j] -= dt / dy * (hAbove[i-1][j - 1] + hBelow[i-1][j]);
				hv[i][j] -= dt / dy * (hvAbove[i-1][j - 1] + hvBelow[i-1][j]);
			}
		}

		//Asserts
		for (int i = 0; i < 3; i++) {
			for (int k = 0; k < 3; k++) {
				TS_ASSERT_DELTA(h[i][k], h_res[i][k], delta);
				TS_ASSERT_DELTA(hu[i][k], hu_res[i][k], delta);
				TS_ASSERT_DELTA(hv[i][k], hv_res[i][k], delta);
			}
		}
	}
};


#endif /* DIMENSIONALSPLITTINGTEST_H_ */
