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
#include "../examples/swe_dimensionalsplitting.cpp"
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

	float getWaterHeight(float x, float y) {
		return 50.0f;
	}

	float getDischargeU(float x, float y) {
		return 0.0f;
	}

	float getDischargeV(float x, float y) {
		return 0.0f;
	}

	void testDimensionalSplitting() {
		swe_dimensionalsplitting dimSplitting(3, 3, 0.5f, 0.5f);

		float (*_h)(float, float) = getWaterHeight;
		float (*_hu)(float, float) = getDischargeU;
		float (*_hv)(float, float) = getDischargeV;

		dimSplitting.setWaterHeight(_h);
		dimSplitting.setBathymetry(-10);
		dimSplitting.setDischarge(_hu, _hv);

		dimSplitting.computeNumericalFluxes();
		dimSplitting.updateUnknowns(dimSplitting.getMaxTimestep());

		float *h;
		float *hu;
		float *hv;

		h = dimSplitting.getWaterHeight();
		hu = dimSplitting.getDischarge_hu();
		hv = dimSplitting.getDischarge_hv();

		float maxSpeed = 0.0f;
		float uphl[3];
		float uphr[3];
		float uphul[3];
		float uphur[3];

		// x-sweep
		for(int i = 0; i < 2; i++) {
			for(int k = 0; k < 3; k++) {
				fwave.computeNetUpdates(getWaterHeight(i, k), getWaterHeight(i+1, k), zero, zero, uphl[k], uphr[k], uphul[k], uphur[k], maxSpeed);
			}
		}

		// y-sweep
		float up2hl[3];
		float up2hr[3];
		float up2hul[3];
		float up2hur[3];
		for(int i = 0; i < 2; i++) {
			for(int k = 0; k < 3; k++) {
				fwave.computeNetUpdates(getWaterHeight(k, i), getWaterHeight(k, i+1), zero, zero, uphl[k], uphr[k], uphul[k], uphur[k], maxSpeed);
			}
		}

		float h_res[3][3];
		float hu_res[3][3];
		float hv_res[3][3];
		float dt = dimSplitting.getMaxTimestep();
		float dx = dimSplitting.dx;
		float dy = dimSplitting.dy;

		for (int i = 1; i < 3; i++) {
			for (int k = 1; k < 3; k++) {
				h_res[i][k] -= dt / dx * (uphr[i - 1][k - 1] + uphr[i][k - 1]) + dt / dy * (up2hl[i - 1][k - 1] + up2hr[i - 1][k]);
				hu_res[i][k] -= dt / dx * (uphur[i - 1][k - 1] + uphul[i][k - 1]);
				hv_res[i][k] -= dt / dy * (up2hul[i - 1][k - 1] + up2hur[i - 1][k]);
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
