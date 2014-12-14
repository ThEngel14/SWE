/*
 * swedimensionalsplitting.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: kyu
 */

#ifndef SWEDIMENSIONALSPLITTING_CPP_
#define SWEDIMENSIONALSPLITTING_CPP_

#include <cassert>
#include <string>
#include <limits>
#include <omp.h>
#include "swe_dimensionalsplitting.hh"

swe_dimensionalsplitting::swe_dimensionalsplitting(int l_nx, int l_ny, float l_dx, float l_dy)
		:SWE_Block(l_nx, l_ny, l_dx, l_dy),
		zeroTol(0.000001f),
		dryTol(0.01f),
		hLeft (nx + 1, ny),
		hRight (nx + 1, ny),
		huLeft (nx + 1, ny),
		huRight (nx + 1, ny),

		hBelow (nx, ny + 1),
		hAbove (nx, ny + 1),
		hvBelow (nx, ny + 1),
		hvAbove (nx, ny + 1)
		{};


void swe_dimensionalsplitting::computeNumericalFluxes(){
	float maxWaveSpeed = (float) 0.0f;

//x-sweep (vertical edges)

#pragma omp parallel shared(maxWaveSpeed, h, hu, b, hLeft, huLeft, hRight, huRight)
	{
	#pragma omp for
		{
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
						#pragma omp critical
						{
						maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
						}
				}
			}
		}
	}


	float dt;
	if(maxWaveSpeed > zeroTol){
		dt = 0.4 *(dx/maxWaveSpeed);
	}else{
		dt = std::numeric_limits<float>::max();
	}


#pragma omp parallel shared(h, hu, hLeft, huLeft, hRight, huRight)
	{
	#pragma omp for
		{
			for (int i = 1; i < nx+1; i++) {
				for (int j=1; j < ny+1; j++) {
					h[i][j] -= (dt / dx) * (hRight[i - 1][j - 1] + hLeft[i][j-1]);
					hu[i][j] -= (dt / dx) * (huRight[i - 1][j - 1] + huLeft[i][j-1]);
				}
			}
		}
	}


//y-sweep (horizontal edges)
#pragma omp parallel shared(maxWaveSpeed, h, hu, b, hvBelow, huLeft, hAbove, hvAbove)
	{
		#pragma omp for
		{
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
						#pragma omp critical
						{
						maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
						}
				}
			}
		}
	}



#pragma omp parallel shared(h, hu, hAbove, hvAbove, hBelow, hvBelow)
	{
	#pragma omp for
		{
			for (int i = 1; i < nx+1; i++) {
				for (int j=1; j < ny+1; j++) {
					h[i][j] -= (dt / dy) * (hAbove[i-1][j - 1] + hBelow[i-1][j]);
					hv[i][j] -= (dt / dy) * (hvAbove[i-1][j - 1] + hvBelow[i-1][j]);
				}
			}
		}
	}


	maxTimestep = dt;

#ifndef NDEBUG
	//Check CFL-condition for y-sweep
    if(maxTimestep >= 0.5 * dy/maxWaveSpeed){
    	std::cerr << "CFL-condition for y-sweep is not satisfied" << std::endl;
    }
#endif // NDEBUG

};

void swe_dimensionalsplitting::updateUnknowns(float dt){

	};

#endif /* SWEDIMENSIONALSPLITTING_CPP_ */
