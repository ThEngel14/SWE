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

void swe_dimensionalsplitting::computeNumericalFluxes()
{

	float maxWaveSpeed_x = 0.0f;
	float maxWaveSpeed_y = 0.0f;
	float maxEdgeSpeed = 0.0f;

	#ifdef OPENMP
		float t_maxWaveSpeed_x = 0.0f;
		float t_maxWaveSpeed_y = 0.0f;

		#pragma omp parallel firstprivate(maxEdgeSpeed, t_maxWaveSpeed_x)
		{
			#pragma omp for
	#endif

	//x-sweep (vertical edges)
	for(int i = 1; i < nx+2 ; i++)
	{
		for (int j = 0; j < ny + 2; j++)
		{
			fwave.computeNetUpdates (
				h[i - 1][j], h[i][j],
				hu[i - 1][j], hu[i][j],
				b[i - 1][j], b[i][j],
				hLeft[i - 1][j - 1], hRight[i - 1][j - 1],
				huLeft[i - 1][j - 1], huRight[i - 1][j - 1],
				maxEdgeSpeed
			);

			// Update maxWaveSpeed
			#ifdef OPENMP
				t_maxWaveSpeed_x = std::max(t_maxWaveSpeed_x,maxEdgeSpeed);
			#else
				maxWaveSpeed_x = std::max(maxWaveSpeed_x,maxEdgeSpeed);
			#endif
		}
	}

	//update max wave speed in x-direction from all threads
	#ifdef OPENMP
			if(maxWaveSpeed_x< t_maxWaveSpeed_x)
			{
				#pragma omp critical
				{
					maxWaveSpeed_x = t_maxWaveSpeed_x;
				}
			}
		}
	#endif




	maxTimestep = 0.4f * (dx / maxWaveSpeed_x);


	//update height here
	#ifdef OPENMP
		#pragma omp parallel for
	#endif
	for (int i = 1; i < nx+1; i++) {
		for (int j=1; j < ny+1; j++) {
			h[i][j] -= (maxTimestep / dx) * (hRight[i - 1][j - 1] + hLeft[i][j-1]);
			hu[i][j] -= (maxTimestep / dx) * (huRight[i - 1][j - 1] + huLeft[i][j-1]);

			if(b[i][j]>0)
			{
				h[i][j]=0.0f;
				hv[i][j]=0.0f;
			}
		}
	}

	maxEdgeSpeed = 0.0f;

	//y-sweep (horizontal edges)
	#ifdef OPENMP
		#pragma omp parallel firstprivate(maxEdgeSpeed, t_maxWaveSpeed_y)
		{
			#pragma omp for
	#endif
	for(int i = 1; i < nx+1 ; i++)
	{
		for (int j = 1; j < ny + 2; j++)
		{
			fwave.computeNetUpdates (
				h[i][j - 1], h[i][j],
				hv[i][j - 1], hv[i][j],
				b[i][j - 1], b[i][j],
				hBelow[i - 1][j - 1], hAbove[i - 1][j - 1],
				hvBelow[i - 1][j - 1], hvAbove[i - 1][j - 1],
				maxEdgeSpeed
			);

		// udpate max wave speed
		#ifdef OPENMP
			t_maxWaveSpeed_y = std::max(t_maxWaveSpeed_y,maxEdgeSpeed);
		#else
			maxWaveSpeed_y = std::max(maxWaveSpeed_y,maxEdgeSpeed);
		#endif
		}
	}

	// update max wave speed in y-direction from all threads
	#ifdef OPENMP
			if(maxWaveSpeed_y< t_maxWaveSpeed_y)
			{
				#pragma omp critical
				{
					maxWaveSpeed_y = t_maxWaveSpeed_y;
				}
			}
		}
	#endif

#ifndef NDEBUG
	//Check CFL-condition for y-sweep
    if(maxTimestep >= 0.5 * dy/maxWaveSpeed_y){
    	std::cerr << "CFL-condition for y-sweep is not satisfied" << std::endl;
    }
#endif // NDEBUG


}

void swe_dimensionalsplitting::updateUnknowns(float dt)
{
	#ifdef OPENMP
		#pragma omp parallel for
	#endif
	for(int i=1;i<nx+1;i++)
	{
		for(int j=1;j<ny+1;j++)
		{
			h[i][j] -= (dt / dy) * (hAbove[i-1][j - 1] + hBelow[i-1][j]);
			hv[i][j] -= (dt / dy) * (hvAbove[i-1][j - 1] + hvBelow[i-1][j]);

			if(b[i][j]>0)
			{
				h[i][j]=0.0f;
				hv[i][j]=0.0f;
			}
		}

	}
}

#endif /* SWEDIMENSIONALSPLITTING_CPP_ */
