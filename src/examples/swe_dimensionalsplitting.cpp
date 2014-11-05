/*
 * swedimensionalsplitting.h
 *
 *  Created on: Nov 5, 2014
 *      Author: kyu
 */

#ifndef SWEDIMENSIONALSPLITTING_CPP_
#define SWEDIMENSIONALSPLITTING_CPP_

#include <limits>
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

		for (int i = 1; i < nx+2; i++) {
	    	for (int j=1; j < ny+1; ++j) {
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
 		if (maxWaveSpeed > zeroTol) {
	    	maxTimestep = 0.4f * std::min (dx / maxWaveSpeed, dy / maxWaveSpeed);
	    } else {
	    	//might happen in dry cells
	    	maxTimestep = std::numeric_limits<float>::max ();
	    }
	};

void swe_dimensionalsplitting::updateUnknowns(float dt){
	  	for (int i = 1; i < nx+1; i++) {
	  		for (int j = 1; j < ny + 1; j++) {
	   			//Formula (8)
	   			h[i][j] -= dt / dx * (hRight[i - 1][j - 1] + hLeft[i][j - 1]) + dt / dy * (hAbove[i - 1][j - 1] + hBelow[i - 1][j]);
	   			hu[i][j] -= dt / dx * (huRight[i - 1][j - 1] + huLeft[i][j - 1]);
	   			hv[i][j] -= dt / dy * (hvAbove[i - 1][j - 1] + hvBelow[i - 1][j]);

	   			if (h[i][j] < dryTol) {
	   				h[i][j] = hu[i][j] = hv[i][j] = 0.0f;
	   			} else if (h[i][j] < 0.1f)
	   				hu[i][j] = hv[i][j] = 0.0f;
	   			}
	   		}
	};


#endif /* SWEDIMENSIONALSPLITTING_CPP_ */
