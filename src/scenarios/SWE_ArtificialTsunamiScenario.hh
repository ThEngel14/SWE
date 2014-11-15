/*
 * SWEArtificialTsunamiScenario.hh
 *
 *  Created on: Nov 13, 2014
 *      Author: kyu
 */

#ifndef SWEARTIFICIALTSUNAMISCENARIO_HH_
#define SWEARTIFICIALTSUNAMISCENARIO_HH_

#include "SWE_Scenario.hh"
#include <cmath>

class SWE_ArtificialTsunamiScenario: public SWE_Scenario {
private:
	static const float scale = 1/100.0f; //use the denominator as cell number

	float computeDisplacement(float x, float y){
		float fhm = 500*scale; 			//FiveHundertMeter * scale
		if(std::fabs(x) <= fhm && std::fabs(y) <=fhm){
			float dx = std::sin(((x/fhm)+1)*M_PI);
			float dy = - ((y/fhm)*(y/fhm)) + 1;
			return 5 * scale * dx * dy;
		}
		return 0.0f;
	}
public:
	SWE_ArtificialTsunamiScenario(){

	};

	float getWaterHeight(float x, float y){
		return 100.0f * scale + computeDisplacement(x,y);
	};

	float getBathymetry(float x, float y){
		float bath = -100.0f * scale + computeDisplacement(x,y);
		if(bath >= -20 *scale && bath <= 20 *scale){
			if(bath > 0){
				return 20 *scale;
			}else{
				return -20*scale;
			}
		}
		return bath;
	};
	virtual float endSimulation() { return (float) 15; };

	/**
	 * @param edge BoundaryEdge to request the BoundaryType
	 *
	 * @return the BoundaryType for the given edge
	 */
	   virtual BoundaryType getBoundaryType(BoundaryEdge edge) { return WALL; };

	/**
	 * Get the boundary positions
	 *
	 * @param i_edge which edge
	 * @return value in the corresponding dimension
	 */
	float getBoundaryPos(BoundaryEdge i_edge) {
	    if ( i_edge == BND_LEFT )
	      return (float)-5000 * scale;
	    else if ( i_edge == BND_RIGHT)
	      return (float)5000 * scale;
	    else if ( i_edge == BND_BOTTOM )
	      return (float)-5000 * scale;
	    else
	      return (float)5000 * scale;
	};
};

#endif /* SWEARTIFICIALTSUNAMISCENARIO_HH_ */
