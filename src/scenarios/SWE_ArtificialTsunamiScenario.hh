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

	float computeDisplacement(float x, float y){
		if(std::fabs(x) <= 500.0f && std::fabs(y) <= 500.0f){
			float dx = std::sin(((x/500.0f)+1)*M_PI);
			float dy = - ((y/500)*(y/500.0f)) + 1;
			return 5.0f  * dx * dy;
		}
		return 0.0f;
	}

	float getBathymetryBefore(float x, float y){
		float bath = -100.0f;
		return bath;
	}
public:
	SWE_ArtificialTsunamiScenario(){

	};

	float getWaterHeight(float x, float y){
		float bath = getBathymetryBefore(x,y);
		float absBath = std::fabs(bath);
		if( absBath < 20.0f){
			bath = bath / absBath * 20.0f;
		}
		return -std::min(bath,0.0f);
	};

	float getBathymetry(float x, float y){
		float bath = getBathymetryBefore(x,y) + computeDisplacement(x,y);

		float absBath = std::fabs(bath);
		if( absBath < 20.0f){
			bath = bath / absBath * 20.0f;
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
	      return (float)-5000;
	    else if ( i_edge == BND_RIGHT)
	      return (float)5000 ;
	    else if ( i_edge == BND_BOTTOM )
	      return (float)-5000 ;
	    else
	      return (float)5000 ;
	};
};

#endif /* SWEARTIFICIALTSUNAMISCENARIO_HH_ */
