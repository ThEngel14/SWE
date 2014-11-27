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
		if(std::fabs(x) <= 500 && std::fabs(y) <=500){
			float dx = std::sin(((x/500)+1)*M_PI);
			float dy = - ((y/500)*(y/500)) + 1;
			return 5  * dx * dy;
		}
		return 0.0f;
	}
public:
	SWE_ArtificialTsunamiScenario(){

	};

	float getWaterHeight(float x, float y){
		return 100.0f; //-std::min(-100.0f,0.0f);
	};

	float getBathymetry(float x, float y){
		float bath = -100.0f + computeDisplacement(x,y);
		if(bath >= -20  && bath <= 20 ){
			if(bath > 0){
				return 20 ;
			}else{
				return -20;
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
