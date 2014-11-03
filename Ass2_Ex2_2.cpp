/**
 * This is a scenario to solve Exercise 2.1 of Assignment 2 
 * SUBCRITICAL FLOW
 */
class SWE_RadialDamBreakScenario : public SWE_Scenario {

  public:

  	// get the bathymerty at Position x 


  	// Formula (6)
    float getBathymetry(float x) {
    	if(x <= 12 && x >= 8)
       		return (-1.8-0.05*((x-10)^2)) ;
       	else
       		return -2;
    };

    // get the waterheight at Position x

    float getWaterHeight(float x) { 

      	 if(x <= 12 && x >= 8)
       		return (-(-1.8-0.05*((x-10)^2))) ;
       	else
       		return 2;

    };


    // get the momentum at Position x

	 float getMomentum(float x) {

	 	return 4.42;

	 }


// get the number of timesteps for the current scenario

	 float getTimesteps() {

	 	return 200;

	 }


    
   