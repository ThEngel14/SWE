
#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>


#include "blocks/swe_dimensionalsplitting.hh"
#include "blocks/SWE_WavePropagationBlock.cpp"
#include "scenarios/SWE_Scenario.hh"
//#include "scenarios/SWE_TsunamiScenario.hh"
#include "scenarios/SWE_ArtificialTsunamiScenario.hh"
//#include "scenarios/SWE_Checkpoint.hh"


/*
#ifdef WRITENETCDF
#include "writer/NetCdfWriter.hh"
#else
#include "writer/VtkWriter.hh"
#endif
*/

#include "writer/Writer.hh"


#include "scenarios/SWE_simple_scenarios.hh"


#include "tools/args.hh"
#include "tools/help.hh"
#include "tools/Logger.hh"
#include "tools/ProgressBar.hh"

/**
 * Main program for the simulation on a single SWE_WavePropagationBlock.
 */
int main( int argc, char** argv ) {
  /**
   * Initialization.
   */
  // Parse command line parameters
  tools::Args args;

  args.addOption("grid-size-x", 'x', "Number of cells in x direction");
  args.addOption("grid-size-y", 'y', "Number of cells in y direction");
  args.addOption("simulation-time", 't', "Number of seconds of simulation", args.Required ,false);
  args.addOption("boundary-condition", 'b', "0: OUTFLOW ,1:WALL ,2:INFLOW, 3:CONNECT ,4:PASSIVE", args.Required, false);
  args.addOption("output-basepath", 'o', "Output base file name");
  args.addOption("scenario", 's', "0: TsunamiScenario, 1: CheckPointsScenario ,2: ArtificialTsunamiScenario ,3:RadialDambreakScenario", args.Required, false);


  tools::Args::Result ret = args.parse(argc, argv);

  switch (ret)
  {
  case tools::Args::Error:
	  return 1;
  case tools::Args::Help:
	  return 0;
  }

  //! number of grid cells in x- and y-direction.
  int l_nX, l_nY;

  int l_time, l_boundary, l_scen;

  //! l_baseName of the plots.
  std::string l_baseName;

  //scenario
  l_scen = args.getArgument<int>("scenario", 0);


  SWE_Scenario *s;


  bool isCheckpointScenario = false;
  switch(l_scen){
  case 1://{
	  //s = new SWE_CheckpointScenario("_00.nc");
	  //isCheckpointScenario = true;}break;
	  exit(-1);
  case 2:{
	  s = new SWE_ArtificialTsunamiScenario;}break;
  case 3:{
	  s = new SWE_RadialDamBreakScenario;}break;
  default:
	  //s = new SWE_TsunamiScenario;

  }

 SWE_Scenario &l_scenario = *s;

 // read command line parameters
 if(isCheckpointScenario) {
	  l_nX = l_scenario.getxDim();
	  l_nY = l_scenario.getyDim();
 }else{
	  l_nX = args.getArgument<int>("grid-size-x");
	  l_nY = args.getArgument<int>("grid-size-y");
 }
 l_time = args.getArgument<int>("simulation-time", l_scenario.endSimulation());
 //boundary type
 l_boundary = args.getArgument<int>("boundary-condition", -1);
 l_baseName = args.getArgument<std::string>("output-basepath");


  //! size of a single cell in x- and y-direction
  float l_dX, l_dY;

  // compute the size of a single cell
  l_dX = (l_scenario.getBoundaryPos(BND_RIGHT) - l_scenario.getBoundaryPos(BND_LEFT) )/l_nX;
  l_dY = (l_scenario.getBoundaryPos(BND_TOP) - l_scenario.getBoundaryPos(BND_BOTTOM) )/l_nY;

  // create a single wave propagation block
  swe_dimensionalsplitting l_dimensionalsplitting(l_nX,l_nY,l_dX,l_dY);


  //! origin of the simulation domain in x- and y-direction
  float l_originX, l_originY;

  // get the origin from the scenario
  l_originX = l_scenario.getBoundaryPos(BND_LEFT);
  l_originY = l_scenario.getBoundaryPos(BND_BOTTOM);

  // initialize the wave propagation block
  l_dimensionalsplitting.initScenario(l_originX, l_originY, l_scenario);

  //! number of checkpoints for visualization (at each checkpoint in time, an output file is written).
    int l_numberOfCheckPoints = 20;
    if(isCheckpointScenario) {
    	/*
    	cout << "Calculated Steps: " << l_scenario.calculatedSteps() << endl;
    	cout << "Continue At: " << l_scenario.continueSimulationAt() << endl;
    	cout << "alpha: " << (l_scenario.calculatedSteps()/l_scenario.continueSimulationAt()) << endl;
    	*/
  	    l_numberOfCheckPoints = (int) ((l_time - l_scenario.continueSimulationAt())*(l_scenario.calculatedSteps()/l_scenario.continueSimulationAt()));
    }

  //! checkpoints when output files are written.
  float* l_checkPoints = new float[l_numberOfCheckPoints+1];

  // compute the checkpoints in time
  for(int cp = 0; cp <= l_numberOfCheckPoints; cp++) {
     l_checkPoints[cp] = l_scenario.continueSimulationAt() + cp*((l_time-l_scenario.continueSimulationAt())/l_numberOfCheckPoints);
  }

  // Init fancy progressbar
  tools::ProgressBar progressBar(l_time);

  // write the output at time zero
  tools::Logger::logger.printOutputTime((float) 0.);
  progressBar.update(0.);

  std::string l_fileName = generateBaseFileName(l_baseName,0,0);
  //boundary size of the ghost layers
  io::BoundarySize l_boundarySize = {{1, 1, 1, 1}};



  BoundaryType l_boundaryType[4];
  if(l_boundary == -1){
	  l_boundaryType[0] = l_scenario.getBoundaryType(BND_LEFT);
	  l_boundaryType[1] = l_scenario.getBoundaryType(BND_RIGHT);
	  l_boundaryType[2] = l_scenario.getBoundaryType(BND_BOTTOM);
	  l_boundaryType[3] = l_scenario.getBoundaryType(BND_TOP);
  }else{
	  BoundaryType b;
	  switch(l_boundary){
	  case 1: b = WALL; break;
	  case 2: b = INFLOW; break;
	  case 3: b = CONNECT; break;
	  case 4: b = PASSIVE; break;
	  default: b = OUTFLOW;
	  }
	  l_boundaryType[0] = l_boundaryType[1] = l_boundaryType[2] = l_boundaryType[3] = b;
	  l_dimensionalsplitting.setBoundaryType(BND_LEFT, b);
	  l_dimensionalsplitting.setBoundaryType(BND_RIGHT, b);
	  l_dimensionalsplitting.setBoundaryType(BND_BOTTOM, b);
	  l_dimensionalsplitting.setBoundaryType(BND_TOP, b);
  }

  float l_boundaryPos[4];
  l_boundaryPos[0] = l_scenario.getBoundaryPos(BND_TOP);
  l_boundaryPos[1] = l_scenario.getBoundaryPos(BND_BOTTOM);
  l_boundaryPos[2] = l_scenario.getBoundaryPos(BND_LEFT);
  l_boundaryPos[3] = l_scenario.getBoundaryPos(BND_RIGHT);


/*
#ifdef WRITENETCDF
  //construct a NetCdfWriter
  io::NetCdfWriter l_writer( l_fileName,
		  l_dimensionalsplitting.getBathymetry(),
		  l_boundarySize,
		  l_boundaryType,
		  l_boundaryPos,
		  l_time,
		  isCheckpointScenario,
		  l_nX, l_nY,
		  l_dX, l_dY,
		  l_originX, l_originY,
		  1);
#else
  // consturct a VtkWriter
  io::VtkWriter l_writer( l_fileName,
		  l_dimensionalsplitting.getBathymetry(),
		  l_boundarySize,
		  l_nX, l_nY,
		  l_dX, l_dY );
#endif
*/

  /*
  if(!isCheckpointScenario){
  // Write zero time step
  l_writer.writeTimeStep( l_dimensionalsplitting.getWaterHeight(),
                          l_dimensionalsplitting.getDischarge_hu(),
                          l_dimensionalsplitting.getDischarge_hv(),
                          (float) 0.);
  }
  */

  /**
   * Simulation.
   */
  // print the start message and reset the wall clock time
  progressBar.clear();
  tools::Logger::logger.printStartMessage();
  tools::Logger::logger.initWallClockTime(time(NULL));

  //! simulation time.
  float l_t = 0.0;
  int beginCount = 1;

  //iterate to the next checkpoint time marker
  if(isCheckpointScenario) {
	  l_t = l_scenario.continueSimulationAt();
	  int i=0;
	  while(l_t > l_checkPoints[i]){
		  i++;
	  }
	  beginCount = i;
  }

  progressBar.update(l_t);

  unsigned int l_iterations = 0;

  // loop over checkpoints
  for(int c=beginCount; c<=l_numberOfCheckPoints; c++) {

    // do time steps until next checkpoint is reached
    while( l_t < l_checkPoints[c] ) {
      // set values in ghost cells:
      l_dimensionalsplitting.setGhostLayer();

      // reset the cpu clock
      tools::Logger::logger.resetClockToCurrentTime("Cpu");


      // compute numerical flux on each edge
      l_dimensionalsplitting.computeNumericalFluxes();

      //! maximum allowed time step width.
      float l_maxTimeStepWidth = l_dimensionalsplitting.getMaxTimestep();

      // update the cell values
      l_dimensionalsplitting.updateUnknowns(l_maxTimeStepWidth);

      // update the cpu time in the logger
      tools::Logger::logger.updateTime("Cpu");

      // update simulation time with time step width.
      l_t += l_maxTimeStepWidth;
      l_iterations++;

      // print the current simulation time
      progressBar.clear();
      tools::Logger::logger.printSimulationTime(l_t);
      progressBar.update(l_t);
    }

    // print current simulation time of the output
    progressBar.clear();
    tools::Logger::logger.printOutputTime(l_t);
    progressBar.update(l_t);

    //cout << "**********************\n*****************" << endl;

    /*
    // write output
    l_writer.writeTimeStep( l_dimensionalsplitting.getWaterHeight(),
                            l_dimensionalsplitting.getDischarge_hu(),
                            l_dimensionalsplitting.getDischarge_hv(),
                            l_t);
    */
  }

  /**
   * Finalize.
   */
  // write the statistics message
  progressBar.clear();
  tools::Logger::logger.printStatisticsMessage();
  // print the cpu time
  tools::Logger::logger.printTime("Cpu", "CPU time");
  // print the wall clock time (includes plotting)
  tools::Logger::logger.printWallClockTime(time(NULL));
  // printer iteration counter
  tools::Logger::logger.printIterationsDone(l_iterations);
  return 0;
}
