/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *         Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Univ.-Prof._Dr._Michael_Bader)
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * Basic setting of SWE, which uses a wave propagation solver and an artificial or ASAGI scenario on a single block.
 */

#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>


#include "blocks/swe_dimensionalsplitting.hh"
#include "scenarios/SWE_Scenario.hh"
#include "scenarios/SWE_TsunamiScenario.hh"
#include "scenarios/SWE_ArtificialTsunamiScenario.hh"
#include "scenarios/SWE_Checkpoint.hh"

//#ifndef CUDA
//#include "blocks/SWE_WavePropagationBlock.hh"
//#else
//#include "blocks/cuda/SWE_WavePropagationBlockCuda.hh"
//#endif

#ifdef WRITENETCDF
#include "writer/NetCdfWriter.hh"
#else
#include "writer/VtkWriter.hh"
#endif

#ifdef ASAGI
#include "scenarios/SWE_AsagiScenario.hh"
#else
#include "scenarios/SWE_simple_scenarios.hh"
#endif

#ifdef READXML
#include "tools/CXMLConfig.hpp"
#endif

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
  #ifndef READXML
  args.addOption("grid-size-x", 'x', "Number of cells in x direction");
  args.addOption("grid-size-y", 'y', "Number of cells in y direction");
  args.addOption("simulation-time", 't', "Number of seconds of simulation", args.Required ,false);
  args.addOption("boundary-condition", 'b', "1: OUTFLOW ,2:WALL ,3:INFLOW, 4:CONNECT ,5:PASSIVE", args.Required, false);
  args.addOption("output-basepath", 'o', "Output base file name");
  #endif

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

  int l_time, l_boundary;

  //! l_baseName of the plots.
  std::string l_baseName;

  // read command line parameters
  #ifndef READXML
  l_nX = args.getArgument<int>("grid-size-x");
  l_nY = args.getArgument<int>("grid-size-y");
  l_time = args.getArgument<int>("simulation-time", 0);
  l_baseName = args.getArgument<std::string>("output-basepath");
  #endif

  // read xml file
  #ifdef READXML
  assert(false); //TODO: not implemented.
  if(argc != 2) {
    s_sweLogger.printString("Aborting. Please provide a proper input file.");
    s_sweLogger.printString("Example: ./SWE_gnu_debug_none_augrie config.xml");
    return 1;
  }
  s_sweLogger.printString("Reading xml-file.");

  std::string l_xmlFile = std::string(argv[1]);
  s_sweLogger.printString(l_xmlFile);

  CXMLConfig l_xmlConfig;
  l_xmlConfig.loadConfig(l_xmlFile.c_str());
  #endif

  #ifdef ASAGI
  /* Information about the example bathymetry grid (tohoku_gebco_ucsb3_500m_hawaii_bath.nc):
   *
   * Pixel node registration used [Cartesian grid]
   * Grid file format: nf = GMT netCDF format (float)  (COARDS-compliant)
   * x_min: -500000 x_max: 6500000 x_inc: 500 name: x nx: 14000
   * y_min: -2500000 y_max: 1500000 y_inc: 500 name: y ny: 8000
   * z_min: -6.48760175705 z_max: 16.1780223846 name: z
   * scale_factor: 1 add_offset: 0
   * mean: 0.00217145586762 stdev: 0.245563641735 rms: 0.245573241263
   */

  //simulation area
  float simulationArea[4];
  simulationArea[0] = -450000;
  simulationArea[1] = 6450000;
  simulationArea[2] = -2450000;
  simulationArea[3] = 1450000;

  SWE_AsagiScenario l_scenario( ASAGI_INPUT_DIR "tohoku_gebco_ucsb3_500m_hawaii_bath.nc",
                                ASAGI_INPUT_DIR "tohoku_gebco_ucsb3_500m_hawaii_displ.nc",
                                (float) 28800., simulationArea);
  #else

  SWE_Scenario *s;


  bool isCheckpointScenario = false;
  int switchScenario = 2;     //edit this to switch between scenarios
  switch(switchScenario){
  case 0:{
	  s = new SWE_RadialDamBreakScenario;}break;
  case 1:{
	  s = new SWE_ArtificialTsunamiScenario;;}break;
  case 2:{
	  s = new SWE_TsunamiScenario;;}break;
  case 3:{
	  s = new SWE_CheckpointScenario;;
	  isCheckpointScenario = true;}break;

  }

 SWE_Scenario &l_scenario = *s;


  #endif

  if(isCheckpointScenario) {
	  l_nX = l_scenario.getxDim();
	  l_nY = l_scenario.getyDim();
  }

  //! size of a single cell in x- and y-direction
  float l_dX, l_dY;

  // compute the size of a single cell
  l_dX = (l_scenario.getBoundaryPos(BND_RIGHT) - l_scenario.getBoundaryPos(BND_LEFT) )/l_nX;
  l_dY = (l_scenario.getBoundaryPos(BND_TOP) - l_scenario.getBoundaryPos(BND_BOTTOM) )/l_nY;

  // create a single wave propagation block
  #ifndef CUDA
  swe_dimensionalsplitting l_wavePropgationBlock(l_nX,l_nY,l_dX,l_dY);
  #else
  SWE_WavePropagationBlockCuda l_wavePropgationBlock(l_nX,l_nY,l_dX,l_dY);
  #endif

  //! origin of the simulation domain in x- and y-direction
  float l_originX, l_originY;

  // get the origin from the scenario
  l_originX = l_scenario.getBoundaryPos(BND_LEFT);
  l_originY = l_scenario.getBoundaryPos(BND_BOTTOM);

  // initialize the wave propagation block
  l_wavePropgationBlock.initScenario(l_originX, l_originY, l_scenario);


  //! time when the simulation ends.
  float l_endSimulation = l_time > 0 ? l_time : l_scenario.endSimulation();

  //! number of checkpoints for visualization (at each checkpoint in time, an output file is written).
    int l_numberOfCheckPoints = 20;
    if(isCheckpointScenario) {
    	/*
    	cout << "Calculated Steps: " << l_scenario.calculatedSteps() << endl;
    	cout << "Continue At: " << l_scenario.continueSimulationAt() << endl;
    	cout << "alpha: " << (l_scenario.calculatedSteps()/l_scenario.continueSimulationAt()) << endl;
    	*/
  	    l_numberOfCheckPoints = (int) ((l_endSimulation - l_scenario.continueSimulationAt())*(l_scenario.calculatedSteps()/l_scenario.continueSimulationAt()));
    }

  //! checkpoints when output files are written.
  float* l_checkPoints = new float[l_numberOfCheckPoints+1];

  // compute the checkpoints in time
  for(int cp = 0; cp <= l_numberOfCheckPoints; cp++) {
     l_checkPoints[cp] = l_scenario.continueSimulationAt() + cp*((l_endSimulation-l_scenario.continueSimulationAt())/l_numberOfCheckPoints);
  }

  // Init fancy progressbar
  tools::ProgressBar progressBar(l_endSimulation);

  // write the output at time zero
  tools::Logger::logger.printOutputTime((float) 0.);
  progressBar.update(0.);

  std::string l_fileName = generateBaseFileName(l_baseName,0,0);
  //boundary size of the ghost layers
  io::BoundarySize l_boundarySize = {{1, 1, 1, 1}};

  //boundary type
  #ifndef READXML
  l_boundary = args.getArgument<int>("boundary-condition", -1);
  #endif
  BoundaryType l_boundaryType[4];
  if(l_boundary == -1){
	  l_boundaryType[0] = l_scenario.getBoundaryType(BND_LEFT);
	  l_boundaryType[1] = l_scenario.getBoundaryType(BND_RIGHT);
	  l_boundaryType[2] = l_scenario.getBoundaryType(BND_BOTTOM);
	  l_boundaryType[3] = l_scenario.getBoundaryType(BND_TOP);
  }else{
	  BoundaryType b;
	  switch(l_boundary){
	  case 2: b = WALL; break;
	  case 3: b = INFLOW; break;
	  case 4: b = CONNECT; break;
	  case 5: b = PASSIVE; break;
	  default: b = OUTFLOW;
	  }
	  l_boundaryType[0] = l_boundaryType[1] = l_boundaryType[2] = l_boundaryType[3] = b;
	  l_wavePropgationBlock.setBoundaryType(BND_LEFT, b);
	  l_wavePropgationBlock.setBoundaryType(BND_RIGHT, b);
	  l_wavePropgationBlock.setBoundaryType(BND_BOTTOM, b);
	  l_wavePropgationBlock.setBoundaryType(BND_TOP, b);
  }

  float l_boundaryPos[4];
  l_boundaryPos[0] = l_scenario.getBoundaryPos(BND_TOP);
  l_boundaryPos[1] = l_scenario.getBoundaryPos(BND_BOTTOM);
  l_boundaryPos[2] = l_scenario.getBoundaryPos(BND_LEFT);
  l_boundaryPos[3] = l_scenario.getBoundaryPos(BND_RIGHT);


#ifdef WRITENETCDF
  //construct a NetCdfWriter
  io::NetCdfWriter l_writer( l_fileName,
		  l_wavePropgationBlock.getBathymetry(),
		  l_boundarySize,
		  l_boundaryType,
		  l_boundaryPos,
		  l_endSimulation,
		  isCheckpointScenario,
		  l_nX, l_nY,
		  l_dX, l_dY,
		  l_originX, l_originY,
		  1);
#else
  // consturct a VtkWriter
  io::VtkWriter l_writer( l_fileName,
		  l_wavePropgationBlock.getBathymetry(),
		  l_boundarySize,
		  l_nX, l_nY,
		  l_dX, l_dY );
#endif

  if(!isCheckpointScenario){
  // Write zero time step
  l_writer.writeTimeStep( l_wavePropgationBlock.getWaterHeight(),
                          l_wavePropgationBlock.getDischarge_hu(),
                          l_wavePropgationBlock.getDischarge_hv(),
                          (float) 0.);
  }

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
      l_wavePropgationBlock.setGhostLayer();

      // reset the cpu clock
      tools::Logger::logger.resetClockToCurrentTime("Cpu");


      // compute numerical flux on each edge
      l_wavePropgationBlock.computeNumericalFluxes();

      //! maximum allowed time step width.
      float l_maxTimeStepWidth = l_wavePropgationBlock.getMaxTimestep();

      // update the cell values
      l_wavePropgationBlock.updateUnknowns(l_maxTimeStepWidth);

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

    // write output
    l_writer.writeTimeStep( l_wavePropgationBlock.getWaterHeight(),
                            l_wavePropgationBlock.getDischarge_hu(),
                            l_wavePropgationBlock.getDischarge_hv(),
                            l_t);
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
