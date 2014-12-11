/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
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
 * A writer for the netCDF-format: http://www.unidata.ucar.edu/software/netcdf/
 */

#include "NetCdfWriter.hh"
#include <string>
#include <vector>
#include <iostream>
#include <cassert>

#define CHECKPOINT_FILE "_00.nc"

/**
 * Create a netCdf-file
 * Any existing file will be replaced.
 *
 * @param i_baseName base name of the netCDF-file to which the data will be written to.
 * @param i_nX number of cells in the horizontal direction.
 * @param i_nY number of cells in the vertical direction.
 * @param i_dX cell size in x-direction.
 * @param i_dY cell size in y-direction.
 * @param i_originX
 * @param i_originY
 * @param i_flush If > 0, flush data to disk every i_flush write operation
 * @param i_dynamicBathymetry
 */
io::NetCdfWriter::NetCdfWriter( const std::string &i_baseName,
		const Float2D &i_b,
		const BoundarySize &i_boundarySize,
		const enum BoundaryType *i_boundaryType,
		const float *i_boundaryPos,
		const float endSimulation,
		bool isCheckPoint,
		int delta,
		int i_nX, int i_nY,
		float i_dX, float i_dY,
		float i_originX, float i_originY,
		unsigned int i_flush) :
		//const bool  &i_dynamicBathymetry) : //!TODO
  io::Writer(i_baseName + ".nc", i_b, i_boundarySize, i_nX, i_nY),
  flush(i_flush)
{

	checkPoint = isCheckPoint;

	if(isCheckPoint) {
		int status;
		status = nc_open(CHECKPOINT_FILE, NC_WRITE, &dataFile);

		//check if the netCDF-file creation constructor succeeded.
		if (status != NC_NOERR) {
			assert(false);
			return;
		}

		status = nc_inq_varid(dataFile, "time", &timeVar);
		status = nc_inq_varid(dataFile, "h", &hVar);
		status = nc_inq_varid(dataFile, "hu", &huVar);
		status = nc_inq_varid(dataFile, "hv", &hvVar);
		status = nc_inq_varid(dataFile, "Boundary", &boundaryVar);
		status = nc_inq_varid(dataFile, "BoundaryPos", &boundaryPosVar);
		status = nc_inq_varid(dataFile, "EndSimulation", &endSimulationVar);

		int timedimid;
		size_t timeDim;
		status = nc_inq_dimid(dataFile, "time", &timedimid);
		status = nc_inq_dimlen(dataFile, timedimid, &timeDim);

		timeStep = timeDim;

	} else {
		int status;

		deltaX = delta;

		//create a netCDF-file, an existing file will be replaced
		status = nc_create(fileName.c_str(), NC_NETCDF4, &dataFile);

		//check if the netCDF-file creation constructor succeeded.
		if (status != NC_NOERR) {
			assert(false);
			return;
		}

	#ifdef PRINT_NETCDFWRITER_INFORMATION
		std::cout << "   *** io::NetCdfWriter::createNetCdfFile" << std::endl;
		std::cout << "     created/replaced: " << fileName << std::endl;
		std::cout << "     dimensions(nx, ny): " << nX << ", " << nY << std::endl;
		std::cout << "     cell width(dx,dy): " << i_dX << ", " << i_dY << std::endl;
		std::cout << "     origin(x,y): " << i_originX << ", " << i_originY << std::endl;
	#endif

		// Calculate cell size
		ndX_single = (int) ((1.0f * nX)/deltaX);
		ndX_extended = (nX % deltaX != 0 ? 1 : 0);
		ndX = ndX_single + ndX_extended;

		ndY_single = (int) ((1.0f * nY)/deltaX);
		ndY_extended = (nY % deltaX != 0 ? 1 : 0);
		ndY = ndY_single + ndY_extended;

		//dimensions
		int l_timeDim, l_xDim, l_yDim, l_boundaryDim, l_boundaryPosDim, l_endSimulationDim;
		nc_def_dim(dataFile, "time", NC_UNLIMITED, &l_timeDim);
		nc_def_dim(dataFile, "x", ndX, &l_xDim);
		nc_def_dim(dataFile, "y", ndY, &l_yDim);
		nc_def_dim(dataFile, "EndSimulation", 1, &l_endSimulationDim);
		nc_def_dim(dataFile, "Boundary", 4, &l_boundaryDim);
		nc_def_dim(dataFile, "BoundaryPos", 4, &l_boundaryPosDim);

		//variables (TODO: add rest of CF-1.5)
		int l_xVar, l_yVar;

		nc_def_var(dataFile, "time", NC_FLOAT, 1, &l_timeDim, &timeVar);
		ncPutAttText(timeVar, "long_name", "Time");
		ncPutAttText(timeVar, "units", "seconds since simulation start"); // the word "since" is important for the paraview reader

		nc_def_var(dataFile, "x", NC_FLOAT, 1, &l_xDim, &l_xVar);
		nc_def_var(dataFile, "y", NC_FLOAT, 1, &l_yDim, &l_yVar);

		//variables, fastest changing index is on the right (C syntax), will be mirrored by the library
		int dims[] = {l_timeDim, l_yDim, l_xDim};
		int boundarydim[] = {l_boundaryDim};
		int boundaryPosdim[] = {l_boundaryPosDim};
		int esdim[] = {l_endSimulationDim};
		nc_def_var(dataFile, "h",  NC_FLOAT, 3, dims, &hVar);
		nc_def_var(dataFile, "hu", NC_FLOAT, 3, dims, &huVar);
		nc_def_var(dataFile, "hv", NC_FLOAT, 3, dims, &hvVar);
		nc_def_var(dataFile, "b",  NC_FLOAT, 2, &dims[1], &bVar);
		nc_def_var(dataFile, "Boundary", NC_INT, 1, boundarydim, &boundaryVar);
		nc_def_var(dataFile, "BoundaryPos", NC_FLOAT, 1,boundaryPosdim, &boundaryPosVar);
		nc_def_var(dataFile, "EndSimulation", NC_FLOAT, 1, esdim, &endSimulationVar);


		//set attributes to match CF-1.5 convention
		ncPutAttText(NC_GLOBAL, "Conventions", "CF-1.5");
		ncPutAttText(NC_GLOBAL, "title", "Computed tsunami solution");
		ncPutAttText(NC_GLOBAL, "history", "SWE");
		ncPutAttText(NC_GLOBAL, "institution", "Technische Universitaet Muenchen, Department of Informatics, Chair of Scientific Computing");
		ncPutAttText(NC_GLOBAL, "source", "Bathymetry and displacement data.");
		ncPutAttText(NC_GLOBAL, "references", "http://www5.in.tum.de/SWE");
		ncPutAttText(NC_GLOBAL, "comment", "SWE is free software and licensed under the GNU General Public License. Remark: In general this does not hold for the used input data.");

		//setup grid size
		float gridPosition = i_originX;
		float offset = 0;
		size_t counter = 0;
		for(size_t i = 1; i <= nX; i++) {
			offset += i_dX;

			if(i % deltaX == 0) {
				gridPosition += offset / 2.0f;

				nc_put_var1_float(dataFile, l_xVar, &counter, &gridPosition);
				gridPosition += offset / 2.0f;
				offset = 0;
				counter++;
			}
		}
		if(counter < ndX) {
			gridPosition += offset / 2.0f;
			nc_put_var1_float(dataFile, l_xVar, &counter, &gridPosition);
		}

		gridPosition = i_originY;
		offset = 0;
		counter = 0;
		for(size_t j = 1; j <= nY; j++) {
			offset += i_dY;

			if(j % deltaX == 0) {
				gridPosition += offset / 2.0f;

				nc_put_var1_float(dataFile, l_yVar, &counter, &gridPosition);
				gridPosition += offset / 2.0f;
				offset = 0;
				counter++;
			}

		}
		if(counter < ndY) {
			gridPosition += offset / 2.0f;
			nc_put_var1_float(dataFile, l_yVar, &counter, &gridPosition);
		}

		//setup boundary type
		int boundary[4] = {0,0,0,0};
		for(int i = 0; i < 4; i++) {
			switch (i_boundaryType[i]) {
			case OUTFLOW: boundary[i] = 1; break;
			case WALL: boundary[i] = 2; break;
			case INFLOW: boundary[i] = 3; break;
			case CONNECT: boundary[i] = 4; break;
			case PASSIVE: boundary[i] = 5; break;
			}
		}

		const size_t startbound[] = {0};
		const size_t countbound[] = {4};
		const size_t countend[] = {1};
		nc_put_vara_int(dataFile, boundaryVar,startbound, countbound, boundary);

		nc_put_vara_float(dataFile, boundaryPosVar, startbound, countbound,i_boundaryPos);

		nc_put_vara_float(dataFile, endSimulationVar,startbound, countend, &endSimulation);
	}
}

/**
 * Destructor of a netCDF-writer.
 */
io::NetCdfWriter::~NetCdfWriter() {
	nc_close(dataFile);
}

/**
 * Writes time dependent data to a netCDF-file (-> constructor) with respect to the boundary sizes.
 *
 * boundarySize[0] == left
 * boundarySize[1] == right
 * boundarySize[2] == bottom
 * boundarySize[3] == top
 *
 * @param i_matrix array which contains time dependent data.
 * @param i_boundarySize size of the boundaries.
 * @param i_ncVariable time dependent netCDF-variable to which the output is written to.
 */
void io::NetCdfWriter::writeVarTimeDependent( const Float2D &i_matrix,
                                              int i_ncVariable ) {
	//write col wise, necessary to get rid of the boundary
	//storage in Float2D is col wise
	//read carefully, the dimensions are confusing
	size_t start[] = {timeStep, 0, 0};
	size_t count[] = {1, ndY, 1};


	for(unsigned int col = 0; col < ndX; col++) {
		start[2] = col; //select col (dim "x")

		float temp[ndY];
		for(int k = 0; k < ndY; k++) {
			temp[k] = getAverage(i_matrix, col*deltaX + boundarySize[0], k*deltaX + boundarySize[2]);
		}


		nc_put_vara_float(dataFile, i_ncVariable, start, count,
				&temp[0]); //write col
	}
}

/**
 * Write time independent data to a netCDF-file (-> constructor) with respect to the boundary sizes.
 * Variable is time-independent
 * boundarySize[0] == left
 * boundarySize[1] == right
 * boundarySize[2] == bottom
 * boundarySize[3] == top
 *
 * @param i_matrix array which contains time independent data.
 * @param i_boundarySize size of the boundaries.
 * @param i_ncVariable time independent netCDF-variable to which the output is written to.
 */
void io::NetCdfWriter::writeVarTimeIndependent( const Float2D &i_matrix,
                                                int i_ncVariable ) {
	//write col wise, necessary to get rid of the boundary
	//storage in Float2D is col wise
	//read carefully, the dimensions are confusing
	size_t start[] = {0, 0};
	size_t count[] = {ndY, 1};

	for(unsigned int col = 0; col < ndX; col++) {
		start[1] = col; //select col (dim "x")

		float temp[ndY];
		for(int k = 0; k < ndY; k++) {
			temp[k] = getAverage(i_matrix, col*deltaX + boundarySize[0], k*deltaX + boundarySize[2]);
		}


		nc_put_vara_float(dataFile, i_ncVariable, start, count,
				&temp[0]); //write col
	}
}

/**
 * Calculates and returns the average cell value of the given matrix for the rectangle ndX x ndY starting at the given coordinates
 *
 * @param matrix
 * @param x_start
 * @param y_start
 *
 * @return the average value
 */
float io::NetCdfWriter::getAverage(const Float2D &matrix, int x_start, int y_start) {
	float value = 0;
	int cell_count = 0;

	for(int i = x_start; i < x_start + ndX_single && i <= nX; i++) {
		for(int k = y_start; k < y_start + ndY_single && k <= nY; k++) {
			cell_count++;
			value += matrix[i][k];
		}
	}


	return value / cell_count;
}

/**
 * Writes the unknwons to a netCDF-file (-> constructor) with respect to the boundary sizes.
 *
 * boundarySize[0] == left
 * boundarySize[1] == right
 * boundarySize[2] == bottom
 * boundarySize[3] == top
 *
 * @param i_h water heights at a given time step.
 * @param i_hu momentums in x-direction at a given time step.
 * @param i_hv momentums in y-direction at a given time step.
 * @param i_boundarySize size of the boundaries.
 * @param i_time simulation time of the time step.
 */
void io::NetCdfWriter::writeTimeStep( const Float2D &i_h,
                                      const Float2D &i_hu,
                                      const Float2D &i_hv,
                                      float i_time) {
	if (timeStep == 0 && !checkPoint) {
		// Write bathymetry
		writeVarTimeIndependent(b, bVar);
	}

	/*
	printf("Timestep: %d\n", timeStep);
	printf("xDim: %d yDim: %d\n", nX, nY);
	printf("DataFile: %d\n", dataFile);
	*/

	//write i_time
	nc_put_var1_float(dataFile, timeVar, &timeStep, &i_time);

	//write water height
	writeVarTimeDependent(i_h, hVar);

	//write momentum in x-direction
	writeVarTimeDependent(i_hu, huVar);

	//write momentum in y-direction
	writeVarTimeDependent(i_hv, hvVar);

	// Increment timeStep for next call
	timeStep++;

	if (flush > 0 && timeStep % flush == 0) {
		nc_sync(dataFile);
	}
}

