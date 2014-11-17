/*
 * SWE_Checkpoint.hh
 *
 *  Created on: 17.11.2014
 *      Author: vmuser
 */

#ifndef SWE_CHECKPOINTSCENARIO_HH_
#define SWE_CHECKPOINTSCENARIO_HH_

#include "SWE_Scenario.hh"
#include <cmath>
#include <limits>
#include <netcdf.h>

#define CHECKPOINT_FILE "_00.nc"

class SWE_CheckpointScenario: public SWE_Scenario {
private:
	float *bathymetry;
	float *water;
	float *hu;
	float *hv;
	int* boundary;
	/**
	 * Print error message.
	 */
	void handle_error(int status) {
		if (status != NC_NOERR) {
			cout<< nc_strerror(status)<<"\n";
			exit(-1);
		}
	};

public:
	/**
	 * Open and preparing a set of nc-files for bathymetry and displacement data.
	 */
	SWE_CheckpointScenario(){
		int ncid;
		int status;
		status = nc_open(CHECKPOINT_FILE, NC_NOWRITE, ncid);
		if(status != NC_NOERR){handle_error(status);}

		// get Dimensions
		int timeDim, xDim, yDim, boundaryDim;

		status = nc_inq_dimlen(*ncid, 0, timeDim);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_dimlen(*ncid, 1, xDim);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_dimlen(*ncid, 2, yDim);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_dimlen(*ncid, 3, boundaryDim);
		if (status != NC_NOERR) handle_error(status);

		// get Ids
		int hId, huId, hvId, bId, boundaryId;

		status = nc_inq_varid(*ncid, "h", &hId);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_varid(*ncid, "hu", &huId);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_varid(*ncid, "hv", &hvId);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_varid(*ncid, "b", &bId);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_varid(*ncid, "Boundary", &boundaryId);
		if (status != NC_NOERR) handle_error(status);

		// init arrays
		int *boundary = (int*) malloc(boundaryDim*sizeof(int));
		float **water = (float*) malloc(xDim*sizeof(float));
		float **hu = (float*) malloc(xDim*sizeof(float));
		float **hv = (float*) malloc(xDim*sizeof(float));
		float **bathymetry = (float*) malloc(xDim*sizeof(float));

		for(int k = 0; k < xDim; k++) {
			water[k] = (float*) malloc(yDim*sizeof(float));
			hu[k] = (float*) malloc(yDim*sizeof(float));
			hv[k] = (float*) malloc(yDim*sizeof(float));
			bathymetry[k] = (float*) malloc(yDim*sizeof(float));
		}

		// get data
		status = nc_get_var_int(ncid, boundaryId, &boundary[0]);
		if(status != NC_NOERR){handle_error(status);}

		status = nc_get_var_float(ncid, bId, &bathymetry[0][0]);
		if(status != NC_NOERR) {handle_error(status);}

		const size_t start = {timeDim-1, 0, 0};
		const size_t count = {timeDim-1, xDim, yDim};

		status = nc_get_vara_float(ncid, hId, start, count, &water);
		if(status != NC_NOERR){handle_error(status);}

		status = nc_get_vara_float(ncid, huId, start, count, &hu);
		if(status != NC_NOERR){handle_error(status);}

		status = nc_get_vara_float(ncid, hvId, start, count, &hv);
		if(status != NC_NOERR){handle_error(status);}

		// close file
		status = nc_close(ncid);
		if (status != NC_NOERR) handle_error(status);
	};

	float getWaterHeight(float x, float y){
		return water[x][y];
	};

	float getBathymetry(float x, float y){
		return bathymetry[x][y];
	};

	float getVeloc_u(float x, float y) {
		return hu[x][y]/water[x][y];
	}

	float getVeloc_v(float x, float y) {
		return hv[x][y]/water[x][y];
	}

	virtual float endSimulation() {
		//TODO - have to write this into the checkpoint file
		return (float) 15;
	};

	/**
	 *
	 * @param edge BoundaryEdge to request the BoundaryType
	 * @return the BoundaryType for the given edge
	 */
	virtual BoundaryType getBoundaryType(BoundaryEdge edge) {
		if(edge == BND_TOP) {
			switch (boundary[0]) {
			case 1: return OUTFLOW;
			case 2: return WALL;
			case 3: return INFLOW;
			case 4: return CONNECT;
			case 5: return PASSIVE;
			}
		}

		if(edge == BND_BOTTOM) {
			switch (boundary[1]) {
			case 1: return OUTFLOW;
			case 2: return WALL;
			case 3: return INFLOW;
			case 4: return CONNECT;
			case 5: return PASSIVE;
			}
		}

		if(edge == BND_LEFT) {
			switch (boundary[2]) {
			case 1: return OUTFLOW;
			case 2: return WALL;
			case 3: return INFLOW;
			case 4: return CONNECT;
			case 5: return PASSIVE;
			}
		}

		if(edge == BND_RIGHT) {
			switch (boundary[3]) {
			case 1: return OUTFLOW;
			case 2: return WALL;
			case 3: return INFLOW;
			case 4: return CONNECT;
			case 5: return PASSIVE;
			}
		}

		return OUTFLOW;
	};

	/**
	 * Get the boundary positions
	 *
	 * @param i_edge which edge
	 * @return value in the corresponding dimension
	 */
	float getBoundaryPos(BoundaryEdge i_edge) {
		//TODO
		return 0;
	};
};


#endif /* SWE_CHECKPOINTSCENARIO_HH_ */
