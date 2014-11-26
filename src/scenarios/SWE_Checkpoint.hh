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
	float* boundaryPos;
	float l_endSimulation;
	float l_time;
	size_t xDim, yDim, timeDim;
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
		int ncid, timedimid, xdimid, ydimid, boundaryid, boundaryPosid, endSimulationid;
		int status;
		status = nc_open(CHECKPOINT_FILE, NC_NOWRITE, &ncid);
		if(status != NC_NOERR){handle_error(status);}

//######## get Dimensions
		size_t boundaryDim, boundaryPosDim, endSimulationDim;

		status = nc_inq_dimid(ncid, "time", &timedimid);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_dimid(ncid, "x", &xdimid);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_dimid(ncid, "y", &ydimid);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_dimid(ncid, "Boundary", &boundaryid);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_dimid(ncid, "BoundaryPos", &boundaryPosid);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_dimid(ncid, "EndSimulation", &endSimulationid);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_dimlen(ncid, timedimid, &timeDim);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_dimlen(ncid, xdimid, &xDim);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_dimlen(ncid, ydimid, &yDim);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_dimlen(ncid, boundaryid, &boundaryDim);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_dimlen(ncid, boundaryPosid, &boundaryPosDim);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_dimlen(ncid, endSimulationid, &endSimulationDim);
		if (status != NC_NOERR) handle_error(status);

//######### get Ids
		int hId, huId, hvId, bId, boundaryId, boundaryPosId, endSimulationId, timeId;

		status = nc_inq_varid(ncid, "h", &hId);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_varid(ncid, "hu", &huId);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_varid(ncid, "hv", &hvId);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_varid(ncid, "b", &bId);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_varid(ncid, "Boundary", &boundaryId);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_varid(ncid, "BoundaryPos", &boundaryPosId);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_varid(ncid, "EndSimulation", &endSimulationId);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_varid(ncid, "time", &timeId);
		if (status != NC_NOERR) handle_error(status);

		// init arrays
		boundary = (int*) malloc(4*sizeof(int));
		boundaryPos = (float*) malloc(4*sizeof(float));
		water = (float*) malloc(xDim*yDim*sizeof(float));
		hu = (float*) malloc(xDim*yDim*sizeof(float));
		hv = (float*) malloc(xDim*yDim*sizeof(float));
		bathymetry = (float*) malloc(xDim*yDim*sizeof(float));

		// get data

		const size_t startboundary[] = {0};
		const size_t countboundary[] = {4};
		status = nc_get_vara_int(ncid, boundaryId,startboundary, countboundary, boundary);
		if(status != NC_NOERR){handle_error(status);}

		status = nc_get_vara_float(ncid, boundaryPosId,startboundary, countboundary, boundaryPos);
		if(status != NC_NOERR){handle_error(status);}

		const size_t esIndex[] = {0};
		status = nc_get_var1_float(ncid, endSimulationId, esIndex, &l_endSimulation);
		if(status != NC_NOERR){handle_error(status);}

		const size_t tIndex[] = {timeDim-1};
		status = nc_get_var1_float(ncid, timeId, tIndex , &l_time);
		if(status != NC_NOERR){handle_error(status);}

		const size_t startbathymetry[] = {0 ,0};
		const size_t countbathymetry[] = {yDim, xDim};
		status = nc_get_vara_float(ncid, bId,startbathymetry, countbathymetry, bathymetry);
		if(status != NC_NOERR) {handle_error(status);}

		const size_t start[] = {timeDim-1, 0, 0};
		const size_t count[] = {1, yDim, xDim};

		status = nc_get_vara_float(ncid, hId, start, count,water);
		if(status != NC_NOERR){handle_error(status);}

		status = nc_get_vara_float(ncid, huId, start, count, hu);
		if(status != NC_NOERR){handle_error(status);}

		status = nc_get_vara_float(ncid, hvId, start, count, hv);
		if(status != NC_NOERR){handle_error(status);}

		// close file
		status = nc_close(ncid);
		if (status != NC_NOERR) handle_error(status);


		//cout << "yDim: " << yDim << " xDim: " << xDim << endl;
		//cout << "Checkpoint water height:" << endl;
		/*for(int i = 0; i < yDim; i++) {
			for(int j = 0; j < xDim; j++) {
				cout << water[i*xDim + j] << " ";
			}
			cout << endl;
		}*/

	};

	float getWaterHeight(float x, float y){
		int xPos = (int) (((x-getBoundaryPos(BND_LEFT))*xDim)/(getBoundaryPos(BND_RIGHT) - getBoundaryPos(BND_LEFT)));
		int yPos = (int) (((y-getBoundaryPos(BND_BOTTOM))*yDim)/(getBoundaryPos(BND_TOP) - getBoundaryPos(BND_BOTTOM)));
		//cout<<xPos<<" "<<yPos<<" "<<water[yPos*xDim + xPos]<<endl;

		return water[yPos*xDim + xPos];
	};

	float getBathymetry(float x, float y){
		int xPos = (int) (((x-getBoundaryPos(BND_LEFT))*xDim)/(getBoundaryPos(BND_RIGHT) - getBoundaryPos(BND_LEFT)));
		int yPos = (int) (((y-getBoundaryPos(BND_BOTTOM))*yDim)/(getBoundaryPos(BND_TOP) - getBoundaryPos(BND_BOTTOM)));

		return bathymetry[yPos*xDim + xPos];
	};

	float getVeloc_u(float x, float y) {
		int xPos = (int) (((x-getBoundaryPos(BND_LEFT))*xDim)/(getBoundaryPos(BND_RIGHT) - getBoundaryPos(BND_LEFT)));
		int yPos = (int) (((y-getBoundaryPos(BND_BOTTOM))*yDim)/(getBoundaryPos(BND_TOP) - getBoundaryPos(BND_BOTTOM)));

		return hu[yPos*xDim + xPos]/water[yPos*xDim + xPos];
	}

	float getVeloc_v(float x, float y) {
		int xPos = (int) (((x-getBoundaryPos(BND_LEFT))*xDim)/(getBoundaryPos(BND_RIGHT) - getBoundaryPos(BND_LEFT)));
		int yPos = (int) (((y-getBoundaryPos(BND_BOTTOM))*yDim)/(getBoundaryPos(BND_TOP) - getBoundaryPos(BND_BOTTOM)));

		return hv[yPos*xDim + xPos]/water[yPos*xDim + xPos];
	}

	virtual float endSimulation() {
		return l_endSimulation;
	};

	virtual float continueSimulationAt() {
		return l_time;
	}

	virtual int calculatedSteps() {
		return timeDim;
	}

	virtual int getxDim() {
		return xDim;
	}

	virtual int getyDim() {
		return yDim;
	}

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
		if ( i_edge == BND_LEFT )
		  return boundaryPos[2];
		else if ( i_edge == BND_RIGHT)
		  return boundaryPos[3];
		else if ( i_edge == BND_BOTTOM )
		  return boundaryPos[1];
		else
		  return boundaryPos[0];
	};
};


#endif /* SWE_CHECKPOINTSCENARIO_HH_ */
