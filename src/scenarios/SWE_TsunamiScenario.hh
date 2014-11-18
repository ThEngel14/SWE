/*
 * SWEArtificialTsunamiScenario.h
 *
 *  Created on: Nov 13, 2014
 *      Author: kyu
 */

#ifndef SWETSUNAMISCENARIO_HH_
#define SWETSUNAMISCENARIO_HH_

#include "SWE_Scenario.hh"
#include <cmath>
#include <limits>
#include <netcdf.h>
//See https://www.unidata.ucar.edu/software/netcdf/docs/netcdf-c.pdf for netcdf instructions

#define BATHFILE "NCScenario/artificialtsunami_bathymetry_1000.nc"
#define DISPFILE "NCScenario/artificialtsunami_displ_1000.nc"


class SWE_TsunamiScenario: public SWE_Scenario {
private:
	static const float scale = 1/100.0f; //use the denominator as cell number
	struct Position{
		int x;
		int y;
	};

	enum Mode{
		BATHYMETRY, DISPLACEMENT
	};

	size_t xBathSize;
	size_t yBathSize;

	size_t xDisplSize;
	size_t yDisplSize;

	float *xBathVals;
	float *yBathVals;
	float *zBathVals;

	float *xDisplVals;
	float *yDisplVals;
	float *zDisplVals;

	/**
	 * Get the index of the closest position in x/y####Vals (e.g. xBathVals, yDisplVals etc.) depending on
	 * the mode parameter.
	 * @param x x-coordinate
	 * @param y y-coordinate
	 * @param m BATHYMETRY for x/yBathVals or DISPLACEMENT for x/yDisplVals.
	 * @return (index(x####Vals) , index(y####Vals))-tuple which is closest to the parameter x and y
	 */
	Position getClosestPosition(float x, float y, Mode m){
		Position ret;
		float dif;
		int i;
		if(m == BATHYMETRY){
			i=0;
			do{
				dif = std::abs(xBathVals[i] - x);
				i++;
			}while(i<xBathSize && std::abs(xBathVals[i] - x) < dif);
			ret.x = i-1;

			i=0;
			do{
				dif = std::abs(yBathVals[i] - y);
				i++;
			}while(i<yBathSize && std::abs(yBathVals[i] - y) < dif);
			ret.y = i-1;
		}else{
			i=0;
			do{
				dif = std::abs(xDisplVals[i] - x);
				i++;
			}while(i<xDisplSize && std::abs(xDisplVals[i] - x) < dif);
			ret.x = i-1;

			i=0;
			do{
				dif = std::abs(yDisplVals[i] - y);
				i++;
			}while(i<yDisplSize && std::abs(yDisplVals[i] - y) < dif);
			ret.y = i-1;
		}
		return ret;
	};

	/**
	 * Open a file and determine the dimension on x and y axis.
	 * @param file path of the file
	 * @param ncid output variable to store the ID of the file
	 * @param xDim dimension on x-axis
	 * @param yDim dimension on y-axis
	 */
	void openNetcdf(const char *file, int *ncid, size_t *xDim, size_t *yDim){
		int status;
		status = nc_open(file, 0, ncid);
		if(status != NC_NOERR){handle_error(status);}

		status = nc_inq_dimlen(*ncid, 0, xDim);
		if (status != NC_NOERR) handle_error(status);

		status = nc_inq_dimlen(*ncid, 1, yDim);
		if (status != NC_NOERR) handle_error(status);
	}

	/**
	 * Reads data from a given nc-file containing variables x,y and z and saves into output arrays.
	 * Before calling this function, be sure to open the file with openNetcdf(...).
	 * @param ncid ID of the file
	 * @param xsize Dimension of the x-axis
	 * @param ysize Dimension of the y-axis
	 * @param xVals output array for all x-values
	 * @param yVals output array for all y-values
	 * @param zVals output array for all z-values
	 */
	void readcloseNetcdf(int ncid, size_t xsize, size_t ysize, float *xVals, float *yVals, float *zVals){
		int status;
		int zid;
		const size_t start[] = {0, 0}; 			// start at first value

		const size_t count[] = {xsize, ysize};
		status = nc_inq_varid(ncid, "x", &zid);
		if(status != NC_NOERR){handle_error(status);}
		status = nc_get_vara_float(ncid, zid, start, count, xVals);
		if(status != NC_NOERR){handle_error(status);}

		status = nc_inq_varid(ncid, "y", &zid);
		if(status != NC_NOERR){handle_error(status);}
		status = nc_get_vara_float(ncid, zid, start, count, yVals);
		if(status != NC_NOERR){handle_error(status);}

		status = nc_inq_varid(ncid, "z", &zid);
		if(status != NC_NOERR){handle_error(status);}
		status = nc_get_vara_float(ncid, zid, start, count, zVals);
		if(status != NC_NOERR){handle_error(status);}

		status = nc_close(ncid);
		if (status != NC_NOERR) handle_error(status);


	};

	/**
	 * Print error message.
	 */
	void handle_error(int status) {
		if (status != NC_NOERR) {
			cout<< nc_strerror(status)<<"\n";
			exit(-1);
		}
	};

	/**
	 * Determine the value of the displacement dependend of the given position.
	 * @param x x-coordinate
	 * @param y y-coordinate
	 *
	 */
	float computeDisplacement(float x, float y){
		if(std::fabs(x) <= xDisplVals[xDisplSize-1] && std::fabs(y) <= yDisplVals[yDisplSize-1]){ //true if x and y are in the displ.-square
			Position posDisp = getClosestPosition(x,y,DISPLACEMENT);
			return zDisplVals[posDisp.y * yDisplSize + posDisp.x];
		}
		return 0.0f;
	}

	/**
	 * Scales all values of the arrays with a predefined number to ensure the correct proportion of the scenario.
	 */
	void scaleAllVals(){

		for(int i=0;i < xBathSize;i++){
			xBathVals[i] *= scale;
		}

		for(int i=0;i < yBathSize;i++){
			yBathVals[i] *= scale;
		}

		int zBathSize = xBathSize * yBathSize;
		for(int i=0;i < zBathSize;i++){
			zBathVals[i] *= scale;
		}

		for(int i=0;i < xDisplSize;i++){
			xDisplVals[i] *= scale;
		}
		for(int i=0;i < yDisplSize;i++){
			yDisplVals[i] *= scale;
		}

		int zDisplSize = xDisplSize * yDisplSize;
		for(int i=0;i < zDisplSize;i++){
			zDisplVals[i] *= scale;
		}


	}
public:
	/**
	 * Open and preparing a set of nc-files for bathymetry and displacement data.
	 */
	SWE_TsunamiScenario(){
		int ncBathid;
		openNetcdf(BATHFILE,&ncBathid, &xBathSize, &yBathSize);
		xBathVals = (float*)malloc(xBathSize * sizeof(float));
		yBathVals = (float*)malloc(yBathSize * sizeof(float));
		zBathVals = (float*)malloc(xBathSize * yBathSize * sizeof(float));
		readcloseNetcdf(ncBathid,xBathSize, yBathSize, xBathVals, yBathVals, zBathVals);

		int ncDisplid;
		openNetcdf(DISPFILE,&ncDisplid, &xDisplSize, &yDisplSize);
		xDisplVals = (float*)malloc(xDisplSize * sizeof(float));
		yDisplVals = (float*)malloc(yDisplSize * sizeof(float));
		zDisplVals = (float*)malloc(xDisplSize * yDisplSize * sizeof(float));
		readcloseNetcdf(ncDisplid, xDisplSize, yDisplSize, xDisplVals, yDisplVals, zDisplVals);

		scaleAllVals();

	};

	float getWaterHeight(float x, float y){
		Position pos = getClosestPosition(x,y,BATHYMETRY);
		return -std::min(zBathVals[pos.y * yBathSize + pos.x], 0.0f) + computeDisplacement(x,y);
	};

	float getBathymetry(float x, float y){
		Position pos = getClosestPosition(x,y,BATHYMETRY);
		float bath = zBathVals[pos.y * yBathSize + pos.x] + computeDisplacement(x,y);
		if(bath >= -20 * scale && bath <= 20 * scale){
			if(bath > 0){
				return 20 * scale;
			}else{
				return -20 * scale;
			}
		}
		return bath;
	};

	virtual float endSimulation() { return (float) 15; };

	/**
	 *
	 * @param edge BoundaryEdge to request the BoundaryType
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
		int xboundary = (xBathSize/2) * std::floor( xBathSize / (xBathVals[xBathSize-1]*2) + 0.5f);
										//compute conversionrate from dimension to meter
		int yboundary = (yBathSize/2) * std::floor( yBathSize / (yBathVals[yBathSize-1]*2) + 0.5f);
	    if ( i_edge == BND_LEFT )
	      return (float)-xboundary * scale;
	    else if ( i_edge == BND_RIGHT)
	      return (float)xboundary * scale;
	    else if ( i_edge == BND_BOTTOM )
	      return (float)-yboundary * scale;
	    else
	      return (float)yboundary * scale;
	};
};

#endif /* SWEARTIFICIALTSUNAMISCENARIO_HH_ */
