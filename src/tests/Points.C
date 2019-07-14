/*
 * Points.C
 *
 *  Created on: Jun 26, 2019
 *      Author: mrilee
 *
 *  Copyright (C) 2019 Rilee Systems Technologies LLC
 */

#include "STARE.h"
#include "HstmRange.h"
#include "tests.h"

#include <mfhdf.h>
#include <hdf.h>
#include <HdfEosDef.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
using namespace std;

void plotBlockingSphere(VizHTM* viz, float r, float g, float b, float radius);
void testShapeFiles(VizHTM *viz, float r, float g, float b, float deltaZ);
void testTenDegreeGrid(VizHTM *viz,float r0, float g0, float b0, float rgbScale);

//void pointField(VizHTM *viz,
//		float lat0, float lon0, float dlat, float dlon, int nlat, int nlon
//		) {
//
//}

//class PointField {
//
//};



bool Points1( VizHTM *viz ) {
	bool ok   = false;
	plotBlockingSphere(viz,0.2,0.2,0.2,0.999);
	testTenDegreeGrid(viz,0.5,0.5,0.5,0);
	testShapeFiles(viz,0.5,1,1,0.0);

	float64
	lat = 60,
	lon = 123.4;

	// int level = 20;
	// int level = 10;
	// int level = 8;
	int level = 4;

	STARE index;
	HstmRange* range = new HstmRange;
	float64 z_offset = 0.03;

	/*
	cout << 200 << endl << flush;
	cout
	<< setprecision(17)
	<< setw(20)
	<< scientific;
	*/

	range->purge();

	double g = 1.0;
	int
	i0 = 0,
	i1 = 7;
	double delta = 0.0001;
	for(int i=i0; i<i1; ++i) {
		z_offset += delta;
		level = i;
		SpatialIndex sIndex = index.getIndex(level);
		STARE_ArrayIndexSpatialValue a = index.ValueFromLatLonDegrees(lat,lon,level);
		EmbeddedLevelNameEncoding leftJustified;
		leftJustified.setIdFromSciDBLeftJustifiedFormat(a);
		uint64 b = leftJustified.getId();
		// range->reset();
		range->addRange(b,b);
		g = 0.5*(1+((1+i-i0)/(i1-i0)));
		viz->addHstmRange(range,0.0,g,0.0,0.0,1.0,true,z_offset,&sIndex);
		range->purge();
	}

	SpatialVector x; x.setLatLonDegrees(lat, lon);
	viz->addSphere(x, 1.0, 0.0, 1.0, 0.005);

	ok = true;

	return ok;
}

void addSIndex(
		VizHTM *viz,
		STARE_ArrayIndexSpatialValue idx,
		float64 r, float64 g, float64 b, float64 a, float64 scale, bool arcFlag, float64 z_offset) {
	STARE index;
	HstmRange* range = new HstmRange;
	int level = index.ResolutionLevelFromValue(idx);
	SpatialIndex sIndex = index.getIndex(level);
	EmbeddedLevelNameEncoding leftJustified;
	leftJustified.setIdFromSciDBLeftJustifiedFormat(idx);
	uint64 idx_hstm = leftJustified.getId();
	range->addRange(idx_hstm, idx_hstm);
	viz->addHstmRange(range,r,g,b,a,scale,arcFlag,z_offset,&sIndex);
	delete range;
}

void cmpSIndices(
		VizHTM *viz,
		STARE_ArrayIndexSpatialValue idx0,
		STARE_ArrayIndexSpatialValue idx1,
		float64 scale,
		float64 z_offset
		) {

	STARE_ArrayIndexSpatialValue idx[2];
	idx[0] = idx0; idx[1] = idx1;
	double r[2],g[2],b[2],a = 0;
	for( int i = 0; i < 2; ++i ) {
		r[i] = 1-i; g[i] = 0; b[i] = i;
	}
	bool arcFlag = true;

	int overlap = cmpSpatial(idx0,idx1);
//	cout << " idx 0,1: "
//			<< setw(16) << setfill('0') << hex << idx0 << " "
//			<< setw(16) << setfill('0') << hex << idx1 << " "
//			<< ", overlap: " << overlap
//			<< endl << flush;

	if( overlap == -1 ) {
		// idx0 contains idx1
		for( int i = 0; i < 2; ++i ) {
			r[i] = 0.75; g[i] = 1.0; b[i] = 0.0;
		}
	} else if( overlap == 1 ) {
		// idx1 contains idx0
		for( int i = 0; i < 2; ++i ) {
			r[i] = 0.0; g[i] = 1.0; b[i] = 0.75;
		}
	}
//	else {
//		r[0] = 1.0;
//		r[1] = 1.0;
//	}

	for( int i = 0; i < 2; ++i ) {
		addSIndex(viz,idx[i],r[i],g[i],b[i],a,scale,arcFlag,z_offset);
	}
}


bool Points2( VizHTM *viz ) {
	bool ok   = false;
	plotBlockingSphere(viz,0.2,0.2,0.2,0.999);
	testTenDegreeGrid(viz,0.5,0.5,0.5,0);
	testShapeFiles(viz,0.5,1,1,0.0);

	STARE index;

	float64
	lat = 60,
	lon = 123.4;

	int level = 4;

	double
	r = 0.0,
	g = 1.0,
	b = 0.0,
	a = 0.0;

	bool
	arcFlag = true;
	int
	i0 = 0,
	i1 = 7;
	double
	z_offset = 0.0001,
	delta    = 0.0001,
	scale    = 1.0;

	for(int i=i0; i<i1; ++i) {
		z_offset += delta;
		level = i;
		STARE_ArrayIndexSpatialValue idx = index.ValueFromLatLonDegrees(lat,lon,level);
		addSIndex(
				viz,
				idx,
				r, g, b, a, scale, arcFlag, z_offset
				);
	}

	SpatialVector x; x.setLatLonDegrees(lat, lon);
	viz->addSphere(x, 1.0, 0.0, 1.0, 0.005);

	ok = true;
	return ok;
}

bool Points3( VizHTM *viz ) {
	bool ok   = false;
	plotBlockingSphere(viz,0.2,0.2,0.2,0.999);
	testTenDegreeGrid(viz,0.5,0.5,0.5,0);
	testShapeFiles(viz,0.5,1,1,0.0);

	STARE index;
	SpatialVector x;

	// float64 radius = 0.00025;
	float64 radius = 0.0025;

	float64 lat, lon;
	int level;

	STARE_ArrayIndexSpatialValue idx0, idx1;

	if(true) {

	lat =  38.0387;	lon = -76.3221; level = 6;
	idx0 = index.ValueFromLatLonDegrees(lat,lon,level);
	x.setLatLonDegrees(lat, lon);
	viz->addSphere(x,1.0,1.0,0.0,radius);

	lat = 38; lon = -75; level = 4;
	idx1 = index.ValueFromLatLonDegrees(lat,lon,level);
	x.setLatLonDegrees(lat, lon);
	viz->addSphere(x,1.0,0.0,1.0,radius);

	cmpSIndices(viz,idx0,idx1,1.0,0.0001);


	lat = 40; lon = -75; level = 6;
	idx0 = index.ValueFromLatLonDegrees(lat,lon,level);
	x.setLatLonDegrees(lat, lon);
	viz->addSphere(x,1.0,1.0,0.0,radius);

	lat = 40.5; lon = -75; level = 6;
	idx1 = index.ValueFromLatLonDegrees(lat,lon,level);
	x.setLatLonDegrees(lat, lon);
	viz->addSphere(x,1.0,0.0,1.0,radius);

	cmpSIndices(viz,idx0,idx1,1.0,0.0001);


	lat = 43.5; lon = -75; level = 6;
	idx0 = index.ValueFromLatLonDegrees(lat,lon,level);
	x.setLatLonDegrees(lat, lon);
	viz->addSphere(x,1.0,1.0,0.0,radius);

	lat = 44.0; lon = -74.5; level = 6;
	idx1 = index.ValueFromLatLonDegrees(lat,lon,level);
	x.setLatLonDegrees(lat, lon);
	viz->addSphere(x,1.0,0.0,1.0,radius);

	cmpSIndices(viz,idx0,idx1,1.0,0.0001);


	lat = 47.5; lon = -75.5; level = 6;
	idx0 = index.ValueFromLatLonDegrees(lat,lon,level);
	x.setLatLonDegrees(lat, lon);
	viz->addSphere(x,1.0,1.0,0.0,radius);

	lat = 48.0; lon = -74.5; level = 7;
	idx1 = index.ValueFromLatLonDegrees(lat,lon,level);
	x.setLatLonDegrees(lat, lon);
	viz->addSphere(x,1.0,0.0,1.0,radius);

	cmpSIndices(viz,idx0,idx1,1.0,0.0001);
	}

	// cout << "false" << endl << flush;
	lat = 39.5; lon = -85.5; level = 4;
	idx0 = index.ValueFromLatLonDegrees(lat,lon,level);
	x.setLatLonDegrees(lat, lon);
	viz->addSphere(x,1.0,1.0,0.0,radius);

	lat = 39.5; lon = -84.5; level = 6;
	idx1 = index.ValueFromLatLonDegrees(lat,lon,level);
	x.setLatLonDegrees(lat, lon);
	viz->addSphere(x,1.0,0.0,1.0,radius);

	cmpSIndices(viz,idx0,idx1,1.0,0.0001);

	if(false) {
	cout << "--" << endl << flush;
	lat = 39.5; lon = -84.5; level = 6;
	idx0 = index.ValueFromLatLonDegrees(lat,lon,level);
	x.setLatLonDegrees(lat, lon);
	viz->addSphere(x,1.0,1.0,0.0,radius);

	lat = 39.5; lon = -85.5; level = 4;
	idx1 = index.ValueFromLatLonDegrees(lat,lon,level);
	x.setLatLonDegrees(lat, lon);
	viz->addSphere(x,1.0,0.0,1.0,radius);

	cmpSIndices(viz,idx0,idx1,1.0,0.0001);
	}

	if(false) {
	cout << "true" << endl << flush;
	lat = 39.5; lon = -85.5; level = 6;
	idx0 = index.ValueFromLatLonDegrees(lat,lon,level);
	x.setLatLonDegrees(lat, lon);
	viz->addSphere(x,1.0,1.0,0.0,radius);

	lat = 39.5; lon = -84.5; level = 4;
	idx1 = index.ValueFromLatLonDegrees(lat,lon,level);
	x.setLatLonDegrees(lat, lon);
	viz->addSphere(x,1.0,0.0,1.0,radius);

	cmpSIndices(viz,idx0,idx1,1.0,0.0001);
	cout << "true" << endl << flush;

	cout << "--" << endl << flush;
	lat = 39.5; lon = -84.5; level = 4;
	idx0 = index.ValueFromLatLonDegrees(lat,lon,level);
	x.setLatLonDegrees(lat, lon);
	viz->addSphere(x,1.0,1.0,0.0,radius);

	lat = 39.5; lon = -85.5; level = 6;
	idx1 = index.ValueFromLatLonDegrees(lat,lon,level);
	x.setLatLonDegrees(lat, lon);
	viz->addSphere(x,1.0,0.0,1.0,radius);

	cmpSIndices(viz,idx0,idx1,1.0,0.0001);
	}

	ok = true;
	return ok;
}

bool Points4( VizHTM *viz ) {
	bool ok   = false;
	plotBlockingSphere(viz,0.2,0.2,0.2,0.999);
	testTenDegreeGrid(viz,0.5,0.5,0.5,0);
	testShapeFiles(viz,0.5,1,1,0.0);

	STARE index;
	SpatialVector x;

	// float64 radius = 0.00025;
	float64 radius = 0.0025;

	float64 lat, lon;
	int level;

	STARE_ArrayIndexSpatialValue idx0, idx1;






	ok = true;
	return ok;
}
