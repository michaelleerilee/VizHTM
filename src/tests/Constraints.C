/*
 * Constraints.C
 *
 *  Created on: Sep 14, 2019
 *      Author: mrilee
 *
 *  Copyright (C) 2019 Rilee Systems Technologies LLC
 */

#include "STARE.h"
#include "HstmRange.h"
#include "tests.h"

#include "SpatialConstraint.h"
#include "SpatialRange.h"
#include "SpatialInterface.h"

#include <mfhdf.h>
#include <hdf.h>
#include <HdfEosDef.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

void plotBlockingSphere(VizHTM* viz, float r, float g, float b, float radius);
void testShapeFiles(VizHTM *viz, float r, float g, float b, float deltaZ);
void testTenDegreeGrid(VizHTM *viz,float r0, float g0, float b0, float rgbScale);


SpatialRange* addConstraint( VizHTM *viz, SpatialVector a, float64 aa, int resolution, float64 r, float64 g, float64 b, float64 z) {
		SpatialConstraint c0(a,aa);
		RangeConvex rc; rc.add(c0);
		// void intersect(const SpatialIndex * index, HtmRange *hr, bool varlen, HtmRange *hrInterior = 0, HtmRange *hrBoundary = 0);

		// int resolution = 8; // really expensive if varlen is false!
		STARE index;
		SpatialIndex sIndex = index.getIndex(resolution);
		HtmRange hr, hrInterior, hrBoundary;

		// varlen false - rc.intersect(&sIndex,&hr,false,&hrInterior,&hrBoundary); // all's at one resolution
		rc.intersect(&sIndex,&hr,true,&hrInterior,&hrBoundary); // variable resolution coolness

		hr.reset();
		KeyPair kp;
		STARE_SpatialIntervals intervals;
		while( hr.getNext(kp) ) {
			STARE_SpatialIntervals interval = spatialIntervalFromHtmIDKeyPair(kp);
			// cout << setw(16) << setfill(' ') << hex << kp.lo << " " << kp.hi << " : " << interval[0] << " " << interval[1] << endl << flush;
			intervals.push_back(interval[0]);
			if(kp.lo != kp.hi) {
				intervals.push_back(interval[1]);
			}
		}
		hr.reset();

		SpatialRange *sRange = new SpatialRange(intervals);
		viz->addSpatialRange(sRange,r,g,b,0,1,true,z);

		return sRange;
}

SpatialRange* addRangeConvex( VizHTM *viz, RangeConvex rc, int resolution, float64 r, float64 g, float64 b, float64 z) {
		// SpatialConstraint c0(a,aa);
		// RangeConvex rc; rc.add(c0);
		// void intersect(const SpatialIndex * index, HtmRange *hr, bool varlen, HtmRange *hrInterior = 0, HtmRange *hrBoundary = 0);

		// int resolution = 8; // really expensive if varlen is false!
		STARE index;
		SpatialIndex sIndex = index.getIndex(resolution);
		HtmRange hr, hrInterior, hrBoundary;

		// varlen false - rc.intersect(&sIndex,&hr,false,&hrInterior,&hrBoundary); // all's at one resolution
		rc.intersect(&sIndex,&hr,true,&hrInterior,&hrBoundary); // variable resolution coolness

		hr.reset();
		KeyPair kp;
		STARE_SpatialIntervals intervals;
		while( hr.getNext(kp) ) {
			STARE_SpatialIntervals interval = spatialIntervalFromHtmIDKeyPair(kp);
			// cout << setw(16) << setfill(' ') << hex << kp.lo << " " << kp.hi << " : " << interval[0] << " " << interval[1] << endl << flush;
			intervals.push_back(interval[0]);
			if(kp.lo != kp.hi) {
				intervals.push_back(interval[1]);
			}
		}
		hr.reset();

		SpatialRange *sRange = new SpatialRange(intervals);
		viz->addSpatialRange(sRange,r,g,b,0,1,true,z);

		return sRange;
}

bool Constraints1( VizHTM *viz ) {
	bool ok   = false;
	bool once = false;

	// bool faces_color = false;
	bool faces_color = true;

	// float faces_transparency = 0.05;
	// float faces_transparency = 0.05;
	// float faces_transparency = 0.25;
	// float faces_transparency = 0.35;
	float faces_transparency = 0.5;
	// float faces_transparency = 0.6;

	float64 z_offset = 0.034;
	// float64 face_z_offset = 0.005;
	float64 face_z_offset = 0.031;

	cout << "Constraint1" << endl;
	if( !viz->notestream() ) {
		viz->addNotes(new stringstream);
	}
	*(viz->notestream()) << endl << "Constraints1 Start" << endl;

	if(true) {
		if(viz->getProjection() == "None") {
			plotBlockingSphere(viz,0.2,0.2,0.2,0.98);
			*(viz->notestream()) << endl << "Constraints1 plotBlockingSphere enabled" << endl;
		}
	}

	*(viz->notestream()) << endl << "Constraints1 faces_color " << faces_color << endl;
	*(viz->notestream()) << endl << "Constraints1 faces_transparency " << faces_transparency << endl;

	if( true ){
		float
		r0       = 0.5,
		g0       = 0.5,
		b0       = 0.5,
		rgbScale = 0
		;
		testTenDegreeGrid(viz,r0,g0,b0,rgbScale);
		*(viz->notestream()) << endl << "Constraints1 testTenDegreeGrid enabled" << endl;
	}

	if( true ) {
		testShapeFiles(viz,0.5,1,1,0.03);
		*(viz->notestream()) << endl << "Constraints1 testShapeFiles enabled" << endl;
	}

	if( false ) {
		SpatialVector a;
		SpatialRange *sr0, *sr1, *sr2;

		a.setLatLonDegrees(0, 0);
		sr0 = addConstraint(viz
				,a,0
				,5
				,1,0,0
				,0.001);

		a.setLatLonDegrees(0, 60);
		sr1 = addConstraint(viz
				,a,0
				,5
				,1,1,0
				,0.0012);

		a.setLatLonDegrees(45-90, 30+180);
		sr2 = addConstraint(viz
				,a,0
				,5
				,0,1,1
				,0.0013);

		SpatialRange *sr3 = (*sr0) & (*sr1);
		SpatialRange *sr4 = (*sr3) & (*sr2);

		viz->addSpatialRange(sr4,1,1,1,0,1,true,0.002);

	}

	if( false ) {

		// int resolution = 5;
		int resolution = 8;

		RangeConvex rc;

		float64 lat[6] = {  10, 0,  0, 10, 20, 20 };
		float64 lon[6] = { -10, 0, 20, 30, 20,  0 };
		SpatialVector vs[6], cs[6];

		for( int i = 0; i < 6; ++i ) {
			vs[i].setLatLonDegrees(lat[i], lon[i]);
		}
		for( int i = 0; i < 6; ++i ) {
			cs[i] = vs[i] ^ vs[i < 5 ? i+1 : 0];
			SpatialRange *sr = addConstraint(viz,cs[i],0,resolution
					,1-0.2*i,0.2*i,0.2*i
					,0.0025);
			SpatialConstraint sc(cs[i],0);
			rc.add(sc);
			viz->addArc(vs[i], vs[i < 5 ? i+1 : 0], 1, 1, 1, -1, 1.00275, 0, 16);
		}
		addRangeConvex(viz,rc,resolution,0,1,0,0.003);
	}

	if( false ) {
		RangeConvex rc;

		float64 lat[6] = {  10, 0,  0, 10, 20, 20 };
		float64 lon[6] = { -10, 0, 20, 30, 20,  0 };
		SpatialVector vs[6], cs[6], vs1[6];

		for( int i = 0; i < 6; ++i ) {
			vs[i].setLatLonDegrees(lat[i], lon[i]);
		}

		// int resolution = 4;
		// int resolution = 8;
		for(int resolution=0; resolution <= 10; ++resolution ) {
			// STARE_ArrayIndexSpatialValue ValueFromLatLonDegrees(float64 latDegrees, float64 lonDegrees, int resolutionLevel = STARE_HARDWIRED_RESOLUTION_LEVEL_MAX);
			STARE index;

			RangeConvex rc;
			STARE_ArrayIndexSpatialValue siv[6];
			STARE_SpatialIntervals sIntervals;

			for( int i = 0; i < 6; ++i ) {
				siv[i] = index.ValueFromLatLonDegrees(lat[i], lon[i], resolution);
				sIntervals.push_back(siv[i]);
				vs1[i] = index.SpatialVectorFromValue(siv[i]);
			}

			SpatialRange *sRange = new SpatialRange(sIntervals);

			float64 r=1.0, g=resolution/10.0, b=resolution/10.0, z=0.0035;

			viz->addSpatialRange(sRange,r,g,b,0,1,true,z);

			//		for( int i = 0; i < 6; ++i ) {
			//			cs[i] = vs[i] ^ vs[i < 5 ? i+1 : 0];
			//			SpatialRange *sr = addConstraint(viz,cs[i],0,resolution
			//					,1-0.2*i,0.2*i,0.2*i
			//					,0.0025);
			//			SpatialConstraint sc(cs[i],0);
			//			rc.add(sc);
			//			viz->addArc(vs[i], vs[i < 5 ? i+1 : 0], 1, 1, 1, -1, 1.00275, 0, 16);
			//		}
			//		addRangeConvex(viz,rc,resolution,0,1,0,0.003);
			for( int i = 0; i < 6; ++i ) {
				viz->addArc(vs1[i], vs1[i < 5 ? i+1 : 0], resolution/10.0, 1, 0, -1, 1.003, 0, 16);
			}
		}
		if(false) {
			for( int i = 0; i < 6; ++i ) {
				viz->addArc(vs[i], vs[i < 5 ? i+1 : 0], 1, 1, 1, -1, 1.004, 0, 16);
			}
		}

	}

	if( false ) {

		STARE index;
		RangeConvex rc;

		int resolution = 8;

		float64 lat[6] = {  10, 0,  0, 10, 20, 20 };
		float64 lon[6] = { -10, 0, 20, 30, 20,  0 };
		SpatialVector vs[6], cs[6], vs1[6];
		LatLonDegrees64ValueVector latlon;

		for( int i = 0; i < 6; ++i ) {
			// vs[i].setLatLonDegrees(lat[i], lon[i]);
			// latlon.push_back(LatLonDegrees64(lat[i], lon[i]));
			vs[5-i].setLatLonDegrees(lat[i], lon[i]);
			latlon.push_back(LatLonDegrees64(lat[5-i], lon[5-i]));
		}

		STARE_SpatialIntervals intervals = index.ConvexHull(latlon, resolution);

		SpatialRange cover_hull_sr(intervals);
		viz->addSpatialRange(&cover_hull_sr,0.125,1.0,1.0,0,1,true,0.002);

		// STARE_SpatialIntervals CoverCircleFromLatLonRadiusDegrees(float64 latDegrees, float64 lonDegrees, float64 radius_degrees, int force_resolution_level = -1);
		STARE_SpatialIntervals circle = index.CoverCircleFromLatLonRadiusDegrees(10,35,20,8);
		SpatialRange circle_sr(circle);
		viz->addSpatialRange(&circle_sr,1.0,0.0,0.0,0,1,true,0.003);

		SpatialRange *intersect_sr = cover_hull_sr & circle_sr;
		viz->addSpatialRange(intersect_sr,1.0,1.0,0.0,0,1,true,0.0032);

	}

	if( true ) {

		STARE index;
		RangeConvex rc;

		int resolution = 8;

		// Try a grid...
		int nLat = 101,nLon = 101, nLatLon=nLat*nLon;
		float64 lat[nLat], lon[nLon];
		SpatialVector vs[nLatLon], cs[nLatLon], vs1[nLatLon];
		LatLonDegrees64ValueVector latlon;
		for( int i = 0; i < nLon; ++i ) {
			lon[i] = i*0.5;
		}
		for( int j = 0; j < nLat; ++j ) {
			lat[j] = -40.0 + j*1.0;
		}
		int k = 0;
		for( int i = 0; i < nLon; ++i ) {
			for( int j = 0; j < nLat; ++j ) {
				vs[k++].setLatLonDegrees(lat[j], lon[i]);
				latlon.push_back(LatLonDegrees64(lat[j], lon[i]));
				viz->addSphere(vs[k-1], 0, 1, 0, 0.001);
			}
		}

		STARE_SpatialIntervals intervals = index.ConvexHull(latlon, resolution);

		SpatialRange cover_hull_sr(intervals);
		viz->addSpatialRange(&cover_hull_sr,0.125,1.0,1.0,0,1,true,0.002);

		// STARE_SpatialIntervals CoverCircleFromLatLonRadiusDegrees(float64 latDegrees, float64 lonDegrees, float64 radius_degrees, int force_resolution_level = -1);
		STARE_SpatialIntervals circle = index.CoverCircleFromLatLonRadiusDegrees(40,35,20,8);
		SpatialRange circle_sr(circle);
		viz->addSpatialRange(&circle_sr,1.0,0.0,0.0,0,1,true,0.003);

		// SpatialRange *intersect_sr = cover_hull_sr & circle_sr;
		SpatialRange *intersect_sr = circle_sr & cover_hull_sr;
		// SpatialRange *intersect_sr = sr_intersect(circle_sr,cover_hull_sr,true);
		// viz->addSpatialRange(intersect_sr,0.5,0.5,0.0,0,1,true,0.0035);

		intersect_sr->range->range->CompressionPass();
		viz->addSpatialRange(intersect_sr,1.0,1.0,0.0,0,1,true,0.005);

	}

	cout << "Done...";
	ok = true;
	return ok;
}

