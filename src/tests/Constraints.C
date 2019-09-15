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
		viz->addSpatialRange(sRange->range,r,g,b,0,1,true,z);

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
	*(viz->notestream()) << endl << "Granule1 Start" << endl;

	if(true) {
		if(viz->getProjection() == "None") {
			plotBlockingSphere(viz,0.2,0.2,0.2,0.98);
			*(viz->notestream()) << endl << "Granule1 plotBlockingSphere enabled" << endl;
		}
	}

	*(viz->notestream()) << endl << "Granule1 faces_color " << faces_color << endl;
	*(viz->notestream()) << endl << "Granule1 faces_transparency " << faces_transparency << endl;

	if( true ){
		float
		r0       = 0.5,
		g0       = 0.5,
		b0       = 0.5,
		rgbScale = 0
		;
		testTenDegreeGrid(viz,r0,g0,b0,rgbScale);
		*(viz->notestream()) << endl << "Granule1 testTenDegreeGrid enabled" << endl;
	}

	if( true ) {
		testShapeFiles(viz,0.5,1,1,0.03);
		*(viz->notestream()) << endl << "Granule1 testShapeFiles enabled" << endl;
	}

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

	cout << 100 << endl << flush;
	SpatialRange *sr3 = (*sr0) & (*sr1);
	SpatialRange *sr4 = (*sr3) & (*sr2);
	cout << 200 << endl << flush;
	viz->addSpatialRange(sr4->range,1,1,1,0,1,true,0.002);

	if(false) {
		cout << dec << 100 << endl << flush;

		SpatialConstraint c0(a,0);
		RangeConvex rc; rc.add(c0);
		// void intersect(const SpatialIndex * index, HtmRange *hr, bool varlen, HtmRange *hrInterior = 0, HtmRange *hrBoundary = 0);

		int resolution = 8; // really expensive if varlen is false!
		STARE index;
		SpatialIndex sIndex = index.getIndex(8);
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

		SpatialRange sRange(intervals);
		viz->addSpatialRange(sRange.range,1,1,1,0,1,true,0.001);
	}

	if(false) {
		cout << dec << 200 << endl << flush;
		SpatialVector a; a.setLatLonDegrees(45, 45);
		SpatialConstraint c0(a,0);
		RangeConvex rc; rc.add(c0);
		// void intersect(const SpatialIndex * index, HtmRange *hr, bool varlen, HtmRange *hrInterior = 0, HtmRange *hrBoundary = 0);

		int resolution = 8; // really expensive if varlen is false!
		STARE index;
		SpatialIndex sIndex = index.getIndex(8);
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

		SpatialRange sRange(intervals);
		viz->addSpatialRange(sRange.range,1,0,1,0,1,true,0.002);
	}
	try {
	if(false) {
		cout << dec << 300 << endl << flush;

		SpatialVector a; a.setLatLonDegrees(0, 45);
		SpatialConstraint c0(a,0);
		RangeConvex rc; rc.add(c0);
		// void intersect(const SpatialIndex * index, HtmRange *hr, bool varlen, HtmRange *hrInterior = 0, HtmRange *hrBoundary = 0);

		int resolution = 8; // really expensive if varlen is false!
		STARE index;
		SpatialIndex sIndex = index.getIndex(8);
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

		SpatialRange sRange(intervals);
		viz->addSpatialRange(sRange.range,0,0,1,0,1,true,0.0025);
	}
	} catch ( SpatialException e ) {
		cout << "e " << e.what() << endl << flush;
	}

	cout << "Done...";
	ok = true;
	return ok;
}

