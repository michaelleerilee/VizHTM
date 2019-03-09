/*
 * PolePosition.C
 *
 *  Created on: Feb 27, 2019
 *      Author: mrilee
 *
 *  Copyright (C) 2019 Rilee Systems Technologies LLC
 */

#include "tests.h"

#include <SpatialRotation.h>


SpatialVector nudge(SpatialVector v0, SpatialVector v1, float64 a1) {
	SpatialVector u = (v0*(1-a1)) + (v1*a1); // u.normalize();
	return u;
}

bool PolePosition1(VizHTM *viz) {
	bool ok = false;

	// float64 sphere_scale = 0.025;
	float64 sphere_scale = 0.05;
	// float64 corner_nudge = 0.25;
	float64 corner_nudge = 0.05;
	// float64 corner_nudge = 0.0;

	STARE index;
	// STARE index1;

	if( true ){
		SpatialVector axis     = 0.5*xhat + 0.5*yhat; axis.normalize();
		float64       theta    = 0.25*gPi;
		// theta = 0.0;
		// theta    = 0.125*gPi;
		// theta    = 0.5*gPi + 0.00001;
		theta    = 0.5*gPi;


		SpatialRotation rotate_root_octahedron = SpatialRotation(axis,theta);
		int search_level = 27, build_level = 5;

		index = STARE(search_level, build_level, rotate_root_octahedron);
	}

	int level = 27; // too deep
	// int level = 25; // too deep
	//+ int level = 24; // Floating point gets really bad on interactive visualization
	// int level = 23; // Interactive viz not great
	// int level = 22; // Interactive viz still rough, but better
	// int level = 21; // Interactive viz a little rough, but much better
	// int level = 20; // Interactive viz a little rough, but much better
	// int level = 19; // ~ 20m Interactive viz pretty good
	// int level = 10;
	// int level = 6;

	float64 lat0 = 90, lon0 = 0.0;
	STARE_ArrayIndexSpatialValue north_pole_sid = index.ValueFromLatLonDegrees(lat0,lon0,level);
	// STARE_ArrayIndexSpatialValue north_pole_sid = index.ValueFromLatLonDegrees(15.0,15.0,level);

	SpatialVector graphicsOrigin(0.0,0.0,0.0); // graphicsOrigin.setLatLonDegrees(lat0, lon0);
	float64       graphicsScale = 1.0;

	// float64 scale_delta = 1.0e-6, scale = 1.0;
	float64 scale_delta = 0.0, scale = 1.0;
	float64 a = 0;

	cout << "level: " << index.ResolutionLevelFromValue(north_pole_sid)
			<< ", scale(m): "
			<< setprecision(17)
			<< setw(20)
			<< scientific
			<< index.lengthMeterScaleFromEdgeFromLevel(level)
			<< endl << flush;
	cout << "north_pole_sid: " << hex << "0x" << north_pole_sid << dec << endl << flush;
	STARE_ArrayIndexSpatialValues neighbors = index.NeighborsOfValue(north_pole_sid);

	if(false) {
		Triangle tr0 = index.TriangleFromValue(north_pole_sid,level);
		graphicsOrigin = tr0.centroid;
		graphicsScale  = 1.0e+8;
	}
	cout << "graphicsOrigin: " << graphicsOrigin << endl << flush;

	level = index.ResolutionLevelFromValue(north_pole_sid);
	Triangle ta[12];
	for(int i=0; i<12; ++i) {
		ta[i] = index.TriangleFromValue(neighbors[i],level);
	}
	cout << ".." << endl << flush;
	for(int i=0; i<12; ++i) {
		for(int j=0; j <= i; ++j ) {
			cout << i << "," << j << " delta = "
					<< (ta[i].centroid-ta[j].centroid) << ", length = "
					<< (ta[i].centroid-ta[j].centroid).length()
					<< endl << flush;
		}
	}

	// for(int i = 8; i >= 0; --i ) {
	// for(int i = 0; i < 9; ++i ) {
	// for(int i = 4; i < 8; ++i ) {
	for(int i = 0; i < 12; ++i ) {
		Triangle tr0 = ta[i];
		if( true ) {
				for( int k=0; k < 3; ++k ) {
					tr0.vertices[k] = (tr0.vertices[k]-graphicsOrigin) * graphicsScale;
				}
				tr0.centroid = (tr0.centroid-graphicsOrigin) * graphicsScale;
		}
		cout << "-- neighbor = " << i << " --" << endl << flush;
		for(int j = 0; j < 3; ++j ) {
			cout << j
					<< " tr0.v: "
					<< setprecision(17)
					<< setw(20)
					<< scientific
					<< tr0.vertices[j] << endl << flush;
		}
		cout << "c" << " tr0.c: " << tr0.centroid << endl << flush;

		// SpatialVector *v = new SpatialVector(tr0.centroid*1.0);
		SpatialVector *v = new SpatialVector(tr0.centroid);

		/*
		char *str = new char[256]; char tmpName[64];
		index->nameById(numericId,tmpName);
		sprintf(str,"%s\nid: %llx\nix: %llu\n",tmpName,numericId,nodeIndex);
		viz->addAnnotation((new SpatialVector(x)),str,size,r,g,b);
		*/

		char *str = new char[256];
		sprintf(str,"%s\n",to_string(i).c_str());
		// viz->makeText(a,str,0.25*0.5*3.14*pow(2.0,-level),0.0,1.0,1.0);
		// viz->addAnnotation(a,"test",0.25*0.5*3.14*pow(2.0,-level),0.0,1.0,1.0);
		viz->addAnnotation(v,str,0.25*0.5*3.14*pow(2.0,-level),0.0,1.0,1.0);

		if(true) {
			// SpatialVector *v = new SpatialVector(tr0.centroid*(scale+scale_delta));
			// viz->makeText(v,to_string(i).c_str(),0.25*0.5*3.14*pow(2.0,-level),0.0,1.0,1.0);

			if(true) {
				viz->addEdge(
						nudge(tr0.vertices[0],tr0.centroid,corner_nudge),
						nudge(tr0.vertices[1],tr0.centroid,corner_nudge),
						1.0, 0.0, 0.0, a, scale);
			}
			if(true) {
				viz->addEdge(
						nudge(tr0.vertices[1],tr0.centroid,corner_nudge),
						nudge(tr0.vertices[2],tr0.centroid,corner_nudge),
						0.0, 1.0, 0.0, a, scale);
			}
			if(true) {
				viz->addEdge(
						nudge(tr0.vertices[2],tr0.centroid,corner_nudge),
						nudge(tr0.vertices[0],tr0.centroid,corner_nudge),
						0.0, 0.0, 1.0, a, scale);
			}

			if(false) {
				viz->addSphere(nudge(tr0.vertices[0],tr0.centroid,corner_nudge)
						,1.0,0.0,0.0,0.125*0.5*3.14*pow(2.0,-level)*graphicsScale*sphere_scale);
				viz->addSphere(nudge(tr0.vertices[1],tr0.centroid,corner_nudge)
						,0.0,1.0,0.0,0.125*0.5*3.14*pow(2.0,-level)*graphicsScale*sphere_scale);
				viz->addSphere(nudge(tr0.vertices[2],tr0.centroid,corner_nudge)
						,0.0,0.0,1.0,0.125*0.5*3.14*pow(2.0,-level)*graphicsScale*sphere_scale);

				// viz->addSphere(tr0.centroid,1.0,1.0,1.0,0.25*0.5*3.14*pow(2.0,-level)*graphicsScale*sphere_scale);
			}

		}
		scale += scale_delta;
	}



	if(true){
		for(int resolution_level = level; resolution_level < level+1; ++resolution_level ) {

			Triangle tr0 = index.TriangleFromValue(north_pole_sid,resolution_level);

			viz->addSphere(tr0.centroid
					,0.5,0.5,1.0,0.125*0.5*3.14*pow(2.0,-level)*graphicsScale*sphere_scale);

			float
			color_scale = 1.0,
			r0 = 1, g0 = 0, b0 = 0,
			r1 = 0, g1 = 1, b1 = 0,
			r2 = 0, g2 = 0, b2 = 1,
			a0 = 0, a1 = 0, a2 = 0;

			cout << ".." << endl << flush;
			for(int j = 0; j < 3; ++j ) {
				cout << j
						<< " tr0.v: "
						<< setprecision(17)
						<< setw(20)
						<< scientific
						<< tr0.vertices[j] << endl << flush;
			}
			cout << "c" << " tr0.c: " << tr0.centroid << endl << flush;

			// viz->addFace3(tr0.vertices[0], tr0.vertices[1], tr0.vertices[2], r0, g0, b0, r1, g1, b1, r2, g2, b2, a0, a1, a2, scale);

			if(true) {
				viz->addEdge(tr0.vertices[0], tr0.vertices[1], 0.0, 1.0, 0.0, a, scale);
				viz->addEdge(tr0.vertices[1], tr0.vertices[2], 0.0, 1.0, 0.0, a, scale);
				viz->addEdge(tr0.vertices[2], tr0.vertices[0], 0.0, 1.0, 0.0, a, scale);
			}

			if(false) {
				viz->addEdge(tr0.vertices[0], tr0.centroid, 0.0, 1.0, 1.0, a, scale);
				viz->addEdge(tr0.vertices[1], tr0.centroid, 0.0, 1.0, 1.0, a, scale);
				viz->addEdge(tr0.vertices[2], tr0.centroid, 0.0, 1.0, 1.0, a, scale);
			}

			scale += scale_delta;
		}
	}

	if(false){
		SpatialVector axis     = 0.5*xhat + 0.5*yhat; axis.normalize();
		float64       theta    = 0.25*gPi;
		theta = 0.0;

		SpatialRotation rotate_root_octahedron = SpatialRotation(axis,theta);
		// int search_level = 27, build_level = 5;
		// SpatialIndex sIndex                    = SpatialIndex(search_level, build_level, rotate_root_octahedron);

		SpatialVector new_pole = rotate_root_octahedron.rotated_from(zhat);

		viz->addSphere(new_pole
				,1.0,1.0,0.0,0.125*0.5*3.14*pow(2.0,-level)*graphicsScale*sphere_scale);
	}


	if(true){
		SpatialIndex sIndex = index.getIndex(level);

		uint64 np_htmid = index.htmIDFromValue(north_pole_sid, level);
		SpatialVector workspace_ev[18];
		uint64 neighbors_e[3];
		sIndex.NeighborsAcrossEdgesFromHtmId(neighbors_e, np_htmid, workspace_ev);

		uint64 neighbors_v[9];
		sIndex.NeighborsAcrossVerticesFromEdges(neighbors_v, neighbors_e, np_htmid, workspace_ev);

		cout << "  np_htmid: 0x" << hex << np_htmid << dec << endl << flush;
		for(int i=0; i<3; ++i) {
			cout << i << " n edge id 0x" << hex << neighbors_e[i] << dec << endl << flush;
		}

		uint64 neighbors_[12];
		for(int i=0; i<9; ++i) {
			neighbors_[i] = neighbors_v[i];
		}
		for(int i=9; i<12; ++i) {
			neighbors_[i] = neighbors_e[i-9];
		}

		cout << 90 << endl << flush;

		for(int i=0; i<3; ++i) {
			cout << i << " v "
					<< setprecision(17)
					<< setw(20)
					<< scientific
					<< workspace_ev[i] << endl << flush;
		}

		cout << endl << flush;
		for(int i=3; i<6; ++i) {
			cout << i << " m "
					<< setprecision(17)
					<< setw(20)
					<< scientific
					<< workspace_ev[i] << endl << flush;
		}

		cout << endl << flush;

		cout << 0 << " q "
				<< setprecision(17)
				<< setw(20)
				<< scientific
				<< workspace_ev[9+3] << " -- " << workspace_ev[8] << endl << flush;

		/*
		viz->addSphere(workspace_ev[1],0.75,1.0,1.0,0.125*0.5*3.14*pow(2.0,-level)*graphicsScale*sphere_scale*2);
		viz->addSphere(workspace_ev[3],0.25,1.0,1.0,0.125*0.5*3.14*pow(2.0,-level)*graphicsScale*sphere_scale*2);
		 */

		if(false) {
			// The 4th vertex neighbor guess
			viz->addSphere(workspace_ev[9+3],1.0,0.75,0.75,0.125*0.5*3.14*pow(2.0,-level)*graphicsScale*sphere_scale*2);
			viz->addEdge(workspace_ev[9+3],workspace_ev[9+3]*1.025,1.0,0.75,0.75,a,1.0);

			// The 3rd edge neighbor guess: v is wev[0:2], m is wev[3:5], edge is wev[6:8], vert is wev[9:17]
			viz->addSphere(workspace_ev[8],1.0,0.25,0.25,0.125*0.5*3.14*pow(2.0,-level)*graphicsScale*sphere_scale*2);
			viz->addEdge(workspace_ev[8],workspace_ev[8]*1.025,1.0,0.25,0.25,a,1.0);
		}

		if(true) {
			viz->addEdge(workspace_ev[0], workspace_ev[1], 0.0, 1.0, 1.0, a, 1.0);
			viz->addEdge(workspace_ev[1], workspace_ev[2], 0.0, 1.0, 1.0, a, 1.0);
			viz->addEdge(workspace_ev[2], workspace_ev[0], 0.0, 1.0, 1.0, a, 1.0);
		}

		for(int i=0; i<1; ++i) {
			viz->addSphere(workspace_ev[6+i],1.0,0.75,0.25,0.125*0.5*3.14*pow(2.0,-level)*graphicsScale*sphere_scale);
			cout << i << " edge sIdx.idByPoint 0x" << hex << sIndex.idByPoint(workspace_ev[6+i]) << dec << endl << flush;
		}

	}

	cout << "PolePosition done" << endl << flush;
	ok = true;
	return ok;
}
