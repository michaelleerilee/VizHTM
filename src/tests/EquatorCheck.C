/*
 * EquatorCheck.C
 *
 *  Created on: Mar 10, 2019
 *      Author: mrilee
 *
 *  Copyright (C) 2019 Rilee Systems Technologies LLC
 */

#include "tests.h"

bool EquatorCheck1( VizHTM* viz ) {
	bool ok = false;

	// float64 sphere_scale = 0.025;
	float64 sphere_scale = 0.05;
	// float64 corner_nudge = 0.25;
	float64 corner_nudge = 0.05;
	// float64 corner_nudge = 0.0;

	// Transparency a;
	float64 a = 0.0;

	// Scale the vectors
	float64 scale = 1.0;

	int level = 27; // Not too bad with really small gEpsilon. Floating point gets really bad on interactive visualization

	STARE index;
	// STARE index1;

	if( true ){
		SpatialVector xhat(1,0,0), yhat(0,1,0);
		SpatialVector axis     = 0.5*xhat + 0.5*yhat; axis.normalize();
		float64       theta    = 0.25*gPi;
		// theta = 0.0;
		// theta    = 0.125*gPi;
		// theta    = 0.5*gPi + 0.00001;
		// theta    = 0.5*gPi; // Okay!
		SpatialRotation rotate_root_octahedron = SpatialRotation(axis,theta);
		int search_level = 27, build_level = 5;
		index = STARE(search_level, build_level, rotate_root_octahedron);

		/*
		float64 dtheta = 1.0e-9;
		for(int i=0; i<21; ++i) {
			SpatialRotation rot = SpatialRotation(axis,i*dtheta);
			SpatialVector z1 = rot.rotated_from(zhat);
			if( i == 0 ) {
				viz->addSphere(z1,1.0,i/20.0,0.0,0.125*0.5*3.14*pow(2.0,-level)*sphere_scale*3);
			} else {
				viz->addSphere(z1,1.0,i/20.0,0.0,0.125*0.5*3.14*pow(2.0,-level)*sphere_scale);
			}
			if( i == 12 ) {
				viz->addSphere(z1,0.0,1.0,0.0,0.125*0.5*3.14*pow(2.0,-level)*sphere_scale*3);
			}
		}
		*/
	}

	SpatialVector xhat(1,0,0), yhat(0,1,0), zhat(0,0,1);

	SpatialVector x0 = xhat;
	SpatialVector x1;
	float64 th=gPio2/900.0;
	SpatialRotation rot0 = SpatialRotation(zhat,th);
	if(true) {
		for(int i=0; i<901; ++i) {
			x1 = rot0.rotated_from(x0);
			viz->addEdge(
					x0,x1,
					0.0, 1.0, 0.0, a, scale);
			x0 = x1;
		}
	}

	float64 dlon = 1.0e-8;
	float64 lat0 = 0.0, lon0 = 0.0-dlon;
	level = 27;
	for(int i=0; i<10; ++i) {
		lon0 += dlon;
		STARE_ArrayIndexSpatialValue sid = index.ValueFromLatLonDegrees(lat0,lon0,level);
		STARE_ArrayIndexSpatialValues neighbors = index.NeighborsOfValue(sid);

		Triangle ta[12];
		for(int i=0; i<12; ++i) {
			ta[i] = index.TriangleFromValue(neighbors[i],level);
			Triangle tr0 = ta[i];
			viz->addEdge(
					tr0.vertices[0],tr0.vertices[1],
					0.0, 0.0, 1.0, a, scale);
			viz->addEdge(
					tr0.vertices[1],tr0.vertices[2],
					0.0, 0.0, 1.0, a, scale);
			viz->addEdge(
					tr0.vertices[2],tr0.vertices[0],
					0.0, 0.0, 1.0, a, scale);
		}
	}

	ok = true;
	return ok;
}

void PlotTriangles( VizHTM* viz, string s, int levelMax, STARE &index, float64 frac,
		float64 r, float64 g, float64 b, float64 a, float64 scale, bool leafOnly ) {

	int level = s.length() - 1;
	// cout << "PlotTriangles: s = " << s << ", level = " << level << endl << flush;

	stringstream ss;
	SpatialVector u0,u1,u2,u3;

	SpatialIndex sIndex = index.getIndex(level);

	for(int ic=0; ic < 4; ++ic) {
		ss.clear(); ss.str("");
		ss << s << ic;

		if( not leafOnly or level == levelMax ) {
			uint64 id = sIndex.idByName(ss.str().c_str());

			sIndex.nodeVertexByHtmId(u1,u2,u3,id);
			// sIndex.nodeVertex(id,u1,u2,u3);

			u1.normalize(); u2.normalize(); u3.normalize();

			u0 = u1 + u2 + u3; u0.normalize();

			int steps = 12;
			// viz->addArc(x0, x1, r, g, b, a, scale, steps)

			SpatialVector
			v1 = nudge(u1,u0,frac),
			v2 = nudge(u2,u0,frac),
			v3 = nudge(u3,u0,frac);

			viz->addArc(v1,v2,r,g,b,a,scale,steps);
			viz->addArc(v2,v3,r,g,b,a,scale,steps);
			viz->addArc(v3,v1,r,g,b,a,scale,steps);
		}

		if(level < levelMax) {
			PlotTriangles(viz,ss.str(),levelMax,index,frac,r,g,b,a,scale,leafOnly);
		}
	}
}

			/*
viz->addArc(u1,u2,1.0, 1.0, 1.0, a, scale, steps);
viz->addArc(u2,u3,1.0, 1.0, 1.0, a, scale, steps);
viz->addArc(u3,u2,1.0, 1.0, 1.0, a, scale, steps);
			 */

bool EquatorCheck2( VizHTM* viz ) {
	bool ok = false;

	// float64 sphere_scale = 0.025;
	float64 sphere_scale = 0.05;
	// float64 corner_nudge = 0.25;
	float64 corner_nudge = 0.01;
	// float64 corner_nudge = 0.0;

	// Color
	float64
	r = 0.25,
	g = 0.75,
	b = 0.25;

	// Transparency a;
	float64 a = 0.0;
	a = 0.5;

	// Scale the vectors
	float64 scale = 1.00001;

	int level = 27; // Not too bad with really small gEpsilon. Floating point gets really bad on interactive visualization

	STARE index;
	// STARE index1;

	SpatialVector xhat(1,0,0), yhat(0,1,0), zhat(0,0,1);
	if( true ){
		SpatialVector axis     = 0.5*xhat + 0.5*yhat; axis.normalize();
		float64       theta    = 0.25*gPi;
		// theta = 0.0;
		// theta    = 0.125*gPi;
		// theta    = 0.5*gPi + 0.00001;
		// theta    = 0.5*gPi; // Okay!
		SpatialRotation rotate_root_octahedron = SpatialRotation(axis,theta);
		int search_level = 27, build_level = 5;
		index = STARE(search_level, build_level, rotate_root_octahedron);
	}

	// int graphLevel = 0;
	// SpatialIndex sIndex = index.getIndex(graphLevel);
	/*
	corner_nudge = 0.0;
	PlotTriangles(viz,"N",3,index,corner_nudge, r, g, b, a, scale, true);
	PlotTriangles(viz,"S",3,index,corner_nudge, r, g, b, a, scale, true);
	*/

	/**/
	corner_nudge=0.0; r=1; g=1; b=0; scale = 1.000005;
	int graphLevel = 4;
	graphLevel = 5;
	// graphLevel = 6;
	cout << "graphLevel = " << graphLevel
			<< ", length(m) ~ " << index.lengthMeterScaleFromEdgeFromLevel(graphLevel) << endl << flush;
	PlotTriangles(viz,"N",graphLevel,index,corner_nudge, r, g, b, a, scale, true);
	PlotTriangles(viz,"S",graphLevel,index,corner_nudge, r, g, b, a, scale, true);
	/**/

	ok = true;
	return ok;
}

bool MeridianCheck1( VizHTM* viz ) {
	bool ok = false;

	// float64 sphere_scale = 0.025;
	float64 sphere_scale = 0.05;
	// float64 corner_nudge = 0.25;
	float64 corner_nudge = 0.05;
	// float64 corner_nudge = 0.0;

	// Transparency a;
	float64 a = 0.0;

	// Scale the vectors
	float64 nudge_eps = 1.0e-8;
	float64 scale = 1.0+nudge_eps;

	int level = 27; // Not too bad with really small gEpsilon. Floating point gets really bad on interactive visualization

	STARE index;
	// STARE index1;

	SpatialVector xhat(1,0,0), yhat(0,1,0), zhat(0,0,1);
	if( true ){
		SpatialVector axis     = 0.5*xhat + 0.5*yhat; axis.normalize();
		float64       theta    = 0.25*gPi;
		// theta = 0.0;
		// theta    = 0.125*gPi;
		// theta    = 0.5*gPi + 0.00001;
		// theta    = 0.5*gPi; // Okay!
		SpatialRotation rotate_root_octahedron = SpatialRotation(axis,theta);
		int search_level = 27, build_level = 5;
		index = STARE(search_level, build_level, rotate_root_octahedron);

		/*
		float64 dtheta = 1.0e-9;
		for(int i=0; i<21; ++i) {
			SpatialRotation rot = SpatialRotation(axis,i*dtheta);
			SpatialVector z1 = rot.rotated_from(zhat);
			if( i == 0 ) {
				viz->addSphere(z1,1.0,i/20.0,0.0,0.125*0.5*3.14*pow(2.0,-level)*sphere_scale*3);
			} else {
				viz->addSphere(z1,1.0,i/20.0,0.0,0.125*0.5*3.14*pow(2.0,-level)*sphere_scale);
			}
			if( i == 12 ) {
				viz->addSphere(z1,0.0,1.0,0.0,0.125*0.5*3.14*pow(2.0,-level)*sphere_scale*3);
			}
		}
		*/
	}

	SpatialVector x0 = xhat;
	SpatialVector x1;
	float64 th=gPio2/900.0;
	SpatialRotation rot0 = SpatialRotation(zhat,th);
	if(false) {
		// equator
		for(int i=0; i<901; ++i) {
			x1 = rot0.rotated_from(x0);
			viz->addEdge(
					x0,x1,
					0.0, 1.0, 0.0, a, scale);
			x0 = x1;
		}
	}

	float64 dlon = 1.0e-3;
	float64 lat0 = 45.0, lon0 = 0.0-dlon;

	SpatialVector graphicsOrigin(0.0,0.0,0.0);
	graphicsOrigin.setLatLonDegrees(lat0, lon0);
	float64       graphicsScale = 1.0;

	level = 27;
	for(int i=0; i<1; ++i) {
		lon0 += dlon;
		STARE_ArrayIndexSpatialValue sid = index.ValueFromLatLonDegrees(lat0,lon0,level);
		STARE_ArrayIndexSpatialValues neighbors = index.NeighborsOfValue(sid);

		{
			Triangle tr0 = index.TriangleFromValue(sid,level);

			for( int k=0; k<3; ++k ) {
				tr0.vertices[k] = tr0.vertices[k] - graphicsOrigin;
			}
			tr0.centroid = tr0.centroid - graphicsOrigin;

			viz->addEdge(
					tr0.vertices[0],tr0.vertices[1],
					1.0, 1.0, 1.0, a, scale);
			viz->addEdge(
					tr0.vertices[1],tr0.vertices[2],
					1.0, 1.0, 1.0, a, scale);
			viz->addEdge(
					tr0.vertices[2],tr0.vertices[0],
					1.0, 1.0, 1.0, a, scale);
			cout << " sid tr0.vert " << endl
					<< setprecision(17)
					<< setw(20)
					<< scientific
					<< tr0.vertices[0] << endl
					<< tr0.vertices[1] << endl
					<< tr0.vertices[2] << endl << flush;
		}

		Triangle ta[12];
		for(int i=0; i<12; ++i) {
			ta[i] = index.TriangleFromValue(neighbors[i],level);
			Triangle tr0 = ta[i];
			for( int k=0; k<3; ++k ) {
				tr0.vertices[k] = tr0.vertices[k] - graphicsOrigin;
			}
			viz->addEdge(
					tr0.vertices[0],tr0.vertices[1],
					0.0, 0.0, 1.0, a, scale);
			viz->addEdge(
					tr0.vertices[1],tr0.vertices[2],
					0.0, 1.0, 1.0, a, scale);
			viz->addEdge(
					tr0.vertices[2],tr0.vertices[0],
					1.0, 0.0, 1.0, a, scale);
			cout << i << " tr0.vert " << endl
					<< setprecision(17)
					<< setw(20)
					<< scientific
					<< tr0.vertices[0] << endl
					<< tr0.vertices[1] << endl
					<< tr0.vertices[2] << endl << flush;
		}
	}
	if(true) {
		SpatialVector x0,x1;
		float64 lat = lat0, lon = 0;
		x0.setLatLonDegrees(lat, lon);
		for(int i=0; i<2; ++i) {
		// for(int i=0; i<101; ++i) {
		// for(int i=0; i<901; ++i) {
			x1.setLatLonDegrees(lat, lon+i*90.0/900.0);
			viz->addEdge(
					x0 - graphicsOrigin,x1 - graphicsOrigin,
					0.0, 1.0, 0.0, a, scale);
			x0 = x1;
		}
	}

	ok = true;
	return ok;
}

