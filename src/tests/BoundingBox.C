/*
 * BoundingBox1.C
 *
 *  Created on: Feb 27, 2019
 *      Author: mrilee
 *
 *  Copyright (C) 2019 Rilee Systems Technologies LLC
 */

#include "tests.h"

bool BoundingBox1(VizHTM *viz) {

	viz->addLatLonBoxEdgesDegrees(0,0,5,5,0.9,0.9,0.9);

	LatLonDegrees64ValueVector latlonbox;
	latlonbox.push_back(LatLonDegrees64(0,0));
	latlonbox.push_back(LatLonDegrees64(5,0));
	latlonbox.push_back(LatLonDegrees64(5,5));
	latlonbox.push_back(LatLonDegrees64(0,5));

	STARE index;
	// STARE_SpatialIntervals intervals = index.BoundingBoxFromLatLonDegrees(latlonbox,6);
	STARE_SpatialIntervals intervals = index.CoverBoundingBoxFromLatLonDegrees(latlonbox);
	EmbeddedLevelNameEncoding leftJustified;
	BitShiftNameEncoding      rightJustified;

	/*
	float
	color_scale = 1.0/intervals.size(),
	r0 = 0, g0 = 1, b0 = 0,
	r1 = 0, g1 = 1, b1 = 0,
	r2 = 0, g2 = 1, b2 = 0,
	a0 = 0, a1 = 0, a2 = 0,
	scale = 1;
	 */

	float
	color_scale = 1.0/intervals.size(),
	r0 = 1, g0 = 0, b0 = 0,
	r1 = 0, g1 = 1, b1 = 0,
	r2 = 0, g2 = 0, b2 = 1,
	a0 = 0, a1 = 0, a2 = 0,
	scale = 1;

	for( STARE_SpatialIntervals::iterator iSid = intervals.begin(); iSid != intervals.end(); ++iSid ) {
		Triangle tr0 = index.TriangleFromValue(*iSid);

		// leftJustified.setIdFromSciDBLeftJustifiedFormat(*iSid);
		// rightJustified.setId(leftJustified.rightJustifiedId());
		// id_line = rightJustified.getId();

		// cout << "centroid: " << tr0.centroid    << endl;
		// cout << "tr0.0:    " << tr0.vertices[0] << endl;
		// cout << "tr0.1:    " << tr0.vertices[1] << endl;
		// cout << "tr0.2:    " << tr0.vertices[2] << endl;

		viz->addFace3(tr0.vertices[0], tr0.vertices[1], tr0.vertices[2], r0, g0, b0, r1, g1, b1, r2, g2, b2, a0, a1, a2, scale);

		/*
		r0 += color_scale;
		r1 += color_scale;
		r2 += color_scale;

		g0 -= color_scale;
		g1 -= color_scale;
		g2 -= color_scale;

		b0 += color_scale;
		b1 += color_scale;
		b2 += color_scale;
		 */
	}

	return true;

}
