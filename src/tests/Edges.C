/*
 * Edges.C
 *
 *  Created on: Feb 28, 2019
 *      Author: mrilee
 *
 *  Copyright (C) 2019 Rilee Systems Technologies LLC
 */

#include "tests.h"

bool Edges1( VizHTM *viz ) {

	SpatialVector origin(0,0,0);
	SpatialVector xhat(1,0,0);
	SpatialVector yhat(0,1,0);
	SpatialVector zhat(0,0,1);

	float64 scale = 1.0;
	float64 a     = -1;

	if(true) {

		viz->addEdge(origin, xhat, 1.0, 0.0, 0.0, a, scale);
		viz->addEdge(origin, yhat, 0.0, 1.0, 0.0, a, scale);
		viz->addEdge(origin, zhat, 0.0, 0.0, 1.0, a, scale);

		viz->addSphere(origin,1.0,1.0,1.0,0.125);
		viz->addSphere(xhat,1.0,0,0,0.125);
		viz->addSphere(yhat,0,1.0,0,0.125);
		viz->addSphere(zhat,0,0,1.0,0.125);

	}


	SpatialVector x0(0.1,0.1,0.1);
	SpatialVector x1(0.1,-0.1,0.1);
	for( int i = 0; i < 10; ++i ) {
		viz->addEdge(x0,x1,1.0,0.0,0.0,a,scale); scale += 0.1;
		viz->addEdge(x0,x1,0.0,1.0,0.0,a,scale); scale += 0.1;
		viz->addEdge(x0,x1,0.0,0.0,1.0,a,scale); scale += 0.1;
	}

	SpatialVector *xplus = new SpatialVector(xhat*1.2);
	viz->addAnnotation(xplus, "x\n", 0.25, 1.0, 0.0, 0.0);

	SpatialVector *zplus = new SpatialVector(zhat*1.2);
	viz->addAnnotation(zplus, "z\n", 0.25, 0.0, 0.0, 1.0);

	return true;
}
