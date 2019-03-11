/*
 * tests_support.C
 *
 *  Created on: Mar 11, 2019
 *      Author: mrilee
 *
 *  Copyright (C) 2019 Rilee Systems Technologies LLC
 */

#include "tests.h"

SpatialVector nudge(SpatialVector v0, SpatialVector v1, float64 a1) {
	SpatialVector u = (v0*(1-a1)) + (v1*a1); // u.normalize();
	return u;
}



