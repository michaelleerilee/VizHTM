/*
 * tests.h
 *
 *  Created on: Feb 27, 2019
 *      Author: mrilee
 *
 *  Copyright (C) 2019 Rilee Systems Technologies LLC
 */

#ifndef SRC_TESTS_TESTS_H_
#define SRC_TESTS_TESTS_H_

#include "VizHTM.h"
#include "STARE.h"
#include <iomanip>
#include <sstream>

bool BoundingBox1 (VizHTM *viz);

bool Edges1       (VizHTM *viz);

bool EquatorCheck1(VizHTM *viz);
bool EquatorCheck2(VizHTM *viz);
bool MeridianCheck1(VizHTM *viz);

bool PolePosition1(VizHTM *viz);
bool PoleCheck1   (VizHTM *viz);

SpatialVector nudge(SpatialVector v0, SpatialVector v1, float64 a1);

#endif /* SRC_TESTS_TESTS_H_ */
