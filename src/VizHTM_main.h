/*
 * VizHTM_main.h
 *
 *  Created on: Dec 24, 2015
 *      Author: mrilee
 */

#ifndef VIZHTM_MAIN_H_
#define VIZHTM_MAIN_H_

void testTriaxis(VizHTM *viz);
void testAnEdge(VizHTM *viz);
void testConstraint(VizHTM *viz, int htmIdDepth);
void testConstraintAD(VizHTM *viz, SpatialVector a, float64 d);
void testConvexHtmRangeIntersection(VizHTM *viz, RangeConvex convex, int htmIdDepth);

void testHTMRange(VizHTM *viz, int htmIdLevel, const char *n0, const char *n1, bool annotation=true);

void testText1(VizHTM *viz,SpatialVector a);
void testShowEdgeProjections(VizHTM *viz, SpatialVector a, float64 d);
void testLatLonSpiral(VizHTM *viz);
void testIJKRGBFace(VizHTM *viz);

void testAddRectangle(VizHTM *viz, int htmIdDepth);

#endif /* VIZHTM_MAIN_H_ */
