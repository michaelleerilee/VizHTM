/*
 * VizHTM.h
 *
 *  Created on: Oct 23, 2015
 *      Author: mrilee
 */

#ifndef VIZHTM_H_
#define VIZHTM_H_

#include "SpatialVector.h"

#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>

// 64M 67108864
//  1M  1048576

#define NARRAY_ 1048576

class VizHTM {

public:
	VizHTM(int nArray);

	void addEdgeColor(float r, float g, float b);
	void addFaceColor(float r, float g, float b);
	void addCoordinate(float x, float y, float z);
	void addCoordinate64(float64 x, float64 y, float64 z);
	void addEdgeIndices(int i0, int i1);
	// Don't assume that adding an edge color requires it be used
	// by the edge list.
	void addEdgeVertexColorIndices(int i0, int i1);
	void addEdgeIndicesTriangle(int i0, int i1, int i2);
	void addFaceIndices3(int i0, int i1, int i2);
	void addFaceVertexColorIndices3(int i0, int i1, int i2);
	void addSphere(int coordianteIndex, float r, float g, float b, float radius);

	void addFace3(
			float x00, float x01, float x02,
			float x10, float x11, float x12,
			float x20, float x21, float x22,
			float r0, float g0, float b0,
			float r1, float g1, float b1,
			float r2, float g2, float b2);

	void debug_dump();

	void triaxis();
	SoSeparator* makeRoot();

	SoSwitch *faceSwitch;
	SoSwitch *edgeSwitch;

	int 	nArray;

	int 	nFaces;
	int   	nFaceColors;
	float 	(*faceColors)[3]; // [nArray][3]
	int		nFaceVertexColorIndices;
	int		*faceVertexColorIndices;

	int   	nEdgeColors;
	float 	(*edgeColors)[3];
	int		nEdgeVertexColorIndices;
	int		*edgeVertexColorIndices;

	int 	nCoordinates;
	float 	(*fCoordinates)[3];

	int   	nSpheres;
	int   	*sphereIndices;
	float 	(*sphereColors)[3];
	float 	*fRadii;

	int 	nEdgeIndices;
	int   	*edgeIndices;

	int 	nFaceIndices;
	int   	*faceIndices;
};

unsigned int lg2(unsigned int v);
float* xyzFromLatLonRadians(float lat,float lon);
float* xyzFromLatLonDegrees(float lat,float lon);

#endif /* VIZHTM_H_ */
