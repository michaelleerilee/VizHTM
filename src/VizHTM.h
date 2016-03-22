/*
 * VizHTM.h
 *
 *  Created on: Oct 23, 2015
 *      Author: mrilee
 */

#ifndef VIZHTM_H_
#define VIZHTM_H_

#include "SpatialIndex.h"
#include "SpatialVector.h"
#include "SpatialConstraint.h"

#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>

// 64M 67108864
//  1M  1048576

#define NARRAY_ 1048576

class VizHTM {

	/**
	 * Provide visualization support for HTM or other related geometries.
	 */

public:
	VizHTM(int nArray);


	void addEdge(
			const SpatialVector x0,
			const SpatialVector x1,
			float r, float g, float b, float a=-1.);
	void addEdgeAndSphere(SpatialVector x0, SpatialVector x1, float r, float g, float b,
						  SpatialVector s0, float rs, float gs, float bs, float radius);

	void addEdgeColor(float r, float g, float b, float a=-1.);
	void addFaceColor(float r, float g, float b);
	void addCoordinate(float x, float y, float z);
	void addCoordinate64(float64 x, float64 y, float64 z);
	void addCoordinate(SpatialVector c);

	void addEdgeIndices(int i0, int i1);
	// Don't assume that adding an edge color requires it be used
	// by the edge list.
	void addEdgeVertexColorIndices(int i0, int i1);
	void addEdgeIndicesTriangle(int i0, int i1, int i2);
	void addEdgeProjections(SpatialVector x1);

	void addFaceIndices3(int i0, int i1, int i2);
	void addFaceVertexColorIndices3(int i0, int i1, int i2);
	void addSphere(int coordianteIndex, float r, float g, float b, float radius);
	void addSphere(SpatialVector x, float r, float g, float b, float radius);

	void addFace3(
			float x00, float x01, float x02,
			float x10, float x11, float x12,
			float x20, float x21, float x22,
			float r0, float g0, float b0,
			float r1, float g1, float b1,
			float r2, float g2, float b2);

	void addRectangle(
			const SpatialVector x0,
			const SpatialVector x1,
			const SpatialVector x2,
			const SpatialVector x3,
			float r, float g, float b);

	void addLatLonBoxEdgesDegrees(
			float64 lat0, float64 lon0,
			float64 lat1, float64 lon1,
			float r, float g, float b
			);

	void addArc(
			const SpatialVector x0,
			const SpatialVector x1,
			float r, float g, float b, float a=-1.);
	void addArcAtLatitudeDegrees(float64 lat, float64 lon0, float64 lon1, float r, float g, float b);

	void addEdgesFromIndexAndId(
			const SpatialIndex *index, uint64 htmId,
			float r, float g, float b, float a=-1.);
	void addEdgesFromIndexAndName(
			SpatialIndex *index, const char* htmIdName,
			float r, float g, float b);
	void addEdgesFromIndexAndLatLonDegrees(
			SpatialIndex *index,
			float64 lat, float64 lon,
			float r, float g, float b
			);

	void addArcFromIndexAndId(
			SpatialIndex *index, uint64 htmId,
			float r, float g, float b, float a=-1.);
	void addArcFromIndexAndName(
			SpatialIndex *index, const char* htmIdName,
			float r, float g, float b, float a=-1.);
	void addArcFromIndexAndLatLonDegrees(
			SpatialIndex *index,
			float64 lat, float64 lon,
			float r, float g, float b, float a=-1.
			);

	void addConstraint(SpatialVector a, float64 d, float r, float g, float b);
	void addConstraint(SpatialConstraint c, float r, float g, float b);

	SoSeparator* makeText(SpatialVector *a, const char *text, float size, float r, float g, float b);
	void addAnnotation(SpatialVector *a, const char *annotation, float size, float r, float g, float b);

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
	float   (*edgeTransparencies);
	bool    edgeTransparency = false;
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

	float   lineWidth = -1;
	float   sphereComplexity = -1;

	struct annotation {
		SpatialVector *v;
		const char *text;
		float size;
		float r,g,b;
	};

	int        nAnnotations;
	annotation *annotations;

};

SpatialVector* VectorFromLatLonDegrees(float lat, float lon);
SpatialVector* VectorFromLatLonRadians(float lat, float lon);

int rollDieWithFourSides();
double uniformDouble();
SpatialVector randomVector();
SpatialVector unitVector(SpatialVector x);

#endif /* VIZHTM_H_ */
