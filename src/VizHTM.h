/*
 * VizHTM.h
 *
 *  Created on: Oct 23, 2015
 *      Author: mrilee
 */

#ifndef VIZHTM_H_
#define VIZHTM_H_

#include "HtmRange.h"
#include "HstmRange.h"
#include "SpatialIndex.h"
#include "SpatialRotation.h"
#include "SpatialVector.h"
#include "SpatialConstraint.h"
#include "SpatialInterface.h"

#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSwitch.h>

#include <shapefil.h>

// TODO eliminate NARRAY_ also consider changing storage method for graphic elements.
// 64M 67108864
//  1M  1048576
// #define NARRAY_ 1048576
// #define NARRAY_ 4000000
#define NARRAY_ 100000000
// #define NARRAY_ 1000000000
//#define NARRAY_ 600000000
//#define NARRAY_ 200000000

typedef std::array<float,3>  Array3f;
typedef std::vector<Array3f> VectorArray3f;
// int, float
typedef std::vector<int>   Vector_i;
typedef std::vector<float> Vector_f;

struct annotation {
	SpatialVector *v;
	const char *text;
	float size;
	float r,g,b;
};
typedef std::vector<annotation> Vector_annotation;

class VizHTM {

	/**
	 * Provide visualization support for HTM or other related geometries.
	 */

public:
	VizHTM(int nArray);
	// TODO Finalization?

	void addEdge(
			const SpatialVector x0,
			const SpatialVector x1,
			float r, float g, float b, float a=-1., float scale=1.0,
			float deltaZ = 0.0
			);
	void addEdgeAndSphere(SpatialVector x0, SpatialVector x1, float r, float g, float b,
						  SpatialVector s0, float rs, float gs, float bs, float radius);

	void addEdge2(
			const SpatialVector x0, float r0, float g0, float b0, float a0, float scale0,
			const SpatialVector x1, float r1, float g1, float b1, float a1, float scale1
			);

	void addEdgeColor(float r, float g, float b, float a=-1.);
	void addFaceColor(float r, float g, float b, float a=-1.);
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
			float r2, float g2, float b2,
			float a0 = -1, float a1 = -1, float a2 = -1,
			float deltaZ = 0.0
			);

	void addFace3(
			SpatialVector x0,
			SpatialVector x1,
			SpatialVector x2,
			float r0, float g0, float b0,
			float r1, float g1, float b1,
			float r2, float g2, float b2,
			float a0 = -1, float a1 = -1, float a2 = -1,
			float scale = 1
			);

	void addFace4(
			float x00, float x01, float x02,
			float x10, float x11, float x12,
			float x20, float x21, float x22,
			float x30, float x31, float x32,
			float r0, float g0, float b0,
			float r1, float g1, float b1,
			float r2, float g2, float b2,
			float r3, float g3, float b3,
			float a0 = -1, float a1 = -1, float a2 = -1, float a3 = -1
			);
	void addFace4FromLatLonDegrees(
			float64 lat0, float64 lon0,
			float64 lat1, float64 lon1,
			float64 lat2, float64 lon2,
			float64 lat3, float64 lon3,
			float r0, float g0, float b0,
			float r1, float g1, float b1,
			float r2, float g2, float b2,
			float r3, float g3, float b3,
			float a0 = -1, float a1 = -1, float a2 = -1, float a3 = -1,
			float scale=1
			);

	void addRectangle(
			const SpatialVector x0,
			const SpatialVector x1,
			const SpatialVector x2,
			const SpatialVector x3,
			float r, float g, float b, float a = -1
			);

	void addLatLonBoxEdgesDegrees(
			float64 lat0, float64 lon0,
			float64 lat1, float64 lon1,
			float r, float g, float b
			);

	void addArc(
			const SpatialVector x0,
			const SpatialVector x1,
			float r, float g, float b, float a=-1., float scale=1., float deltaZ = 0.,
			int steps=-1);
	void addArcAtLatitudeDegrees(float64 lat, float64 lon0, float64 lon1, float r, float g, float b);

	void addEdgesFromIndexAndId(
			const SpatialIndex *index, uint64 htmId,
			float r, float g, float b, float a=-1., float scale=1.0);
	void addEdgesFromIndexAndName(
			SpatialIndex *index, const char* htmIdName,
			float r, float g, float b);
	void addEdgesFromIndexAndLatLonDegrees(
			SpatialIndex *index,
			float64 lat, float64 lon,
			float r, float g, float b
			);

	void addFaceFromIndexAndId(
			const SpatialIndex *index, uint64 htmId,
			float r0, float g0, float b0,
			float r1, float g1, float b1,
			float r2, float g2, float b2,
			float a=-1., float scale=1., float deltaZ = 0.
			);

	void addArcFromIndexAndId(
			const SpatialIndex *index, uint64 htmId,
			float r, float g, float b, float a=-1., float scale=1.0, float deltaZ = 0.0);
	void addArcFromIndexAndName(
			SpatialIndex *index, const char* htmIdName,
			float r, float g, float b, float a=-1.);
	void addArcFromIndexAndLatLonDegrees(
			SpatialIndex *index,
			float64 lat, float64 lon,
			float r, float g, float b, float a=-1.
			);

	void addArcsFromLatLonDegrees(
			float64 *lat, float64 *lon, int nPoints, bool close,
			float r, float g, float b, float a=-1., float deltaZ=0, int steps=-1
			);

	void addCellFromIndexAndId(
			const SpatialIndex *index, uint64 htmId,
			float r0, float g0, float b0, float a0, float zScale0,
			float r1, float g1, float b1, float a1, float zScale1,
			float hScale
			);

	void addConstraint(SpatialVector a, float64 d, float r, float g, float b);
	void addConstraint(SpatialConstraint c, float r, float g, float b);

	SoSeparator* makeText(SpatialVector *a, const char *text, float size, float r, float g, float b);
	// SoGroup* makeText(SpatialVector *a, const char *text, float size, float r, float g, float b);
	void addAnnotation(SpatialVector *a, const char *annotation, float size, float r, float g, float b);

	void addHTMInterval(SpatialIndex index, htmRange interval);
	void addHTMInterval(htmRange interval);

	void addHTMRange(
			const SpatialIndex *index, HtmRange *range,
			float r, float g, float b, float a=-1., float scale=1.0,
			bool arcFlag=true
			);
	void addHTMRange(
			HtmRange *range,
			float r, float g, float b, float a=-1., float scale=1.0,
			bool arcFlag=true
			);

	void addHstmRange(
			HstmRange *range,
			float r, float g, float b, float a=-1., float scale=1.0, bool arcFlag = true,
			float deltaZ = 0.0,
			SpatialIndex *index = NULL
	);
	void addCellsFromHstmRange(
			HstmRange *range,
			float r0, float g0, float b0, float a0, float scale0,
			float r1, float g1, float b1, float a1, float scale1,
			float hScale
	);

	void addHstmRangeIDs(
			HstmRange *range,
			float r, float g, float b, float a, float size, float scale, float dScale,
			bool arcFlag = true,
			bool edgeFlag = false

	);

	void addHstmRangeFaces(
			HstmRange *range,
			float r, float g, float b, float a=-1., float scale=1.0,
			float deltaZ = 0.0,
			SpatialIndex *index = NULL
	);

	void addCircleFacet(
			SpatialVector center,
			float64 halfSubtendedAngleInRadians,
			float r, float g, float b, float a=-1., float scale=1.0);

	void addCircleEdges(
			SpatialVector center,
			float64 halfSubtendedAngleInRadians,
			float r, float g, float b, float a=-1., float scale=1.0);

	void addShapeFile(
			string shapeFile, float r, float g, float b, float deltaZ = 0, bool verbose = false,
			int nStart = -1, int nEnd = -1);

	void debug_dump();

	void triaxis();
	SoSeparator* makeRoot();

	SoSwitch *faceSwitch;
	SoSwitch *edgeSwitch;

	int 	nArray;

	int 	nFaces;
	int   	nFaceColors;
	// float 	(*faceColors)[3]; // [nArray][3]
	VectorArray3f faceColors;

	// colors -> SbColor, coordinates -> SoCoordinate3
	// VectorArray3f       edgeColors, fCoordinates, sphereColors;

	// xIndices -> SoMFInt32
	// Vector_i            faceVertexColorIndices, edgeVertexColorIndices, edgeIndices, faceIndices;

	// SoSFFloat
	// Vector_f            fRadii;
	// Vector_annotations  annotations;

	float   (*faceTransparencies);
	bool    faceTransparency = false;
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

	string  projection = "None";
	bool    setProjection(string projection);
	string  getProjection();
	SpatialRotation projectionRotateLon = SpatialRotation(zhat,-90*gPr);

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
