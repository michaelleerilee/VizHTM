/*
 * VizHTM.C
 *
 *  Created on: Dec 15, 2015
 *      Author: mrilee
 */


#include "VizHTM.h"

#include <iostream>
#include <iomanip>

#include "SpatialException.h"
#include "SpatialIndex.h"
#include "SpatialVector.h"
#include "SpatialInterface.h"

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>

#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoIndexedLineSet.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoMaterialBinding.h>
#include <Inventor/nodes/SoSelection.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoSwitch.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/sensors/SoTimerSensor.h>
#include <Inventor/SbBasic.h>

using namespace std;

VizHTM::VizHTM(int nArray) : nArray(nArray) {

	nCoordinates 		= 0;

	nEdgeIndices 		= 0;
	nEdgeVertexColorIndices =0;
	nEdgeColors = 0;

	nFaces 				= 0;
	nFaceIndices 		= 0;
	nFaceVertexColorIndices =0;
	nFaceColors = 0;

	nSpheres    = 0;

	faceColors 	=  new float[nArray][3];
	faceVertexColorIndices = new int[nArray];

	edgeColors 	=  new float[nArray][3];
	edgeVertexColorIndices = new int[nArray];

	fCoordinates 	=  new float[nArray][3];

	sphereIndices 		= new int[nArray];
	sphereColors  	=  new float[nArray][3];
	fRadii        		= new float[nArray];

	for(int i = 0; i < nArray; i++) fRadii[i] = 0.;
	edgeIndices = new int[nArray];
	faceIndices = new int[nArray];

}

void VizHTM::addEdgeColor(float r, float g, float b){
	edgeColors[nEdgeColors][0] = r;
	edgeColors[nEdgeColors][1] = g;
	edgeColors[nEdgeColors][2] = b;
	nEdgeColors++;
}

void VizHTM::addFaceColor(float r, float g, float b){
	faceColors[nFaceColors][0] = r;
	faceColors[nFaceColors][1] = g;
	faceColors[nFaceColors][2] = b;
	nFaceColors++;
}

void VizHTM::addCoordinate(float x, float y, float z){
	fCoordinates[nCoordinates][0] = x;
	fCoordinates[nCoordinates][1] = y;
	fCoordinates[nCoordinates][2] = z;
	nCoordinates++;
}

// TODO float64 to float;  Mea maxima culpa.
void VizHTM::addCoordinate64(float64 x, float64 y, float64 z){
	cout << "  (adding-coord64 " << x << " " << y << " " << z << ") " << endl << flush;
	fCoordinates[nCoordinates][0] = (float) x;
	fCoordinates[nCoordinates][1] = (float) y;
	fCoordinates[nCoordinates][2] = (float) z;
	nCoordinates++;
}

void VizHTM::addEdgeIndices(int i0, int i1) {  // Generally i1 = i0 + 1
	edgeIndices[nEdgeIndices] = i0; nEdgeIndices++;
	edgeIndices[nEdgeIndices] = i1; nEdgeIndices++;
	edgeIndices[nEdgeIndices] = SO_END_LINE_INDEX; nEdgeIndices++;
}

void VizHTM::addEdgeIndicesTriangle(int i0, int i1, int i2) {  // Generally i1 = i0 + 1
	edgeIndices[nEdgeIndices] = i0; nEdgeIndices++;
	edgeIndices[nEdgeIndices] = i1; nEdgeIndices++;
	edgeIndices[nEdgeIndices] = SO_END_LINE_INDEX; nEdgeIndices++;
	edgeIndices[nEdgeIndices] = i1; nEdgeIndices++;
	edgeIndices[nEdgeIndices] = i2; nEdgeIndices++;
	edgeIndices[nEdgeIndices] = SO_END_LINE_INDEX; nEdgeIndices++;
	edgeIndices[nEdgeIndices] = i2; nEdgeIndices++;
	edgeIndices[nEdgeIndices] = i0; nEdgeIndices++;
	edgeIndices[nEdgeIndices] = SO_END_LINE_INDEX; nEdgeIndices++;
}

void VizHTM::addFaceIndices3(int i0, int i1, int i2) {  // Generally i2 = i1 + 1 = i0 + 2
	faceIndices[nFaceIndices] = i0;   nFaceIndices++;
	faceIndices[nFaceIndices] = i1; nFaceIndices++;
	faceIndices[nFaceIndices] = i2; nFaceIndices++;
	faceIndices[nFaceIndices] = SO_END_FACE_INDEX; nFaceIndices++;
}

void VizHTM::addFace3(
		float x00, float x01, float x02,
		float x10, float x11, float x12,
		float x20, float x21, float x22,
		float r0, float g0, float b0,
		float r1, float g1, float b1,
		float r2, float g2, float b2) {
	int indexBase = nCoordinates;
	int colorBase = nFaceVertexColorIndices;
	addFaceColor(r0,g0,b0);
	addFaceColor(r1,g1,b1);
	addFaceColor(r2,g2,b2);
	addFaceVertexColorIndices3(colorBase,colorBase+1,colorBase+2);
	addCoordinate( x00,  x01,  x02);
	addCoordinate( x10,  x11,  x12);
	addCoordinate( x20,  x21,  x22);
	addFaceIndices3(indexBase,indexBase+1,indexBase+2);
}

void VizHTM::addSphere(int coordinateIndex, float r, float g, float b, float radius) {
	sphereIndices[nSpheres]   = coordinateIndex;
	sphereColors[nSpheres][0] = r;
	sphereColors[nSpheres][1] = g;
	sphereColors[nSpheres][2] = b;
	fRadii[nSpheres]          = radius;
	nSpheres++;
}

void VizHTM::triaxis() {
	int indexBase = nCoordinates;
	int edgeColorBase = nEdgeColors;
	addEdgeColor(1.0,0.0,0.0);
	addEdgeColor(1.0,0.0,0.0);
	addEdgeVertexColorIndices(edgeColorBase,edgeColorBase+1);
	addCoordinate(0.0,0.0,0.0);
	addCoordinate(1.0,0.0,0.0);
	addEdgeIndices(indexBase,indexBase+1);

	indexBase = nCoordinates;
	edgeColorBase = nEdgeColors;
	addEdgeColor(1.0,0.0,0.0);
	addEdgeColor(1.0,0.0,0.0);
	addEdgeVertexColorIndices(edgeColorBase,edgeColorBase+1);
	addCoordinate(1.025,0.0,0.0);
	addCoordinate(1.5,0.0,0.0);
	addEdgeIndices(indexBase,indexBase+1);

	indexBase = nCoordinates;
	edgeColorBase = nEdgeColors;
	addEdgeColor(0.0,1.0,0.0);
	addEdgeColor(0.0,1.0,0.0);
	addEdgeVertexColorIndices(edgeColorBase,edgeColorBase+1);
	addCoordinate(0.0,0.0,0.0);
	addCoordinate(0.0,1.0,0.0);
	addEdgeIndices(indexBase,indexBase+1);

	indexBase = nCoordinates;
	edgeColorBase = nEdgeColors;
	addEdgeColor(0.0,1.0,0.0);
	addEdgeColor(0.0,1.0,0.0);
	addEdgeVertexColorIndices(edgeColorBase,edgeColorBase+1);
	addCoordinate(0.0,1.025,0.0);
	addCoordinate(0.0,1.5,0.0);
	addEdgeIndices(indexBase,indexBase+1);

	indexBase = nCoordinates;
	edgeColorBase = nEdgeColors;
	addEdgeColor(0.0,0.0,1.0);
	addEdgeColor(0.0,0.0,1.0);
	addEdgeVertexColorIndices(edgeColorBase,edgeColorBase+1);
	addCoordinate(0.0,0.0,0.0);
	addCoordinate(0.0,0.0,1.0);
	addEdgeIndices(indexBase,indexBase+1);

	indexBase = nCoordinates;
	edgeColorBase = nEdgeColors;
	addEdgeColor(0.0,0.0,1.0);
	addEdgeColor(0.0,0.0,1.0);
	addEdgeVertexColorIndices(edgeColorBase,edgeColorBase+1);
	addCoordinate(0.0,0.0,1.025);
	addCoordinate(0.0,0.0,1.5);
	addEdgeIndices(indexBase,indexBase+1);

	indexBase = nCoordinates;
	addCoordinate(1.05,0.,0.);
	addCoordinate(0.,1.05,0.);
	addCoordinate(0.,0.,1.05);
	addSphere(indexBase,1.,0.,0.,0.025);
	addSphere(indexBase+1,0.,1.,0.,0.025);
	addSphere(indexBase+2,0.,0.,1.,0.025);

//	indexBase = nCoordinates;
//	addCoordinate(0.,0.,1.05);
//	addSphere(indexBase,0.,0.,1.,0.025);
}

template <typename T>
void printA_(const char* prefix, const int n, const T v[], const char* postfix) {
	cout << prefix << flush;
	for(int i=0; i<n; i++) {
		cout << v[i] << " " << flush;
	}
	cout << postfix << flush;
}
template <typename T>
void print_(const char* prefix, const T v, const char* postfix) {
	cout << prefix << flush;
	cout << v << " " << flush;
	cout << postfix << flush;
}
void VizHTM::debug_dump() {
	cout << "debug_dumpStart" << endl << flush;

	cout << "nCoordinates: " << nCoordinates << endl;
	cout << "nFaceIndices: " << nFaceIndices << endl;
	cout << "nFaceColors:  " << nFaceColors  << endl;
	cout << "nEdgeIndices: " << nEdgeIndices << endl;
	cout << "nEdgeColors:  " << nEdgeColors  << endl;
	cout << "nEdgeVertexColorIndices: " << nEdgeVertexColorIndices << endl;

	for(int i=0;
			i<max(nFaces,
				max(nFaceColors,
				max(nEdgeColors,
				max(nCoordinates,
				max(nSpheres,
				max(nEdgeIndices,
				max(nFaceIndices,
				max(nEdgeVertexColorIndices,
					-1))))))));
			i++) {
		cout 	<< " i=" << i << " ";
		cout 	<< scientific;
		if(i<nEdgeIndices)	print_("e=",edgeIndices[i]," ");
		if(i<nEdgeVertexColorIndices) print_("eCi=",edgeVertexColorIndices[i]," ");
		if(i<nFaceIndices)	print_("f=( ",faceIndices[i],") ");
		if(i<nCoordinates)	printA_("c=( ",3,fCoordinates[i],") ");
		if(i<nEdgeColors) 	printA_("eColor=( ",3,edgeColors[i],") ");

		cout << endl << flush;
	}
	cout << "debug_dumpDone" << endl << flush;
}

void VizHTM::addEdgeVertexColorIndices(int i0, int i1) {
	edgeVertexColorIndices[nEdgeVertexColorIndices++]=i0;
	edgeVertexColorIndices[nEdgeVertexColorIndices++]=i1;
	edgeVertexColorIndices[nEdgeVertexColorIndices++]=SO_END_LINE_INDEX;
}

void VizHTM::addFaceVertexColorIndices3(int i0, int i1, int i2) {
	faceVertexColorIndices[nFaceVertexColorIndices++]=i0;
	faceVertexColorIndices[nFaceVertexColorIndices++]=i1;
	faceVertexColorIndices[nFaceVertexColorIndices++]=i2;
	faceVertexColorIndices[nFaceVertexColorIndices++]=SO_END_FACE_INDEX;
}

SoSeparator* VizHTM::makeRoot() {

//	if(true){ triaxis(); }// Try the triaxis. Incrementally add color and geometry...

	SoSeparator *root = new SoSeparator;

	SoCoordinate3 *coordinates = new SoCoordinate3;
	coordinates->point.setValues(0,nCoordinates,fCoordinates);
	root->addChild(coordinates);

	SoSeparator *edgeNode = new SoSeparator;

	SoMaterial *edgeMaterials = new SoMaterial;
	edgeMaterials->diffuseColor.setValues(0,nEdgeColors,edgeColors);
	edgeNode->addChild(edgeMaterials);

	SoMaterialBinding *edgeMaterialBinding = new SoMaterialBinding;
	edgeMaterialBinding->value = SoMaterialBinding::PER_VERTEX_INDEXED;
	edgeNode->addChild(edgeMaterialBinding);

	SoIndexedLineSet *edgeSet = new SoIndexedLineSet;
	edgeSet->coordIndex.setValues(0,nEdgeIndices,edgeIndices);
	edgeSet->materialIndex.setValues(0,nEdgeVertexColorIndices,edgeVertexColorIndices);
	edgeNode->addChild(edgeSet);

	edgeSwitch = new SoSwitch;
	edgeSwitch->whichChild = SO_SWITCH_ALL;
	edgeSwitch->addChild(edgeNode);
	root->addChild(edgeSwitch);


	SoSeparator *faceNode = new SoSeparator;

	SoMaterial *faceMaterials = new SoMaterial;
	faceMaterials->diffuseColor.setValues(0,nFaceColors,faceColors);
	faceNode->addChild(faceMaterials);

	SoMaterialBinding *faceMaterialBinding = new SoMaterialBinding;
	faceMaterialBinding->value = SoMaterialBinding::PER_VERTEX_INDEXED; // _INDEXED or not?
	faceNode->addChild(faceMaterialBinding);

	SoShapeHints* pHints   = new SoShapeHints;
	pHints->faceType       = SoShapeHints::UNKNOWN_FACE_TYPE;
	pHints->vertexOrdering = SoShapeHints::CLOCKWISE;
	faceNode->addChild(pHints);

	SoIndexedFaceSet *faceSet = new SoIndexedFaceSet;
	faceSet->coordIndex.setValues(0,nFaceIndices,faceIndices);
	faceSet->materialIndex.setValues(0,nFaceVertexColorIndices,faceVertexColorIndices);
	faceNode->addChild(faceSet);

	faceSwitch = new SoSwitch;
	faceSwitch->whichChild = SO_SWITCH_ALL;
	faceSwitch->addChild(faceNode);
	root->addChild(faceSwitch);

	if(nSpheres>0){ // If we added spheres, then add them to the graph.
		SoSeparator *sphereNodes = new SoSeparator;
		for(int iSphere = 0; iSphere < nSpheres; iSphere++) {
			SoSeparator *sphereNode = new SoSeparator;

			SoMaterial  *sphereMaterials = new SoMaterial;
			sphereMaterials->diffuseColor.setValue(sphereColors[iSphere]);
			sphereNode->addChild(sphereMaterials);

			SoTranslation *T = new SoTranslation;
			float x = fCoordinates[sphereIndices[iSphere]][0];
			float y = fCoordinates[sphereIndices[iSphere]][1];
			float z = fCoordinates[sphereIndices[iSphere]][2];
			T->translation.setValue(x,y,z);
			sphereNode->addChild(T);

			SoSphere *sp = new SoSphere();
			sp->radius = fRadii[iSphere];
			sphereNode->addChild(sp);

			sphereNodes->addChild(sphereNode);
		}
		root->addChild(sphereNodes);
	}

	root->unrefNoDelete();
	return root;
}



// http://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
// Find the log base 2 of an N-bit integer in O(lg(N)) operations
unsigned int
lg2(unsigned int v) {  // 32-bit value to find the log2 of
	const unsigned int b[] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000};
	const unsigned int S[] = {1, 2, 4, 8, 16};
	int i;

	register unsigned int r = 0; // result of log2(v) will go here
	for (i = 4; i >= 0; i--) // unroll for speed...
	{
		if (v & b[i])
		{
			v >>= S[i];
			r |= S[i];
		}
	}
	return r;
}
float* xyzFromLatLonRadians(float lat,float lon) {
	float *ret = new float[3];
	ret[0] = cos(lat)*cos(lon);
	ret[1] = cos(lat)*sin(lon);
	ret[2] = sin(lat);
	return ret;
}

float* xyzFromLatLonDegrees(float lat,float lon) {
	float piDiv180 = 4*atan(1.0f)/180.;
	return xyzFromLatLonRadians(lat*piDiv180,lon*piDiv180);
}



