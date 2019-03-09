/*
 * VizHTM.C
 *
 *  Created on: Dec 15, 2015
 *      Author: mrilee
 */


#include "VizHTM.h"

#include <iostream>
#include <iomanip>
#include <random>

#include "SpatialConstraint.h"
#include "SpatialException.h"
#include "SpatialIndex.h"
#include "SpatialVector.h"

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>

#include <Inventor/nodes/SoBaseColor.h>
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
#include <Inventor/nodes/SoTransparencyType.h>
#include <Inventor/sensors/SoTimerSensor.h>
#include <Inventor/SbBasic.h>

#include <Inventor/nodes/SoFont.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoText3.h>

#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoComplexity.h>

// #include <Inventor/nodes/.h>



using namespace std;

VizHTM::VizHTM(int nArray) : nArray(nArray) {

	nCoordinates 		    = 0;

	nEdgeIndices 		    = 0;
	nEdgeVertexColorIndices = 0;
	nEdgeColors             = 0;

	nFaces 				    = 0;
	nFaceIndices 		    = 0;
	nFaceVertexColorIndices = 0;
	nFaceColors             = 0;

	nSpheres                = 0;

	// faceColors 	            = new float[nArray][3];
	faceTransparencies      = new float[nArray];
	faceVertexColorIndices  = new int[nArray];

	edgeColors 	            = new float[nArray][3];
	edgeTransparencies      = new float[nArray];
	edgeVertexColorIndices  = new int[nArray];

	fCoordinates 	        = new float[nArray][3];

	sphereIndices 		    = new int[nArray];
	sphereColors  	        = new float[nArray][3];
	fRadii        		    = new float[nArray];

	for(int i = 0; i < nArray; i++) fRadii[i] = 0.;
	edgeIndices             = new int[nArray];
	faceIndices             = new int[nArray];

	nAnnotations            = 0;
	annotations             = new annotation[nArray];


}

void VizHTM::addEdgeColor(float r, float g, float b, float a){
	edgeColors[nEdgeColors][0] = r;
	edgeColors[nEdgeColors][1] = g;
	edgeColors[nEdgeColors][2] = b;
	if(a == -1.) {
		edgeTransparencies[nEdgeColors] = 0.; // The default case.
	} else {
		edgeTransparencies[nEdgeColors] = a;
		edgeTransparency = true; // If we ever see a transparency, turn it on.
	}
	nEdgeColors++;
}

void VizHTM::addFaceColor(float r, float g, float b,float a){
	/*
	faceColors[nFaceColors][0] = r;
	faceColors[nFaceColors][1] = g;
	faceColors[nFaceColors][2] = b;
	*/
	Array3f color; color[0] = r; color[1] = g; color[2] = b;
	faceColors.push_back(color);

	if(a == -1.) {
		faceTransparencies[nFaceColors] = 0.; // The default case.
	} else {
		faceTransparencies[nFaceColors] = a;
		faceTransparency = true; // If we ever see a transparency, turn it on.
	}
	nFaceColors++;
}

void VizHTM::addCoordinate(float x, float y, float z){
	fCoordinates[nCoordinates][0] = x;
	fCoordinates[nCoordinates][1] = y;
	fCoordinates[nCoordinates][2] = z;
	nCoordinates++;
	if(nCoordinates> nArray){
		cout << "VizHTM::addCoordinate overflow ERROR" << endl << flush;
	}
}

// TODO float64 to float;  Mea maxima culpa.
void VizHTM::addCoordinate64(float64 x, float64 y, float64 z){
//	cout << "  (adding-coord64 " << x << " " << y << " " << z << ") " << endl << flush;
	fCoordinates[nCoordinates][0] = (float) x;
	fCoordinates[nCoordinates][1] = (float) y;
	fCoordinates[nCoordinates][2] = (float) z;
	nCoordinates++;
	if(nCoordinates> nArray){
		cout << "VizHTM::addCoordinate64 overflow ERROR" << endl << flush;
	}
}

void VizHTM::addCoordinate(const SpatialVector c) {
	addCoordinate64(c.x(),c.y(),c.z());
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
	if(nEdgeIndices> nArray){
		cout << "VizHTM::addEdgeIndicesTriangle overflow ERROR" << endl << flush;
	}
}

void VizHTM::addFaceIndices3(int i0, int i1, int i2) {  // Generally i2 = i1 + 1 = i0 + 2
	faceIndices[nFaceIndices] = i0; nFaceIndices++;
	faceIndices[nFaceIndices] = i1; nFaceIndices++;
	faceIndices[nFaceIndices] = i2; nFaceIndices++;
	faceIndices[nFaceIndices] = SO_END_FACE_INDEX; nFaceIndices++;
	if(nFaceIndices> nArray){
		cout << "VizHTM::addFaceIndices3 overflow ERROR nFaceIndices=" << nFaceIndices << endl << flush;
		exit(1);
	}
	// Should we have nFaces++ here instead?
}

void VizHTM::addFace3(
		float x00, float x01, float x02,
		float x10, float x11, float x12,
		float x20, float x21, float x22,
		float r0, float g0, float b0,
		float r1, float g1, float b1,
		float r2, float g2, float b2,
		float a0, float a1, float a2
		) {
	int indexBase = nCoordinates;
	int colorBase = nFaceVertexColorIndices;
	addFaceColor(r0,g0,b0,a0);
	addFaceColor(r1,g1,b1,a1);
	addFaceColor(r2,g2,b2,a2);
	addFaceVertexColorIndices3(colorBase,colorBase+1,colorBase+2);
	addCoordinate( x00,  x01,  x02);
	addCoordinate( x10,  x11,  x12);
	addCoordinate( x20,  x21,  x22);
	addFaceIndices3(indexBase,indexBase+1,indexBase+2);
	nFaces++;
}

void VizHTM::addFace3(
		SpatialVector x0,
		SpatialVector x1,
		SpatialVector x2,
		float r0, float g0, float b0,
		float r1, float g1, float b1,
		float r2, float g2, float b2,
		float a0, float a1, float a2,
		float scale
		) {
	x0 = x0*scale;
	addFace3(
			x0.x(), x0.y(), x0.z(),
			x1.x(), x1.y(), x1.z(),
			x2.x(), x2.y(), x2.z(),
			r0, g0, b0,
			r1, g1, b1,
			r2, g2, b2,
			a0, a1, a2
			);
}

void VizHTM::addFace4(
		float x00, float x01, float x02,
		float x10, float x11, float x12,
		float x20, float x21, float x22,
		float x30, float x31, float x32,
		float r0, float g0, float b0,
		float r1, float g1, float b1,
		float r2, float g2, float b2,
		float r3, float g3, float b3,
		float a0, float a1, float a2, float a3
		) {
	addFace3(
			x00,x01,x02,
			x10,x11,x12,
			x20,x21,x22,
			r0, g0, b0,
			r1, g1, b1,
			r2, g2, b2,
			a0, a1, a2
	);
	addFace3(
			x00,x01,x02,
			x20,x21,x22,
			x30,x31,x32,
			r0, g0, b0,
			r2, g2, b2,
			r3, g3, b3,
			a0, a2, a3
			);
}
void VizHTM::addFace4FromLatLonDegrees(
		float64 lat0, float64 lon0,
		float64 lat1, float64 lon1,
		float64 lat2, float64 lon2,
		float64 lat3, float64 lon3,
		float r0, float g0, float b0,
		float r1, float g1, float b1,
		float r2, float g2, float b2,
		float r3, float g3, float b3,
		float a0, float a1, float a2, float a3,
		float scale
		) {
	SpatialVector *x0 = VectorFromLatLonDegrees(lat0,lon0);
	SpatialVector *x1 = VectorFromLatLonDegrees(lat1,lon1);
	SpatialVector *x2 = VectorFromLatLonDegrees(lat2,lon2);
	SpatialVector *x3 = VectorFromLatLonDegrees(lat3,lon3);
	x0->normalize(scale);
	x1->normalize(scale);
	x2->normalize(scale);
	x3->normalize(scale);

	addFace4(
			x0->x(),x0->y(),x0->z(),
			x1->x(),x1->y(),x1->z(),
			x2->x(),x2->y(),x2->z(),
			x3->x(),x3->y(),x3->z(),
			r0, g0, b0,
			r1, g1, b1,
			r2, g2, b2,
			r3, g3, b3,
			a0, a1, a2, a3
			);
	delete x0,x1,x2,x3;
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
	addSphere(indexBase,.5,0.,0.,0.025);
	addSphere(indexBase+1,0.,.5,0.,0.025);
	addSphere(indexBase+2,0.,0.,.5,0.025);

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
	/*
	edgeVertexColorIndices[nEdgeVertexColorIndices++]=i0;
	edgeVertexColorIndices[nEdgeVertexColorIndices++]=i1;
	edgeVertexColorIndices[nEdgeVertexColorIndices++]=SO_END_LINE_INDEX;
	*/
	edgeVertexColorIndices[nEdgeVertexColorIndices]=i0; nEdgeVertexColorIndices++;
	edgeVertexColorIndices[nEdgeVertexColorIndices]=i1; nEdgeVertexColorIndices++;
	edgeVertexColorIndices[nEdgeVertexColorIndices]=SO_END_LINE_INDEX; nEdgeVertexColorIndices++;
}

void VizHTM::addFaceVertexColorIndices3(int i0, int i1, int i2) {
	faceVertexColorIndices[nFaceVertexColorIndices++]=i0;
	faceVertexColorIndices[nFaceVertexColorIndices++]=i1;
	faceVertexColorIndices[nFaceVertexColorIndices++]=i2;
	faceVertexColorIndices[nFaceVertexColorIndices++]=SO_END_FACE_INDEX;
}

void VizHTM::addConstraint(SpatialVector a, float64 d, float r, float g, float blue) {
	SpatialVector ad = a*d;
	SpatialVector o  = SpatialVector(0.,0.,0.);
	SpatialVector as = 1.025*ad;
	addEdgeAndSphere(o,ad,r,g,blue,as,r,g,blue,0.025*sqrt(1.0-d*d));

	SpatialVector i, a_cross_i;
	float aci_norm;
	do {
		i = randomVector();
//		cout << "i: " << i << endl << flush;
		a_cross_i = a^i;
		aci_norm = a_cross_i.length();
	} while (!aci_norm);

	SpatialVector b = a_cross_i; b.normalize();
	SpatialVector c = a^b; c.normalize();
//	viz->addEdgeAndSphere(o,b,0.,1.,0.5*d,b,0.,1.,0.5*d,0.025);
//	viz->addEdgeAndSphere(o,c,0.,0.5*d,1.,c,0.,0.5*d,1.,0.025);

	float64 theta = 0.;
	float64 twopi  = atan2(0.,-0.);
	float64 deltaTheta = 0.05;
	SpatialVector dlast = sin(-twopi*(1.+deltaTheta))*b + cos(-twopi*(1.+deltaTheta))*c; dlast.normalize();
	dlast *= sqrt(1.0-d*d);
	SpatialVector last = ad + dlast;
	for(float64 theta = -1.; theta<=1.; theta += deltaTheta) {
		float64 lambda = sin(twopi*theta);
		float64 mu     = cos(twopi*theta);
		SpatialVector delta_ad = lambda*b + mu*c; delta_ad.normalize();
		delta_ad *= sqrt(1.0-d*d);
		SpatialVector aPlusDelta = ad + delta_ad;
//			viz->addEdge(ad,aPlusDelta,1.,0.5*d,1.);
//			viz->addEdge(o,aPlusDelta,1.,0.5*d,1.);
		addEdge(last,aPlusDelta,r,g,blue);
		last = aPlusDelta;
	}

}

void VizHTM::addConstraint(SpatialConstraint c, float r, float g, float b) {
		addConstraint(c.v(),c.d(),r,g,b);
}

SoSeparator* VizHTM::makeRoot() {

//	if(true){ triaxis(); }// Try the triaxis. Incrementally add color and geometry...

	cout << "makeRoot" << endl
			<< "nArray                  " << nArray << endl
			<< "nFaces                  " << nFaces << endl
			<< "nFaceColors             " << nFaceColors << endl
			<< " faceTransparency       " << faceTransparency << endl
			<< "nFaceVertexColorIndices " << nFaceVertexColorIndices << endl
			<< "nEdgeColors             " << nEdgeColors << endl
			<< " edgeTransparency       " << edgeTransparency << endl
			<< "nEdgeVertexColorIndices " << nEdgeVertexColorIndices << endl
			<< "nCoordinates            " << nCoordinates << endl
			<< "nSpheres                " << nSpheres << endl
			// << "nEdges                  " << nEdges << endl
			<< "nEdgeIndices            " << nEdgeIndices << endl
			<< "nFaceIndices            " << nFaceIndices << endl
			<< "nAnnotations            " << nAnnotations << endl;

	SoSeparator *root = new SoSeparator;

	if(lineWidth != -1) {
		SoDrawStyle *style = new SoDrawStyle;
		style->lineWidth.setValue(lineWidth);
		root->addChild(style);
	}

	if( true ) {
	// Broken
//	if(faceTransparency || edgeTransparency) {
//		SoTransparencyType *transparencyType = new SoTransparencyType;
//		transparencyType->value = SoTransparencyType::DELAYED_BLEND;
//		root->addChild(transparencyType);
//	}

	if( nCoordinates > 0 ) {

	SoCoordinate3 *coordinates = new SoCoordinate3;
	coordinates->point.setValues(0,nCoordinates,fCoordinates);
	root->addChild(coordinates);

	SoSeparator *edgeNode = new SoSeparator;

	SoMaterial *edgeMaterials = new SoMaterial;
	edgeMaterials->diffuseColor.setValues(0,nEdgeColors,edgeColors);
	// edgeMaterials->diffuseColor.setValues(0,nEdgeColors,(SbColor*)&edgeColors[0]);
	if(edgeTransparency) {
		edgeMaterials->transparency.setValues(0,nEdgeColors,edgeTransparencies);
	}
	edgeNode->addChild(edgeMaterials);

	/*
	for( int i=0; i < nEdgeColors; ++i ) {
		const SbColor *color = edgeMaterials->diffuseColor.getValues(i);
		cout << i << " color: " << color->toString().getString() << endl << flush;
	}
	*/

	SoMaterialBinding *edgeMaterialBinding = new SoMaterialBinding;
	// edgeMaterialBinding->value = SoMaterialBinding::PER_PART_INDEXED;
	// edgeMaterialBinding->value = SoMaterialBinding::PER_PART;
	edgeMaterialBinding->value = SoMaterialBinding::PER_VERTEX_INDEXED;
	// edgeMaterialBinding->value = SoMaterialBinding::PER_VERTEX;
	// edgeMaterialBinding->value = SO_END_LINE_INDEX;

	edgeNode->addChild(edgeMaterialBinding);

	SoIndexedLineSet *edgeSet = new SoIndexedLineSet;
	edgeSet->coordIndex.setValues(0,nEdgeIndices,edgeIndices);
	edgeSet->materialIndex.setValues(0,nEdgeVertexColorIndices,edgeVertexColorIndices);
	edgeNode->addChild(edgeSet);

	/*
	cout << " num mat index: " << edgeSet->materialIndex.getNum() << endl << flush;
	for( int i=0; i < edgeSet->materialIndex.getNum(); ++i ) {
		cout << i << ", mat index " << *(edgeSet->materialIndex.getValues(i));
		if( *(edgeSet->materialIndex.getValues(i)) >= 0 ) {
			int idx = *(edgeSet->materialIndex.getValues(i));
			cout << " " << edgeMaterials->diffuseColor.getValues(idx)->toString().getString();
		}
		cout << endl << flush;
	}
	*/

	edgeSwitch = new SoSwitch;
	edgeSwitch->whichChild = SO_SWITCH_ALL;
	edgeSwitch->addChild(edgeNode);
	root->addChild(edgeSwitch);

	SoSeparator *faceNode = new SoSeparator;

	SoMaterialBinding *faceMaterialBinding = new SoMaterialBinding;
	// faceMaterialBinding->value = SoMaterialBinding::PER_VERTEX_INDEXED; // _INDEXED or not?
	faceMaterialBinding->value = SoMaterialBinding::PER_VERTEX; // _INDEXED or not? // Why? Why? Why? // TODO WHY?
	faceNode->addChild(faceMaterialBinding);

	SoMaterial *faceMaterials = new SoMaterial;
	// faceMaterials->diffuseColor.setValues(0,nFaceColors,faceColors);
	faceMaterials->diffuseColor.setValues(0,nFaceColors,(SbColor*)&faceColors[0]);
	if(faceTransparency) {
		faceMaterials->transparency.setValues(0,nFaceColors,faceTransparencies);
	}
	faceNode->addChild(faceMaterials);

//	cout << "1000 " << faceMaterials->transparency.getNum() << endl << flush;

	SoShapeHints* pHints   = new SoShapeHints;
// Working
	pHints->faceType       = SoShapeHints::UNKNOWN_FACE_TYPE;
	pHints->vertexOrdering = SoShapeHints::CLOCKWISE;

// Working? Unknown
//	pHints->vertexOrdering = SoShapeHints::UNKNOWN_ORDERING;
//	pHints->vertexOrdering = SoShapeHints::COUNTERCLOCKWISE;
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

		SoComplexity *complexity = new SoComplexity;
		if(sphereComplexity != -1) {
			complexity->value.setValue(sphereComplexity);
		}
		sphereNodes->addChild(complexity);

		for(int iSphere = 0; iSphere < nSpheres; iSphere++) {
			SoSeparator *sphereNode = new SoSeparator;

			SoMaterial  *sphereMaterials = new SoMaterial;
			sphereMaterials->diffuseColor.setValue(sphereColors[iSphere]);
			sphereMaterials->ambientColor.setValue(sphereColors[iSphere]);
//			sphereMaterials->specularColor.setValue(sphereColors[iSphere]);
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
	} // if nCoordinates > 0

//	{
//		cout << "Starting text test" << flush;
//		SpatialVector a = SpatialVector(0.0,0.0,1.0);
//		float size = 1.0;
//		float r = 0.8, g = 0.8, b = 0.9;
//		addText(a,"Test",size,r,g,b);
//		cout << "...done." << endl << flush;
//	}

	} // top

	if(nAnnotations>0){
		SoSeparator *texts = new SoSeparator;
		for(int ia=0;ia<nAnnotations;ia++) {
			/*
			cout
				<< " adding annotation: " << ia
				<< " text: " << annotations[ia].text
				<< endl << flush;
			 */
			texts->addChild(makeText(
					annotations[ia].v,
					annotations[ia].text,
					annotations[ia].size,
					annotations[ia].r,
					annotations[ia].g,
					annotations[ia].b));
		}
		root->addChild(texts);
	} else {
		cout << " no annotations..." << endl << flush;
	}

	// root->unrefNoDelete();
	return root;
}


SpatialVector* VectorFromLatLonRadians(float lat, float lon) {
	float64 *x = xyzFromLatLonRadians(lat,lon);
	SpatialVector *ret = new SpatialVector(x[0],x[1],x[2]); ret->normalize();
	return ret;
}

SpatialVector* VectorFromLatLonDegrees(float lat, float lon) {
	float piDiv180 = 4*atan(1.0f)/180.;
	return VectorFromLatLonRadians(lat*piDiv180,lon*piDiv180);
}

int rollDieWithFourSides() {
	static std::default_random_engine e{};
	static std::uniform_int_distribution<int> d{0,3};
	return d(e);
}

double uniformDouble() {
	static std::default_random_engine e{};
	static std::uniform_real_distribution<double> d{0,1};
	//uniform_int_distribution<int> d{0,3};
	return d(e);
}

SpatialVector randomVector() {
	SpatialVector r = SpatialVector(
			uniformDouble(),
			uniformDouble(),
			uniformDouble()
			);
	r.normalize();
	return r;
}

void VizHTM::addEdge(
		const SpatialVector x0,
		const SpatialVector x1,
		float r, float g, float b, float a, float scale) {

	// float scale_ = (scale<0) ? 1.0 : scale;

	//cout << " addEdge: start, " << flush;
	int colorBase = nEdgeColors;
	addEdgeColor(r,g,b,a);
	addEdgeColor(r,g,b,a);
	addEdgeVertexColorIndices(colorBase,colorBase+1);

	int coordBase = nCoordinates;
//	SpatialVector a = x0*scale;
//	SpatialVector b = x1*scale;
//	addCoordinate(a); addCoordinate(b);
	addCoordinate(x0*scale);
	addCoordinate(x1*scale);
	addEdgeIndices(coordBase,coordBase+1);
	//cout << " addEdge: done; " << flush;
}

void VizHTM::addEdge2(
		const SpatialVector x0, float r0, float g0, float b0, float a0, float scale0,
		const SpatialVector x1, float r1, float g1, float b1, float a1, float scale1
		) {

	// float scale_ = (scale<0) ? 1.0 : scale;

	//cout << " addEdge: start, " << flush;
	int colorBase = nEdgeColors;
	addEdgeColor(r0,g0,b0,a0);
	addEdgeColor(r1,g1,b1,a1);
	addEdgeVertexColorIndices(colorBase,colorBase+1);

	int coordBase = nCoordinates;
//	SpatialVector a = x0*scale;
//	SpatialVector b = x1*scale;
//	addCoordinate(a); addCoordinate(b);
	addCoordinate(x0*scale0);
	addCoordinate(x1*scale1);
	addEdgeIndices(coordBase,coordBase+1);
	//cout << " addEdge: done; " << flush;
}

void VizHTM::addEdgeProjections(SpatialVector x1) {
	SpatialVector o = SpatialVector(0.,0.,0.);
	SpatialVector p = SpatialVector(x1.x(),x1.y(),0.);
//	SpatialVector z = SpatialVector(0.,0.,x1.z());
	addEdge(o,p,1.,1.,0.);
	addEdge(p,x1,0.,0.,1.);
}

void VizHTM::addEdgeAndSphere(SpatialVector x0, SpatialVector x1, float r, float g, float b, SpatialVector s0, float rs, float gs, float bs,
		float radius) {
	addEdge(x0,x1,r,g,b);
	addSphere(s0,rs,gs,bs,radius);
}

void VizHTM::addSphere(int coordinateIndex, float r, float g, float b, float radius) {
	sphereIndices[nSpheres]   = coordinateIndex;
	sphereColors[nSpheres][0] = r;
	sphereColors[nSpheres][1] = g;
	sphereColors[nSpheres][2] = b;
	fRadii[nSpheres]          = radius;
	nSpheres++;
}

void VizHTM::addSphere(SpatialVector x, float r, float g, float b, float radius) {
	int coordinateBase = nCoordinates;
	addCoordinate(x);
	addSphere(coordinateBase,r,g,b,radius);

}

SoSeparator* VizHTM::makeText(SpatialVector* a, const char* annotation, float size, float r, float g, float b) {
	// cf. http://oivdoc90.vsg3d.com/content/62-three-dimensional-text

	/*
	cout << " makeText - adding " << annotation << endl << flush;
	cout << "       size,r,g,b: " << size << ", " << r << ", " << g << ", " << b << endl << flush;
	cout << "                v: " << (*a) << endl << flush;
	 */

	SoSeparator *root = new SoSeparator;
	// SoGroup *root = new SoGroup;
	SoSeparator *textSeparator = new SoSeparator;

	SoFont *font = new SoFont;
	// font->name.setValue("Times-Roman");
	// font->name.setValue("Cantarell-Regular");
	font->name.setValue("TGS_Complex_Roman");
	font->size.setValue(size);
	root->addChild(font);

	/**/
	SoMaterial *material = new SoMaterial;
	SoMaterialBinding *binding = new SoMaterialBinding;
	material->diffuseColor.set1Value(0,SbColor(r,g,b));
	material->diffuseColor.set1Value(1,SbColor(0.5*r,0.5*g,0.5*b));
	binding->value = SoMaterialBinding::PER_PART;
	root->addChild(material);
	root->addChild(binding);
	/**/

	SoBaseColor *baseColor = new SoBaseColor();
	baseColor->rgb.setValue(SbColor(r,g,b));
	root->addChild(baseColor);

	SoText3 *text = new SoText3;

	SoTransform *rot0 = new SoTransform;
	SoTransform *tra0 = new SoTransform;

	SpatialVector i = SpatialVector(1.0,0.0,0.0);
	SpatialVector k = SpatialVector(0.0,0.0,1.0);

	SpatialVector axis = k^(*a);
	tra0->rotation.setValue(SbVec3f(axis.x(),axis.y(),axis.z()),acos(k*(*a)));
//	transform->rotation.setValue(SbVec3f(a->x(),a->y(),a->z()),0.1);
	tra0->translation.setValue(a->x(),a->y(),a->z());

	rot0->rotation.setValue(SbVec3f(a->x(),a->y(),a->z()),acos(i*(*a)));

	if( true ) {
		text->parts = SoText3::FRONT;
		//  text->parts = SoText3::ALL;
		string *s = new string(annotation);
		string delimiter = "\n";
		int idx=0; int last=0; int next=0; while((next=s->find(delimiter,last)) != string::npos) {
			text->string.set1Value(idx++,(s->substr(last,next-last)).c_str());
			last = next + 1;
		}
	} else {
	// text->string = "Test";
		text->string = annotation;
	}

	/**/
	textSeparator->addChild(tra0);
	textSeparator->addChild(rot0);
	textSeparator->addChild(text);
	/**/

	/*
	root->addChild(tra0);
	root->addChild(rot0);
	root->addChild(text);
	*/

	root->addChild(textSeparator);

	return root;
}

void VizHTM::addAnnotation(SpatialVector *a, const char *annotation, float size, float r, float g, float b){
//	cout << " vhtm::addA::nA: " << nAnnotations << flush;
	annotations[nAnnotations].r = r;
	annotations[nAnnotations].g = g;
	annotations[nAnnotations].b = b;
	annotations[nAnnotations].size = size;
	annotations[nAnnotations].text = annotation;
	annotations[nAnnotations].v = a;
	nAnnotations++;
}

SpatialVector unitVector(SpatialVector x) {
	SpatialVector n = x; n.normalize();
	return n;
}

void VizHTM::addRectangle(
		const SpatialVector x0,
		const SpatialVector x1,
		const SpatialVector x2,
		const SpatialVector x3,
		float r, float g, float b, float a
		) {
	addArc(x0,x1,r,g,b,a);
	addArc(x1,x2,r,g,b,a);
	addArc(x2,x3,r,g,b,a);
	addArc(x3,x0,r,g,b,a);
}

void VizHTM::addLatLonBoxEdgesDegrees(
		float64 lat0, float64 lon0,
		float64 lat1, float64 lon1,
		float r, float g, float b
		) {
	addArcAtLatitudeDegrees(lat0,lon0,lon1,r,g,b);
	addArcAtLatitudeDegrees(lat1,lon0,lon1,r,g,b);
	addArc( *VectorFromLatLonDegrees(lat0,lon0),
			*VectorFromLatLonDegrees(lat1,lon0),
			r,g,b
			);
	addArc( *VectorFromLatLonDegrees(lat0,lon1),
			*VectorFromLatLonDegrees(lat1,lon1),
			r,g,b
			);
}

void VizHTM::addArc(
		const SpatialVector x0,
		const SpatialVector x1,
		float r, float g, float b, float alpha, float scale,
		int steps) {

	SpatialVector dx = x1 - x0;
	SpatialVector *last = new SpatialVector(x0);

	if(steps<0) {
		steps = 20;
		if(dx.length()>0.5){
			steps = 180;
		}
	}

//	cout << "steps: " << steps << flush;

	for(double i=1; i<steps-1; i++) {
		double a = i / (steps - 1);
//		cout << i << ", " << a << "; " << flush;
		SpatialVector *next = new SpatialVector(x0 + a * dx); next->normalize();
		addEdge(*last,*next, r, g, b, alpha, scale );
		last = next;
	}

//	cout << endl << flush;

	addEdge(*last,x1,r,g,b,alpha, scale);
}

void VizHTM::addArcAtLatitudeDegrees(float64 lat, float64 lon0, float64 lon1, float r, float g, float b) {
	SpatialVector *x0 = VectorFromLatLonDegrees(lat,lon0);
	SpatialVector *x1 = VectorFromLatLonDegrees(lat,lon1);

	const int steps = 20;
	double dlon = lon1 - lon0;

	SpatialVector *last = VectorFromLatLonDegrees(lat,lon0);
	for(double i = 0; i<steps-1; i++) {
		double a = i / (steps - 1);
		SpatialVector *next = VectorFromLatLonDegrees(lat,lon0+a*dlon); next->normalize();
		addEdge(*last,*next, r, g, b);
		last = next;
	}
	addEdge(*last,*x1,r,g,b);

}

void VizHTM::addEdgesFromIndexAndId(
		const SpatialIndex *index, uint64 htmId,
		float r, float g, float b, float a, float scale
		) {
	SpatialVector v0, v1, v2;
	uint64 nodeIndex = index->nodeIndexFromId(htmId);
//	cout << "htmId:  " << htmId << ", "<< flush;
//	cout << "nodeId: " << nodeIndex << "; " << flush;
	index->nodeVertex(nodeIndex,v0,v1,v2);
//	cout << "nodeVertex;" << flush;
	if(nodeIndex) {
		addEdge(v0,v1,r,g,b,a,scale);
		addEdge(v1,v2,r,g,b,a,scale);
		addEdge(v2,v0,r,g,b,a,scale);
	}
//	cout << "addedEdges; " << flush;
}

void VizHTM::addFaceFromIndexAndId(
		const SpatialIndex *index, uint64 htmId,
		float r0, float g0, float b0,
		float r1, float g1, float b1,
		float r2, float g2, float b2,
		float a, float scale
		) {
	SpatialVector x1_,x2_,x3_;
	uint64 nodeIndex = index->nodeIndexFromId(htmId);
//	cout << "htmId:  " << htmId << endl << flush;
//	cout << "nodeId: " << nodeIndex << endl << flush;
	index->nodeVertex(nodeIndex,x1_,x2_,x3_);
	if(nodeIndex) {
//		addFace3(
//				v0.x(), v0.y(), v0.z(),
//				v1.x(), v1.y(), v1.z(),
//				v2.x(), v2.y(), v2.z(),
//				r0,g0,b0,
//				r1,g1,b1,
//				r2,g2,b2
//				);

		int colorBase = this->nFaceColors;
		this->addFaceColor(r0,g0,b0,a);
		this->addFaceColor(r1,g1,b1,a);
		this->addFaceColor(r2,g2,b2,a);
		this->addFaceVertexColorIndices3(colorBase,colorBase+1,colorBase+2);

		SpatialVector x1,x2,x3;
		x1 = x1_*scale;
		x2 = x2_*scale;
		x3 = x3_*scale;

		int indexBase = this->nCoordinates;
		this->addCoordinate64(x1.x(),x1.y(),x1.z());
		this->addCoordinate64(x2.x(),x2.y(),x2.z());
		this->addCoordinate64(x3.x(),x3.y(),x3.z());
		this->addFaceIndices3(indexBase,indexBase+1,indexBase+2);

//		addEdge(v0,v1,r,g,b,a);
//		addEdge(v1,v2,r,g,b,a);
//		addEdge(v2,v0,r,g,b,a);
	}
}

void VizHTM::addEdgesFromIndexAndName(
		SpatialIndex *index, const char* htmIdName,
		float r, float g, float b
		) {
	SpatialVector v0, v1, v2;
	uint64 htmId = index->idByName(htmIdName);
	uint64 nodeIndex = index->nodeIndexFromId(htmId);
	index->nodeVertex(nodeIndex,v0,v1,v2);
	addEdge(v0,v1,r,g,b);
	addEdge(v1,v2,r,g,b);
	addEdge(v2,v0,r,g,b);
}

void VizHTM::addEdgesFromIndexAndLatLonDegrees(
		SpatialIndex *index, float64 lat, float64 lon,
		float r, float g, float b
		) {
	SpatialVector v0, v1, v2;
	uint64 htmId = index->idByLatLon(lat,lon);
	uint64 nodeIndex = index->nodeIndexFromId(htmId);
	index->nodeVertex(nodeIndex,v0,v1,v2);
	addEdge(v0,v1,r,g,b);
	addEdge(v1,v2,r,g,b);
	addEdge(v2,v0,r,g,b);
}

void VizHTM::addArcFromIndexAndId(
		const SpatialIndex *index, uint64 htmId,
		float r, float g, float b, float a, float scale
		) {
	SpatialVector v0, v1, v2;
	uint64 nodeIndex = index->nodeIndexFromId(htmId);
	index->nodeVertex(nodeIndex,v0,v1,v2);
	addArc(v0,v1,r,g,b,a,scale);
	addArc(v1,v2,r,g,b,a,scale);
	addArc(v2,v0,r,g,b,a,scale);
}

void VizHTM::addCellFromIndexAndId(
		const SpatialIndex *index, uint64 htmId,
		float r0, float g0, float b0, float a0, float zScale0,
		float r1, float g1, float b1, float a1, float zScale1,
		float hScale
		) {
	SpatialVector v0, v1, v2, vc;
	uint64 nodeIndex = index->nodeIndexFromId(htmId);
	index->nodeVertex(nodeIndex,v0,v1,v2);

	vc = v0 + v1 + v2; vc.normalize();
	v0 = vc * hScale + v0 * ( 1 - hScale ); v0.normalize();
	v1 = vc * hScale + v1 * ( 1 - hScale ); v1.normalize();
	v2 = vc * hScale + v2 * ( 1 - hScale ); v2.normalize();

	addArc(v0,v1,r0,g0,b0,a0,zScale0);
	addArc(v1,v2,r0,g0,b0,a0,zScale0);
	addArc(v2,v0,r0,g0,b0,a0,zScale0);

	addArc(v0,v1,r1,g1,b1,a1,zScale1);
	addArc(v1,v2,r1,g1,b1,a1,zScale1);
	addArc(v2,v0,r1,g1,b1,a1,zScale1);

	addEdge2(
			v0,r0,g0,b0,a0,zScale0,
			v0,r1,g1,b1,a1,zScale1
			);

	addEdge2(
			v1,r0,g0,b0,a0,zScale0,
			v1,r1,g1,b1,a1,zScale1
			);

	addEdge2(
			v2,r0,g0,b0,a0,zScale0,
			v2,r1,g1,b1,a1,zScale1
			);

}

void VizHTM::addArcFromIndexAndName(
		SpatialIndex *index, const char* htmIdName,
		float r, float g, float b, float a
		) {
	SpatialVector v0, v1, v2;
	uint64 htmId = index->idByName(htmIdName);
	uint64 nodeIndex = index->nodeIndexFromId(htmId);
	index->nodeVertex(nodeIndex,v0,v1,v2);
	addArc(v0,v1,r,g,b,a);
	addArc(v1,v2,r,g,b,a);
	addArc(v2,v0,r,g,b,a);
}

void VizHTM::addArcFromIndexAndLatLonDegrees(
		SpatialIndex *index, float64 lat, float64 lon,
		float r, float g, float b, float a
		) {
	SpatialVector v0, v1, v2;
	uint64 htmId = index->idByLatLon(lat,lon);
	uint64 nodeIndex = index->nodeIndexFromId(htmId);
	index->nodeVertex(nodeIndex,v0,v1,v2);
	addArc(v0,v1,r,g,b,a);
	addArc(v1,v2,r,g,b,a);
	addArc(v2,v0,r,g,b,a);
}

void VizHTM::addArcsFromLatLonDegrees(
		float64 *lat, float64 *lon, int nPoints, bool close,
		float r, float g, float b, float a, int steps
		) {
	float scale_ = 1.0;
	SpatialVector v0, v1, v2;
	int i = 0;
	v0.setLatLonDegrees(lat[i],lon[i]);
	v1 = v0;
	for(i=1; i < nPoints; i++) {
		v2.setLatLonDegrees(lat[i],lon[i]);
		addArc(v1,v2,r,g,b,a,scale_,steps);
		v1 = v2;
	}
	if(close) {
		addArc(v2,v0,r,g,b,a,scale_,steps);
	}
}

void VizHTM::addHTMRange(
		const SpatialIndex *index, HtmRange *range,
		float r, float g, float b, float a, float scale,
		bool arcFlag ) {
	if(!range) return;
	range->reset();
	Key lo=0, hi=0;
//	cout << 1100 << "; " << flush;
	int indexp = range->getNext(lo,hi);
	if(indexp) {
		int loLevel = levelOfId(lo);
		int hiLevel = levelOfId(hi);
		if(loLevel!=hiLevel) {
			cout << "addHTMRange error! lo level != hi level" << endl << flush;
			return;
		}
		do {
//			cout << "1101 levels; id: " << flush;
//			cout << loLevel << ", " << hiLevel << "; "
//					<< lo << ", " << hi << "; " << flush;
			for(uint64 numericId=lo; numericId<=hi; numericId++) {
				if(arcFlag){
//					cout << "arc; " << flush;
					addArcFromIndexAndId(index,numericId,r,g,b,a,scale);
				}else{
//					cout << "edge; " << flush;
					addEdgesFromIndexAndId(index,numericId,r,g,b,a,scale);
				}
			}
//			cout << 110 << endl << flush;
		} while (range->getNext(lo,hi));
//		cout << 1103 << flush << endl;
	}
//	cout << 1110 << " addedHTMRange:si " << endl;
}

void VizHTM::addHTMRange(
		HtmRange *range,
		float r, float g, float b, float a, float scale, bool arcFlag) {
	if(!range) return;
	range->reset();
	Key lo=0, hi=0;
	int indexp = range->getNext(lo,hi);
	if(indexp) {
		SpatialIndex index = SpatialIndex(levelOfId(lo));
		do {
			int loLevel = levelOfId(lo);
			int hiLevel = levelOfId(hi);
			if(loLevel!=hiLevel) {
				cout << "addHTMRange error! lo level != hi level" << endl << flush;
				return;
			}
			if(index.getLeafLevel()!=levelOfId(lo)) {
				index.setMaxlevel(levelOfId(lo));
			}
			for(uint numericId=lo; numericId<=hi; numericId++) {
				if(arcFlag) {
					addArcFromIndexAndId(&index,numericId,r,g,b,a,scale);
				}else{
					addEdgesFromIndexAndId(&index,numericId,r,g,b,a,scale);
				}
			}
		} while (range->getNext(lo,hi));
	}
}


void VizHTM::addHstmRange(
		HstmRange *range,
		float r, float g, float b, float a, float scale, bool arcFlag
) {

	int indexLevel = 5, level;
	EmbeddedLevelNameEncoding leftJustified;
	SpatialIndex *index = new SpatialIndex(indexLevel);

	KeyPair kp; int indexp;
	range->reset();
	while((indexp = range->getNext(kp)) > 0) {
		leftJustified.setId(kp.lo);
		level = leftJustified.getLevel();
		string loName = leftJustified.getName();
		uint64 termId = leftJustified.idFromTerminatorAndLevel_NoDepthBit(kp.hi,level);
		leftJustified.setId(termId);
		string hiName = leftJustified.getName();
//		cout << 100 << " lo,hi name: " << loName << " " << hiName << endl << flush;
//		cout << 101 << " level:      " << level << endl << flush;
//		cout << 102 << " kp:         " << hex << kp.lo << " " << kp.hi << dec << endl << flush;
		if( level != indexLevel ) {
			delete index;
			indexLevel = level;
			index = new SpatialIndex(indexLevel);
		}
		uint64 id0 = index->idByName(loName.c_str());
		uint64 id1 = index->idByName(hiName.c_str());
		for(uint64 id=id0; id<=id1; id++) {
			if(arcFlag) {
				addArcFromIndexAndId(
						index,
						id,
						r, g, b,
						a,
						scale
				);
			} else {
				addEdgesFromIndexAndId(
						index,
						id,
						r, g, b,
						a,
						scale
				);
			}
		}
	}
	delete index;
}

void VizHTM::addCellsFromHstmRange(
		HstmRange *range,
		float r0, float g0, float b0, float a0, float scale0,
		float r1, float g1, float b1, float a1, float scale1,
		float hScale
) {

	int indexLevel = 5, level;
	EmbeddedLevelNameEncoding leftJustified;
	SpatialIndex *index = new SpatialIndex(indexLevel);

	KeyPair kp; int indexp;
	range->reset();
	while((indexp = range->getNext(kp)) > 0) {
		leftJustified.setId(kp.lo);
		level = leftJustified.getLevel();
		string loName = leftJustified.getName();
		uint64 termId = leftJustified.idFromTerminatorAndLevel_NoDepthBit(kp.hi,level);
		leftJustified.setId(termId);
		string hiName = leftJustified.getName();
//		cout << 100 << " lo,hi name: " << loName << " " << hiName << endl << flush;
//		cout << 101 << " level:      " << level << endl << flush;
//		cout << 102 << " kp:         " << hex << kp.lo << " " << kp.hi << dec << endl << flush;
		if( level != indexLevel ) {
			delete index;
			indexLevel = level;
			index = new SpatialIndex(indexLevel);
		}
		uint64 id0 = index->idByName(loName.c_str());
		uint64 id1 = index->idByName(hiName.c_str());
		for(uint64 id=id0; id<=id1; id++) {
//			if(arcFlag) {
			addCellFromIndexAndId(
					index,
					id,
					r0, g0, b0, a0, scale0,
					r1, g1, b1, a1, scale1,
					hScale
			);
//			} else {
//				addEdgesFromIndexAndId(
//						index,
//						id,
//						r, g, b,
//						a,
//						scale
//				);
//			}
		}
	}
	delete index;
}






void VizHTM::addHstmRangeIDs(
		HstmRange *range,
		float r, float g, float b, float a, float size, float scale, float dScale, bool arcFlag, bool edgeFlag
) {

	int indexLevel = 5, level;
	EmbeddedLevelNameEncoding leftJustified;
	SpatialIndex *index = new SpatialIndex(indexLevel);

	SpatialVector *v;
	SpatialVector *v0;
	float scale0 = scale;

	KeyPair kp; int indexp;
	range->reset();
	while((indexp = range->getNext(kp)) > 0) {
		leftJustified.setId(kp.lo);
		level = leftJustified.getLevel();
		string loName = leftJustified.getName();
		uint64 termId = leftJustified.idFromTerminatorAndLevel_NoDepthBit(kp.hi,level);
		leftJustified.setId(termId);
		string hiName = leftJustified.getName();
//		cout << 100 << " lo,hi name: " << loName << " " << hiName << endl << flush;
//		cout << 101 << " level:      " << level << endl << flush;
//		cout << 102 << " kp:         " << hex << kp.lo << " " << kp.hi << dec << endl << flush;
		if( level != indexLevel ) {
			delete index;
			indexLevel = level;
			index = new SpatialIndex(indexLevel);
		}

		uint64 id0 = index->idByName(loName.c_str());
		uint64 id1 = index->idByName(hiName.c_str());

		for(uint64 id=id0; id<=id1; id++) {

			SpatialVector vTmp;
			index->pointByHtmId(vTmp,id);
			v = new SpatialVector(vTmp); // Should probably use a shared_ptr or the like.
			v0 = new SpatialVector(vTmp);
			{
				v->normalize(); (*v) *= scale0;
				v0->normalize();

				char *str = new char[256];
				sprintf(str,"0x%llX\n",id);
				cout << "id: " << str;
				cout << " scale0: " << scale0;
				cout << " size: " << size;
				cout << " v: " << (*v);
				cout << " rgb: " << r << " " << g << " " << b;
				cout << endl << flush;
				addAnnotation(v,str,size,r,g,b);

				addEdge(*v0,*v,r,g,b);

			}
			//		index->pointByHtmId(v,id1);
			//		{
			//			SpatialVector *p = new SpatialVector(v); p->normalize(); (*p) *= scale0;
			//			sprintf(str,"%llx\n",id1);
			//			this->addAnnotation(p,str,size,r,g,b);
			//		}
			scale0 += dScale;


			if(arcFlag) {
				addArcFromIndexAndId(
						index,
						id,
						r, g, b,
						a,
						scale
				);
			}
			if(edgeFlag){
				addEdgesFromIndexAndId(
						index,
						id,
						r, g, b,
						a,
						scale
				);
			}
		}
	}
	delete index;
}

void VizHTM::addHstmRangeFaces(
		HstmRange *range,
		float r, float g, float b, float a, float scale
) {

	int indexLevel = 5, level;
	EmbeddedLevelNameEncoding leftJustified;
	SpatialIndex *index = new SpatialIndex(indexLevel);

	KeyPair kp; int indexp;
	range->reset();
	while((indexp = range->getNext(kp)) > 0) {
		leftJustified.setId(kp.lo);
		level = leftJustified.getLevel();
		string loName = leftJustified.getName();
		uint64 termId = leftJustified.idFromTerminatorAndLevel_NoDepthBit(kp.hi,level);
		leftJustified.setId(termId);
		string hiName = leftJustified.getName();
//		cout << 100 << " lo,hi name: " << loName << " " << hiName << endl << flush;
//		cout << 101 << " level:      " << level << endl << flush;
//		cout << 102 << " kp:         " << hex << kp.lo << " " << kp.hi << dec << endl << flush;
		if( level != indexLevel ) {
			delete index;
			indexLevel = level;
			index = new SpatialIndex(indexLevel);
		}
		uint64 id0 = index->idByName(loName.c_str());
		uint64 id1 = index->idByName(hiName.c_str());
		for(uint64 id=id0; id<=id1; id++) {
			addFaceFromIndexAndId(
					index,
					id,
					r, g, b,
					r, g, b,
					r, g, b,
					a, scale
			);
		}
	}
	delete index;
}

void VizHTM::addHTMInterval(SpatialIndex index, htmRange interval) {
	size_t htmIdLevel = index.getMaxlevel();
	Key lo = interval.lo;
	Key hi = interval.hi;
	SpatialVector x1,x2,x3;
	KeyPair adjustedRange = HTMRangeAtLevelFromHTMRange(htmIdLevel,lo,hi);
	lo = adjustedRange.lo;
	hi = adjustedRange.hi;

	for(uint64 numericId=lo; numericId<=hi;numericId++) {
		uint64 nodeIndex = index.nodeIndexFromId(numericId);
		if(nodeIndex!=0){
			index.nodeVertex(nodeIndex,x1,x2,x3);
			if(true) {
				float r=0.; float g=0.; float b=0.;
				switch(numericId % 4) {
				case 0:
					r=1.;
					break;
				case 1:
					g=1.;
					break;
				case 2:
					b=1.;
					break;
				default:
					r=1.; g=1.; b=1.;
					break;
				}
				//	for(int i=0; i<3; i++) addEdgeColor(r,g,b);
				int colorBase = this->nFaceColors;
				for(int i=0; i<3; i++) {
					this->addFaceColor(r,g,b);
				}
				this->addFaceVertexColorIndices3(colorBase,colorBase+1,colorBase+2);

				int indexBase = this->nCoordinates;
				this->addCoordinate64(x1.x(),x1.y(),x1.z());
				this->addCoordinate64(x2.x(),x2.y(),x2.z());
				this->addCoordinate64(x3.x(),x3.y(),x3.z());
				this->addFaceIndices3(indexBase,indexBase+1,indexBase+2);

				//					printf("id: %llx ix: %llu",numericId,nodeIndex);
				//					cout << endl << flush;

				if(true){
					float size = pow(0.5,htmIdLevel+3);
					float r = 0.4, g = 0.4, b = 0.6;
					SpatialVector x = 3.*x1+x2+x3; x.normalize(); x *= 1.0+1.0e-6;
					SpatialVector x_ = x1+x2+x3; x_.normalize();
//					if(false){
//						cout
//						<< " nI: " << nodeIndex
//						<< " x: " << x.x() << " " << x.y() << " " << x.z();
//					}
					char *str = new char[256];
					sprintf(str,"id: %llx\nix: %llu\n",numericId,nodeIndex);
					this->addAnnotation((new SpatialVector(x)),str,size,r,g,b);
					this->addEdge(x,x_,0.5,0.5,0.9);
				}
			}
		}
	}
}

void VizHTM::addHTMInterval(htmRange interval) {
	int loLevel = levelOfId(interval.lo);
	int hiLevel = levelOfId(interval.hi);
	if(loLevel!=hiLevel) {
		cout << "addHTMInterval error! interval.lo level != interval.hi level" << endl << flush;
		return;
	}
	SpatialIndex index = SpatialIndex(loLevel);
	addHTMInterval(index,interval);
}

float64 const PI = atan2(0,-1);
float64 const twoPI = 2.0*PI;

void VizHTM::addCircleFacet(
		SpatialVector center,
		float64 halfSubtendedAngleInRadians,
		float r, float g, float b, float a, float scale) {

	// transparency a is not implemented yet.

	SpatialVector const zhat(0.0,0.0,1.0);
	// Should check to see if center is along zhat.
	SpatialVector c1 = center^zhat; c1.normalize();
	SpatialVector c2 = center^c1;   c2.normalize(); // Now c1^c2 -> center.

	int steps = 20;
	float64 alpha = halfSubtendedAngleInRadians;
	float64 theta = 0.0;
	float64 dTheta = twoPI/steps;

	int iStep = 0;

	SpatialVector x0, x1, x2;

	x0 = center + c1*alpha;
	x1 = x0;

	for( iStep = 0; iStep < steps; iStep ++) {
		theta = iStep*dTheta;
		x2 = center + ( c1 * cos(theta) + c2 * sin(theta) ) * alpha;
		addFace3(
				center*scale,
				x1*scale,
				x2*scale,
				r,g,b,
				r,g,b,
				r,g,b,
				a,a,a
				);
//				r,g,b*iStep/steps,
//				r,g,b*iStep/steps,
//				r,g,b*iStep/steps);
		x1=x2;
	}
	addFace3(
			center*scale,
			x1*scale,
			x0*scale,
			r,g,b,
			r,g,b,
			r,g,b,
			a,a,a
			);
}

void VizHTM::addCircleEdges(
		SpatialVector center,
		float64 halfSubtendedAngleInRadians,
		float r, float g, float b, float a, float scale) {

	// transparency a is not implemented yet.

	SpatialVector const zhat(0.0,0.0,1.0);
	// Should check to see if center is along zhat.
	SpatialVector c1 = center^zhat; c1.normalize();
	SpatialVector c2 = center^c1;   c2.normalize(); // Now c1^c2 -> center.

	int steps = 20;
	float64 alpha = halfSubtendedAngleInRadians;
	float64 theta = 0.0;
	float64 dTheta = twoPI/steps;

	int iStep = 0;

	SpatialVector x0, x1, x2;

	x0 = center + c1*alpha;
	x1 = x0;

	for( iStep = 0; iStep < steps; iStep ++) {
		theta = iStep*dTheta;
		x2 = center + ( c1 * cos(theta) + c2 * sin(theta) ) * alpha;
		addEdge(x1,x2,r,g,b,a,scale);
		x1=x2;
	}
	addEdge(x1,x0,r,g,b,a,scale);

}

void VizHTM::addShapeFile(
		string shapeFile,
		float r, float g, float b,
		bool verbose,
		int nStart, int nEnd) {

	SHPHandle hSHP = SHPOpen(shapeFile.c_str(),"rb");

//	if(hSHP>0) {
	if(hSHP) {
		cout << "file: " << shapeFile << " found." << endl << flush;

		int nShapeType, nEntities;
		double adfMinBound[4], adfMaxBound[4];
		SHPGetInfo(hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );

		cout << "nEntities: " << nEntities << ", nShapeType: " << nShapeType << endl << flush;

		/* Debugging using string shapeFile = "data/ne_50m_coastline/ne_50m_coastline.shp";
		   Asia and Africa
			int nStart = 1387;
			int nEnd   = 1388;
		*/

		if( nStart == -1 ) {
			nStart = 0;
		}
		if( nEnd == -1 ) {
			nEnd = nEntities;
		}
		if( nStart < 0 || nStart > nEntities || nEnd < 0 || nEnd < nStart || nEnd > nEntities ) {
			cout << "loadShapeFile::ERROR::nStart..nEnd" << endl << flush;
			cout << "   nStart..nEnd " << nStart << ".." << nEnd << endl << flush;
			cout << "   nEntities    " << nEntities << endl << flush;
		}

		for(int i=nStart; i < nEnd; i++ ) {
			if(verbose){
				if(i % 100 == 0) {
					cout << "on entity " << i << endl << flush;
				}
			}
			SHPObject* psShape = SHPReadObject(hSHP,i);
			if(verbose){
				cout << "<psShape>" << endl
						<< " nSHPType  " << psShape->nSHPType  << " " << endl
						<< " nVertices " << psShape->nVertices << " " << endl
						<< " nParts    " << psShape->nParts    << " " << endl;
				for(int ipart=0; ipart<psShape->nParts; ipart++) {
					cout << "  part [" << ipart << "] = " << psShape->panPartStart[ipart] << endl;
				}
				cout << "</psShape>" << endl;
			}
			if(psShape->nSHPType == SHPT_ARC) {
				double *yA = psShape->padfY;
				double *xA = psShape->padfX;
				int nVerts = psShape->nVertices;
				int jStart = 0;
				if(psShape->nParts != 0) {
					for(int j=0; j < psShape->nParts; j++) {
						jStart = psShape->panPartStart[j];
						yA = psShape->padfY + psShape->panPartStart[j];
						xA = psShape->padfX + psShape->panPartStart[j];
						if(j < psShape->nParts-1) {
							nVerts = psShape->panPartStart[j+1] - psShape->panPartStart[j];
						} else {
							nVerts = psShape->nVertices - psShape->panPartStart[j];
						}
						if(verbose) cout << "<addArcParts j,start=[" << j << ", " << jStart << "] nVerts= " << nVerts << " />" << endl;
						if( true ) {
//							float r_ = j/(psShape->nParts-1.0);
							float r_ = r;
							addArcsFromLatLonDegrees(
									yA, xA, nVerts,
									false,
									r_,g,b,-1.,3
							);
						}
					}
				} else {
					if(verbose) cout << "<addArc 0,start=[ 0, " << jStart << "] nVerts= " << nVerts << " />" << endl;
					if( true ) {
						addArcsFromLatLonDegrees(
								yA, xA, nVerts,
								false,
								0.,g,b,-1.,3
						);
					}
				}
			}
		}
		if(verbose) cout << endl;
		SHPClose(hSHP);
	}
}
