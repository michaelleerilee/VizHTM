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

	faceColors 	            = new float[nArray][3];
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
//	cout << "  (adding-coord64 " << x << " " << y << " " << z << ") " << endl << flush;
	fCoordinates[nCoordinates][0] = (float) x;
	fCoordinates[nCoordinates][1] = (float) y;
	fCoordinates[nCoordinates][2] = (float) z;
	nCoordinates++;
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

	SoSeparator *root = new SoSeparator;

	if(lineWidth != -1) {
		SoDrawStyle *style = new SoDrawStyle;
		style->lineWidth.setValue(lineWidth);
		root->addChild(style);
	}

	SoCoordinate3 *coordinates = new SoCoordinate3;
	coordinates->point.setValues(0,nCoordinates,fCoordinates);
	root->addChild(coordinates);

	SoSeparator *edgeNode = new SoSeparator;

	SoMaterial *edgeMaterials = new SoMaterial;
	edgeMaterials->diffuseColor.setValues(0,nEdgeColors,edgeColors);
	if(edgeTransparency) {
		edgeMaterials->transparency.setValues(0,nEdgeColors,edgeTransparencies);
	}
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

//	{
//		cout << "Starting text test" << flush;
//		SpatialVector a = SpatialVector(0.0,0.0,1.0);
//		float size = 1.0;
//		float r = 0.8, g = 0.8, b = 0.9;
//		addText(a,"Test",size,r,g,b);
//		cout << "...done." << endl << flush;
//	}

	if(nAnnotations>0){
		SoSeparator *texts = new SoSeparator;
		for(int ia=0;ia<nAnnotations;ia++) {
//			cout
//				<< " adding annotation: " << ia
//				<< " text: " << annotations[ia].text
//				<< endl << flush;
			texts->addChild(makeText(
					annotations[ia].v,
					annotations[ia].text,
					annotations[ia].size,
					annotations[ia].r,
					annotations[ia].g,
					annotations[ia].b));
		}
		root->addChild(texts);
	}

	root->unrefNoDelete();
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
		float r, float g, float b, float a) {
	int colorBase = nEdgeColors;
	addEdgeColor(r,g,b,a);
	addEdgeColor(r,g,b,a);
	addEdgeVertexColorIndices(colorBase,colorBase+1);

	int coordBase = nCoordinates;
	addCoordinate(x0);
	addCoordinate(x1);
	addEdgeIndices(coordBase,coordBase+1);
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
	SoSeparator *root = new SoSeparator;

	SoFont *font = new SoFont;
	font->name.setValue("Times-Roman");
	font->size.setValue(size);
	root->addChild(font);

	SoMaterial *material = new SoMaterial;
	SoMaterialBinding *binding = new SoMaterialBinding;
	material->diffuseColor.set1Value(0,SbColor(r,g,b));
	material->diffuseColor.set1Value(1,SbColor(0.5*r,0.5*g,0.5*b));
	binding->value = SoMaterialBinding::PER_PART;
	root->addChild(material);
	root->addChild(binding);

	SoSeparator *textSeparator = new SoSeparator;
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

	text->parts = SoText3::FRONT;
//	text->parts = SoText3::ALL;

	string *s = new string(annotation);
	string delimiter = "\n";
	int idx=0; int last=0; int next=0; while((next=s->find(delimiter,last)) != string::npos) {
		text->string.set1Value(idx++,(s->substr(last,next-last)).c_str());
		last = next + 1;
	}

	textSeparator->addChild(rot0);
	textSeparator->addChild(tra0);
	textSeparator->addChild(text);

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
		float r, float g, float b) {
	addArc(x0,x1,r,g,b);
	addArc(x1,x2,r,g,b);
	addArc(x2,x3,r,g,b);
	addArc(x3,x0,r,g,b);
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
		float r, float g, float b, float alpha,
		int steps) {
	SpatialVector dx = x1-x0;
	SpatialVector *last = new SpatialVector(x0);
	if(steps<0) {
		steps = 20;
		if(dx.length()>0.5){
			steps = 180;
		}
	}
	for(double i=1; i<steps-1; i++) {
		double a = i / (steps - 1);
		SpatialVector *next = new SpatialVector(x0 + a * dx); next->normalize();
		addEdge(*last,*next, r, g, b, alpha);
		last = next;
	}
	addEdge(*last,x1,r,g,b,alpha);
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
		float r, float g, float b, float a
		) {
	SpatialVector v0, v1, v2;
	uint64 nodeIndex = index->nodeIndexFromId(htmId);
//	cout << "htmId:  " << htmId << endl << flush;
//	cout << "nodeId: " << nodeIndex << endl << flush;
	index->nodeVertex(nodeIndex,v0,v1,v2);
	if(nodeIndex) {
		addEdge(v0,v1,r,g,b,a);
		addEdge(v1,v2,r,g,b,a);
		addEdge(v2,v0,r,g,b,a);
	}
}

void VizHTM::addFaceFromIndexAndId(
		const SpatialIndex *index, uint64 htmId,
		float r0, float g0, float b0,
		float r1, float g1, float b1,
		float r2, float g2, float b2
		) {
	SpatialVector x1,x2,x3;
	uint64 nodeIndex = index->nodeIndexFromId(htmId);
//	cout << "htmId:  " << htmId << endl << flush;
//	cout << "nodeId: " << nodeIndex << endl << flush;
	index->nodeVertex(nodeIndex,x1,x2,x3);
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
		this->addFaceColor(r0,g0,b0);
		this->addFaceColor(r1,g1,b1);
		this->addFaceColor(r2,g2,b2);
		this->addFaceVertexColorIndices3(colorBase,colorBase+1,colorBase+2);

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
		SpatialIndex *index, uint64 htmId,
		float r, float g, float b, float a
		) {
	SpatialVector v0, v1, v2;
	uint64 nodeIndex = index->nodeIndexFromId(htmId);
	index->nodeVertex(nodeIndex,v0,v1,v2);
	addArc(v0,v1,r,g,b,a);
	addArc(v1,v2,r,g,b,a);
	addArc(v2,v0,r,g,b,a);
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

void VizHTM::addHTMRange(
		const SpatialIndex *index, HtmRange *range,
		float r, float g, float b, float a) {
	if(!range) return;
	range->reset();
	Key lo=0, hi=0;
	int indexp = range->getNext(lo,hi);
	if(indexp) {
		int loLevel = levelOfId(lo);
		int hiLevel = levelOfId(hi);
		if(loLevel!=hiLevel) {
			cout << "addHTMRange error! lo level != hi level" << endl << flush;
			return;
		}
		do {
			for(uint numericId=lo; numericId<=hi; numericId++) {
				addEdgesFromIndexAndId(index,numericId,r,g,b,a);
			}
		} while (range->getNext(lo,hi));
	}
}

void VizHTM::addHTMRange(
		HtmRange *range,
		float r, float g, float b, float a) {
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
				addEdgesFromIndexAndId(&index,numericId,r,g,b,a);
			}
		} while (range->getNext(lo,hi));
	}
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
