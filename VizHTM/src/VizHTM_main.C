/*
 * VizHTM_main.C
 *
 *  Created on: Dec 16, 2015
 *      Author: mrilee
 */




#include "VizHTM.h"

#include <iostream>
#include <iomanip>
#include <random>

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

int main(int argc, char *argv[]) {
	const char* mainName = "VizHTM-main";
	cout << mainName << endl;

	float64 PI = atan2(0,-1); cout << "PI=" << PI << endl << flush;

	VizHTM *viz = new VizHTM(NARRAY_);


	if(true) {
		viz->triaxis(); // TODO Breaks if first?
	}

	if(false) {
		viz->addFace3(
				1.0,0.0,0.0,
				0.0,1.0,0.0,
				0.0,0.0,1.0,
				1.0,0.0,0.0,
				0.0,1.0,0.0,
				0.0,0.0,1.0
		);
	}

	SpatialVector o = SpatialVector(0.0,0.0,0.0);

//	SpatialVector a = SpatialVector(1.0,0.0,0.0);
//	SpatialVector as = SpatialVector(1.05,0.0,0.0);

	SpatialVector a = randomVector();
	SpatialVector as = 1.05*a;

//	int colorBase = viz->nEdgeColors;
//	viz->addEdgeColor(0,1,1);
//	viz->addEdgeColor(0,1,1);
//	viz->addEdgeVertexColorIndices(colorBase,colorBase+1);
//
//	int coordBase = viz->nCoordinates;
//	viz->addEdgeIndices(coordBase,coordBase+1);
//
//	viz->addCoordinate(o);
//	viz->addCoordinate(a);

	int coordBase = viz->nCoordinates;
	viz->addCoordinate(as);
	viz->addSphere(coordBase,0,1,1,0.025);

	viz->addEdge(o,a,1.,1.,1.);

	float64 d = 0.3;
	viz->addConstraint(a,d,0.5,0.5,1.0);
	viz->addEdgeProjections(a*(1.-d));

	for(float64 d = 0.; d < 2.0; d += 0.2) {

//		viz->addConstraint(a,d,0.5,0.5,1.0);
//
//		SpatialVector ad = (1.0-d)*a;
//
//		SpatialVector i, a_cross_i;
//		float64 aci_norm;
//
//		do {
//			i = randomVector();
//			cout << "i: " << i << endl << flush;
//			a_cross_i = a^i;
//			aci_norm = a_cross_i.length();
//		} while (!aci_norm);
//
//		SpatialVector b = a_cross_i; b.normalize();
//		SpatialVector c = a^b; c.normalize();
//		viz->addEdgeAndSphere(o,b,0.,1.,0.5*d,b,0.,1.,0.5*d,0.025);
//		viz->addEdgeAndSphere(o,c,0.,0.5*d,1.,c,0.,0.5*d,1.,0.025);
//
//		//	for(int i=0; i<100; i++) cout << "i: " << uniformDouble() << endl << flush;
//
//		float64 theta = 0.;
//		float64 twopi  = atan2(0.,-0.);
//		float64 deltaTheta = 0.05;
//		SpatialVector dlast = sin(-twopi*(1.+deltaTheta))*b + cos(-twopi*(1.+deltaTheta))*c; dlast.normalize();
//		dlast *= sqrt(d*(2.-d));
//		SpatialVector last = ad + dlast;
//		for(float64 theta = -1.; theta<=1.; theta += deltaTheta) {
//			float64 lambda = sin(twopi*theta);
//			float64 mu     = cos(twopi*theta);
//			SpatialVector delta_ad = lambda*b + mu*c; delta_ad.normalize();
//			delta_ad *= sqrt(d*(2.-d));
//			SpatialVector aPlusDelta = ad + delta_ad;
////			viz->addEdge(ad,aPlusDelta,1.,0.5*d,1.);
////			viz->addEdge(o,aPlusDelta,1.,0.5*d,1.);
//			viz->addEdge(last,aPlusDelta,1.0,0.,0.);
//			last = aPlusDelta;
//		}
	}

	if(false) {
		viz->triaxis();
	}

	if(false){ // Try to draw latlon lines.
		int indexBase = viz->nCoordinates;
		int edgeColorBase = viz->nEdgeColors;
		int nPoints = 1001;
//		nPoints = 3;
		for(int i = 0;i < nPoints; i++) {
			float lat = -90.0 + (i*180.)/(nPoints-1);	float lon = -(i*16*4*90.)/(nPoints-1);
			float* xyz = xyzFromLatLonDegrees(lat,lon);
			//		cout << "1000: " << i << " ( " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " ) " << endl;
			viz->addCoordinate(xyz[0],xyz[1],xyz[2]);
		}
		for(int i=0;i<nPoints;i++) {
			edgeColorBase = viz->nEdgeColors;
			float r = (i/(1.0*nPoints));
			float g = 1.0;
			float b = 1.0;
			viz->addEdgeColor(r,g,b);
			viz->addEdgeVertexColorIndices(edgeColorBase,edgeColorBase);
		}
		int nEdges = nPoints-1;
		for(int iPoint = 0; iPoint < nEdges; iPoint++) {
			viz->edgeIndices[viz->nEdgeIndices] = indexBase + iPoint; viz->nEdgeIndices++;
			viz->edgeIndices[viz->nEdgeIndices] = indexBase + iPoint + 1; viz->nEdgeIndices++;
			viz->edgeIndices[viz->nEdgeIndices] = SO_END_LINE_INDEX; viz->nEdgeIndices++;
		}
//		cout << "1000: " << nCoordinates << " " << nEdgeColors << " " << nFaceColors << " " << nEdgeIndices << " " << nFaceIndices << endl << endl;
	}

	if(false) {
		viz->triaxis();
	}

	if(false) viz->debug_dump();

	QWidget *window = SoQt::init(argv[0]);
	if (window == NULL) exit(1);

	SoSelection *selectionRoot = new SoSelection;
	selectionRoot->policy = SoSelection::SINGLE;

	SoSeparator *root = new SoSeparator;
//	root->ref();

	root->addChild(viz->makeRoot());

	selectionRoot->addChild(root);

	SoQtExaminerViewer *viewer = new SoQtExaminerViewer(window);
	viewer->setSceneGraph(selectionRoot);
	viewer->setTitle(mainName);
	viewer->show();

	SoQt::show(window);
	SoQt::mainLoop();

	delete viewer;
//	root->unref();

	return 0;
}

