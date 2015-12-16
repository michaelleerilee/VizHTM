/*
 * VizHTM_main.C
 *
 *  Created on: Dec 16, 2015
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

int main(int argc, char *argv[]) {
	const char* mainName = "VizHTM-main";
	cout << mainName << endl;

	VizHTM *viz = new VizHTM(NARRAY_);


	if(true) {
		viz->triaxis(); // TODO Breaks if first?
	}

	if(true) {
		viz->addFace3(
				1.0,0.0,0.0,
				0.0,1.0,0.0,
				0.0,0.0,1.0,
				1.0,0.0,0.0,
				0.0,1.0,0.0,
				0.0,0.0,1.0
		);
	}

	if(false) {
		viz->triaxis();
	}

	if(true){ // Try to draw latlon lines.
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

