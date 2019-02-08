/*
 * testSTARE.C
 *
 *  Created on: Feb 7, 2019
 *      Author: mrilee
 *
 *  Copyright (C) 2019 Rilee Systems Technologies LLC
 */

#include "VizHTM.h"

#include "STARE.h"

#include <geompack.h>

#include <unistd.h>
#include <getopt.h>

#include <QtOpenGL/QGL>
#include <QtGui/QImage>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <random>
#include <fstream>
#include <vector>
#include <array>

#include "SpatialException.h"
#include "SpatialIndex.h"
#include "SpatialVector.h"
#include "SpatialInterface.h"

#include "BitShiftNameEncoding.h"
#include "EmbeddedLevelNameEncoding.h"

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
//#include <Inventor/Qt/viewers/SoQtFlyViewer.h>

#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoOrthographicCamera.h>

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

#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoRotationXYZ.h>

void testTenDegreeGrid(VizHTM *viz);
void testShapeFiles(VizHTM *viz);

// For offscreen rendering.
#include "misc.h"
#include "OffScreenViz.h"

using namespace std;

void plotBlockingSphere(VizHTM* viz, float r, float g, float b, float radius);
// { viz->addSphere(SpatialVector(0.,0.,0.),r,g,b,radius); }

SpatialVector *cam_vizCenter = VectorFromLatLonDegrees(33.63,-76.575);
// float offscreenPerscamPositionScale = 1.1;
// float examinerCameraPositionScale   = 1.17;

// string lookFromArgs;
void setLookFrom(SoSeparator *root, SoSeparator *content, SbViewportRegion *vpRegion, string lookFromArgs_ = string());
/*
{

	float scale,lon,lat;

	sscanf(lookFromArgs.c_str(),"%f,%f,%f",&scale,&lon,&lat);

    SoDirectionalLight * light = new SoDirectionalLight;

    SbRotation cameraRotation = SbRotation::identity();

    SoPerspectiveCamera *camera = new SoPerspectiveCamera;
    camera->orientation.setValue(cameraRotation);
    camera->nearDistance = 0.0001;

    //    SoCamera *camera = viewer->getCamera();
    //    SoCamera *camera = perscam;
    SoSFVec3f position;
    //		SpatialVector *p = VectorFromLatLonDegrees(33.0,-80.25); // Centers on gridded data
    //		SpatialVector *p = VectorFromLatLonDegrees(33.5,-76.75); // Centers on storm intersection
    //		SpatialVector *p = VectorFromLatLonDegrees(33.63,-76.575); // Centers on storm intersection
    SpatialVector *p_ = VectorFromLatLonDegrees(lat,lon);
    SpatialVector p = (*p_);
    p = p * scale;

    //    p = p * 1.05; // Works for perscam
//    p = p * 1.1; // Works for perscam
//    p = p * 1.17; // Works for perscam
//    p = p * 1.25; // Nice top level view
//    p = p * 1.5; // Nice top level view
//    cout << "p: " << p << endl << flush;
    position.setValue(p.x(),p.y(),p.z());
    camera->position = position;
//    ((SoPerspectiveCamera*)camera)->position = position;
    //		((SoOrthographicCamera*)camera)->position = position;
    camera->pointAt(SbVec3f(0.,0.,0.));
    //		SoSFRotation rotation;
    //		p->normalize();
    //		rotation.setValue(SbVec3f(p->x(),p->y(),p->z()),0.0);
    //		camera->orientation = rotation;
    //
    root->addChild(light);
    root->addChild(camera); // perscam

//    SoCube * cube = new SoCube;
//    root->addChild(cube);

    root->addChild(content);
    // make sure that the cube is visible
//    perscam->viewAll(root, *vpRegion);
}
*/


void loadScene(SoSeparator *root, SoSeparator *content, SbViewportRegion *vpRegion, SpatialVector *cam_vizCEnter);
/*
{

    SoDirectionalLight * light = new SoDirectionalLight;

    SbRotation cameraRotation = SbRotation::identity();

    SoPerspectiveCamera *camera = new SoPerspectiveCamera;
    camera->orientation.setValue(cameraRotation);
    camera->nearDistance = 0.0001;

    //    SoCamera *camera = viewer->getCamera();
    //    SoCamera *camera = perscam;
    SoSFVec3f position;
    //		SpatialVector *p = VectorFromLatLonDegrees(33.0,-80.25); // Centers on gridded data
    //		SpatialVector *p = VectorFromLatLonDegrees(33.5,-76.75); // Centers on storm intersection
    //		SpatialVector *p = VectorFromLatLonDegrees(33.63,-76.575); // Centers on storm intersection
    SpatialVector p = (*cam_vizCenter);
    p = p * offscreenPerscamPositionScale;
    //    p = p * 1.05; // Works for perscam
//    p = p * 1.1; // Works for perscam
//    p = p * 1.17; // Works for perscam
//    p = p * 1.25; // Nice top level view
//    p = p * 1.5; // Nice top level view
//    cout << "p: " << p << endl << flush;
    position.setValue(p.x(),p.y(),p.z());
    camera->position = position;
//    ((SoPerspectiveCamera*)camera)->position = position;
    //		((SoOrthographicCamera*)camera)->position = position;
    camera->pointAt(SbVec3f(0.,0.,0.));
    //		SoSFRotation rotation;
    //		p->normalize();
    //		rotation.setValue(SbVec3f(p->x(),p->y(),p->z()),0.0);
    //		camera->orientation = rotation;
    //
    root->addChild(light);
    root->addChild(camera); // perscam

//    SoCube * cube = new SoCube;
//    root->addChild(cube);

    root->addChild(content);
    // make sure that the cube is visible
//    perscam->viewAll(root, *vpRegion);
}
*/


int main(int argc, char *argv[]) {

	const char* mainName = "testSTARE";

	cout << "viz..." << flush;
	VizHTM *viz = new VizHTM(NARRAY_);
	cout << "allocated." << endl << flush;

	bool
		examiner_viz           = true,
		blockingSphere_flag    = true,
		testTenDegreeGrid_flag = true,
		testShapeFiles_flag    = true;

	int
		lookFrom_flag = 0;

	string
		lookFromArgs;

	bool
		offscreen_viz = not examiner_viz;

	float
		lineWidth = -1;

	string baseName = "TrmmNmq";


//	testTenDegreeGrid_flag = true;

	if(blockingSphere_flag)    plotBlockingSphere(viz,0.2,0.2,0.2,0.999);
	if(testTenDegreeGrid_flag) testTenDegreeGrid(viz);
	if(testShapeFiles_flag)    testShapeFiles(viz);

	viz->addLatLonBoxEdgesDegrees(0,0,5,5,0.9,0.9,0.9);

	LatLonDegrees64ValueVector latlonbox;
	latlonbox.push_back(LatLonDegrees64(0,0));
	latlonbox.push_back(LatLonDegrees64(5,0));
	latlonbox.push_back(LatLonDegrees64(5,5));
	latlonbox.push_back(LatLonDegrees64(0,5));

	STARE index;
	// STARE_Intervals intervals = index.BoundingBoxFromLatLonDegrees(latlonbox,6);
	STARE_Intervals intervals = index.BoundingBoxFromLatLonDegrees(latlonbox);
	EmbeddedLevelNameEncoding leftJustified;
	BitShiftNameEncoding      rightJustified;

	/*
	float
	color_scale = 1.0/intervals.size(),
	r0 = 0, g0 = 1, b0 = 0,
	r1 = 0, g1 = 1, b1 = 0,
	r2 = 0, g2 = 1, b2 = 0,
	a0 = 0, a1 = 0, a2 = 0,
	scale = 1;
	*/

	float
		color_scale = 1.0/intervals.size(),
		r0 = 1, g0 = 0, b0 = 0,
		r1 = 0, g1 = 1, b1 = 0,
		r2 = 0, g2 = 0, b2 = 1,
		a0 = 0, a1 = 0, a2 = 0,
		scale = 1;

	for( STARE_Intervals::iterator iSid = intervals.begin(); iSid != intervals.end(); ++iSid ) {
		Triangle tr0 = index.TriangleFromValue(*iSid);

		// leftJustified.setIdFromSciDBLeftJustifiedFormat(*iSid);
		// rightJustified.setId(leftJustified.rightJustifiedId());
		// id_line = rightJustified.getId();

		// cout << "centroid: " << tr0.centroid    << endl;
		// cout << "tr0.0:    " << tr0.vertices[0] << endl;
		// cout << "tr0.1:    " << tr0.vertices[1] << endl;
		// cout << "tr0.2:    " << tr0.vertices[2] << endl;

		viz->addFace3(tr0.vertices[0], tr0.vertices[1], tr0.vertices[2], r0, g0, b0, r1, g1, b1, r2, g2, b2, a0, a1, a2, scale);

		/*
		r0 += color_scale;
		r1 += color_scale;
		r2 += color_scale;

		g0 -= color_scale;
		g1 -= color_scale;
		g2 -= color_scale;

		b0 += color_scale;
		b1 += color_scale;
		b2 += color_scale;
		*/
	}

	// Last chance change...
	if(lineWidth != -1) {
		viz->lineWidth = lineWidth;
	}

	std::vector<std::shared_ptr<VizHTM>> vizContainer;

	if(false) { // Test vizContainer
		std::shared_ptr<VizHTM> v0 = std::make_shared<VizHTM>(10000);
		vizContainer.push_back(v0);
		v0->triaxis();

		std::shared_ptr<VizHTM> v1 = std::make_shared<VizHTM>(10000);
		vizContainer.push_back(v1);
		v1->addSphere(SpatialVector(0.0,0.0,1.1),0.,1.,0.,0.1);
	}

	/**************** Start Graphics ****************/

	const int width = 2000, height = 1400;


	QWidget *window = SoQt::init(argv[0]);
	window->setMinimumSize(width,height);
	if (window == NULL) exit(1);

	SoSelection *selectionRoot = new SoSelection;
	selectionRoot->policy = SoSelection::SINGLE;

	SoSeparator *root = new SoSeparator;
	//	root->ref(); // TODO Figure out ->ref();

	root->addChild(viz->makeRoot());

	selectionRoot->addChild(root);

//	cout << 10000 << endl << flush;

	SoSeparator *roots = new SoSeparator;
	for(std::vector<std::shared_ptr<VizHTM>>::iterator it = vizContainer.begin();
			it != vizContainer.end();
			++it ) {
		roots->addChild((*it)->makeRoot());
	}
	selectionRoot->addChild(roots);

//	cout << 10100 << endl << flush;

	// Offscreen renderer for viz. TODO selectionRoot

	if(offscreen_viz){
//		OffScreenViz *offscreen = new OffScreenViz(800,600);
		OffScreenViz *offscreen = new OffScreenViz(width,height);
		offscreen->initImageDirectory("tmp/offscreen/"+formattedDateTime()+"/"+baseName+"/",4);
		offscreen->root = new SoSeparator;
		if(!lookFrom_flag) {
			loadScene(offscreen->root,selectionRoot,offscreen->vpRegion,cam_vizCenter);
		} else {
			setLookFrom(offscreen->root,selectionRoot,offscreen->vpRegion);
		}
		offscreen->saveImage(1);
	}

//	cout << 10200 << endl << flush;

	// Examiner viewer
	if(examiner_viz){
		SoQtExaminerViewer *viewer = new SoQtExaminerViewer(window);

		//	SoQtFlyViewer *viewer = new SoQtFlyViewer(window); // Still fails with transparency sorted_*
		// Transparency on examiner viewer is broken on my Mac Pro.
		//	viewer->setTransparencyType(SoGLRenderAction::DELAYED_ADD); // Crash
		//	viewer->setTransparencyType(SoGLRenderAction::NONE); // No crash
		//	viewer->setTransparencyType(SoGLRenderAction::DELAYED_BLEND); // Crash
		//	viewer->setTransparencyType(SoGLRenderAction::ADD); // Weird
		//	viewer->setTransparencyType(SoGLRenderAction::BLEND); // Weird. No transparency!
		//	viewer->setTransparencyType(SoGLRenderAction::SORTED_OBJECT_BLEND);
		// viewer->setTransparencyType(SoGLRenderAction::BLEND);
		viewer->setTransparencyType(SoGLRenderAction::NONE);

		viewer->setSceneGraph(selectionRoot);
		viewer->setTitle(mainName);
		viewer->show();

		if(lookFrom_flag) {
			SoCamera *camera = viewer->getCamera();
			if(camera) {

				float scale,lon,lat;
				cout << "lookFrom scanning " << lookFromArgs.c_str() << endl;
				sscanf(lookFromArgs.c_str(),"%f,%f,%f",&scale,&lon,&lat);
				cout << "Setting interactive lookFrom: (scale,lon,lat) = "
						<< scale << ", " << lon << ", " << lat << endl;
			    SpatialVector *p_ = VectorFromLatLonDegrees(lat,lon);
			    SpatialVector p = (*p_);
			    p = p * scale;

				SoSFVec3f position;
				position.setValue(p.x(),p.y(),p.z());
				//		((SoPerspectiveCamera*)camera)->position = position;
				camera->position = position;
				camera->orientation.setValue(SbRotation::identity());
				camera->pointAt(SbVec3f(0.,0.,0.));

				//		SoSFRotation rotation;
				//		p->normalize();
				//		rotation.setValue(SbVec3f(p->x(),p->y(),p->z()),0.0);
				//		camera->orientation = rotation;
			}
		}

		if(false) {
			SoCamera *camera = viewer->getCamera();
			if(camera) {
				SoSFVec3f position;
				//		SpatialVector *p = VectorFromLatLonDegrees(33.0,-80.25); // Centers on gridded data
				//		SpatialVector *p = VectorFromLatLonDegrees(33.5,-76.75); // Centers on storm intersection
				//		SpatialVector *p = VectorFromLatLonDegrees(33.63,-76.575); // Centers on storm intersection
				SpatialVector *p = cam_vizCenter;
				float examinerCameraPositionScale   = 1.17;
				(*p) = (*p) * examinerCameraPositionScale;
				position.setValue(p->x(),p->y(),p->z());
				//		((SoPerspectiveCamera*)camera)->position = position;
				camera->position = position;
				camera->pointAt(SbVec3f(0.,0.,0.));
				//		SoSFRotation rotation;
				//		p->normalize();
				//		rotation.setValue(SbVec3f(p->x(),p->y(),p->z()),0.0);
				//		camera->orientation = rotation;
			}
		}

		SoQt::show(window);
		SoQt::mainLoop();

		delete viewer;
		//	root->unref();
	}

	cout << mainName << " done." << endl;

	return 0;
}


