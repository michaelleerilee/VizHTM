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

#include "tests/tests.h"


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
// #include <Inventor/Qt/Ui.h>

#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoOrthographicCamera.h>

#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoFont.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoIndexedLineSet.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoMaterialBinding.h>
#include <Inventor/nodes/SoSelection.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoSwitch.h>
#include <Inventor/nodes/SoText3.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/nodes/SoTransparencyType.h>
#include <Inventor/sensors/SoTimerSensor.h>
#include <Inventor/SbBasic.h>

#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoRotationXYZ.h>

#include <Inventor/SbVec2f.h>

#include <Inventor/Qt/SoQtGLWidget.h>

void testTenDegreeGrid(VizHTM *viz,float r0, float g0, float b0, float rgbScale);
void testShapeFiles(VizHTM *viz, float r, float g, float b, float deltaZ);

// For offscreen rendering.
#include "misc.h"
#include "OffScreenViz.h"

using namespace std;


SoGroup* testText() {

	  SoGroup *root = new SoGroup;
	  root->ref();

	  // Choose a font
	  SoFont *myFont = new SoFont();

	  // The Stroke Font !
	  myFont->name.setValue("TGS_Complex_Roman");
	  myFont->size.setValue(0.2f);
	  root->addChild(myFont);
	  SoBaseColor *myBaseColor = new SoBaseColor();
	  myBaseColor->rgb.setValue(SbColor(1.0f, 0.0f, 0.0f));
	  root->addChild(myBaseColor);

	  // Add Text3D with stroke font.
	  SoText3 *strokeFontText = new SoText3();
	  strokeFontText->string = "STROKE FONT !";
	  root->addChild(strokeFontText);

	  return root;
}

void plotBlockingSphere(VizHTM* viz, float r, float g, float b, float radius);
// { viz->addSphere(SpatialVector(0.,0.,0.),r,g,b,radius); }

// SpatialVector *cam_vizCenter = VectorFromLatLonDegrees(33.63,-76.575);
// SpatialVector *cam_vizCenter = VectorFromLatLonDegrees(90.0, 0.0);
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

/*
void loadScene(SoSeparator *root, SoSeparator *content, SbViewportRegion *vpRegion, SpatialVector *cam_vizCEnter);
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
void loadScene(SoSeparator *root, SoSeparator *content, SbViewportRegion *vpRegion, SpatialVector *cam_vizCenter, float offscreenPerscamPositionScale, SoCamera *inputCamera = NULL ) {

    SoDirectionalLight * light = new SoDirectionalLight;
    root->addChild(light);

    SoPerspectiveCamera *camera;
    if( !inputCamera ) {
    	SbRotation cameraRotation = SbRotation::identity();
    	camera = new SoPerspectiveCamera;
    	camera->orientation.setValue(cameraRotation);
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
        camera->nearDistance = 0.0001;
        root->addChild(camera); // perscam
    } else {
        root->addChild(inputCamera); // perscam
    }
//    SoCube * cube = new SoCube;
//    root->addChild(cube);

    root->addChild(content);
    // make sure that the cube is visible
//    perscam->viewAll(root, *vpRegion);
}

void loadTestScene
	(SoSeparator *root, SbViewportRegion *vpRegion) {
    // Init Coin
    // SoDB::init();
    // The root node
    // root = new SoSeparator;
    // root->ref();

	// cout << 100 << endl << flush;
    // Add the light _before_ you add the camera
    SoDirectionalLight * light = new SoDirectionalLight;
    // cout << 110 << endl << flush;
    root->addChild(light);

    // cout << 200 << endl << flush;
    vpRegion->setViewportPixels(0, 0, 450, 450);

    // cout << 300 << endl;
    SoPerspectiveCamera *perscam = new SoPerspectiveCamera();
    root->addChild(perscam);

    // cout << 400 << endl;
    SbRotation cameraRotation = SbRotation::identity();
    cameraRotation *= SbRotation(SbVec3f(0, 1, 0), 0.4f);
    perscam->orientation = cameraRotation;

    // cout << 500 << endl;
    SoCube * cube = new SoCube;
    root->addChild(cube);

    // cout << 600 << endl;
    // make sure that the cube is visible
    perscam->viewAll(root, *vpRegion);

}

int main(int argc, char *argv[]) {

	bool ok = false;

	QWidget *window = SoQt::init(argv[0]);

	const char* mainName = "testSTARE";

	cout << "viz..." << flush;
	VizHTM *viz = new VizHTM(NARRAY_);
	cout << "allocated." << endl << flush;

	// Project... Needs to come first.
	// cout << 100 << endl << flush;
	ok = viz->setProjection("Equirectangular");
	// ok = viz->setProjection("Mercator");
	// cout << 200 << endl << flush;

	bool
		examiner_viz           = false,
		blockingSphere_flag    = false,
		testTenDegreeGrid_flag = false,
		testShapeFiles_flag    = false;

	int
		lookFrom_flag = 0;

	string
		lookFromArgs;

	SpatialVector
		lookFrom;

	bool
		offscreen_viz = not examiner_viz;

	float
		lineWidth = -1;

	// string baseName = "TrmmNmq";
	string baseName = "testSTARE";

//	testTenDegreeGrid_flag = true;
	float
	r0       = 0.5,
	g0       = 0.5,
	b0       = 0.5,
	rgbScale = 0
	;

	// Global diagnostics
	if(false) {
		blockingSphere_flag    = true;
		testTenDegreeGrid_flag = true;
		testShapeFiles_flag    = true;
	}

	// blockingSphere_flag = false;
	// testShapeFiles_flag = false;

	if(blockingSphere_flag)    plotBlockingSphere(viz,0.2,0.2,0.2,0.999);
	if(testTenDegreeGrid_flag) testTenDegreeGrid(viz,r0,g0,b0,rgbScale);
	if(testShapeFiles_flag)    testShapeFiles(viz,0.5,1,1,0.0);

	// ok = BoundingBox1(viz);
	// ok = Edges1(viz);
	// ok = Edges2(viz); // Looking at a case were data lies on an edge.

	// ok = EquatorCheck1(viz); lookFromArgs = "0.0 1.1 0.0"; lookFrom_flag = true;
	// ok = EquatorCheck1(viz); lookFromArgs = "1.00001 0.0 0.0"; lookFrom_flag = true;
	// ok = EquatorCheck2(viz);

	/*
	ok = MeridianCheck1(viz);
	lookFromArgs = "";
	lookFrom.setLatLonDegrees(45.0, 0.0);
	lookFrom_flag = false;
	*/

	// ok = PolePosition1(viz);
	// ok = PoleCheck1(viz);

	// ok = Granule1(viz);
	ok = Chunk1(viz);

	// Last chance changes...
	if(lineWidth != -1) {
		viz->lineWidth = lineWidth;
	}

	/**************** Viz Containers! ****************/

	std::vector<std::shared_ptr<VizHTM>> vizContainer;
	bool enableVizContainers = false;

	if(enableVizContainers) { // Test vizContainer
		std::shared_ptr<VizHTM> v0 = std::make_shared<VizHTM>(10000);
		vizContainer.push_back(v0);
		v0->triaxis();

		std::shared_ptr<VizHTM> v1 = std::make_shared<VizHTM>(10000);
		vizContainer.push_back(v1);
		v1->addSphere(SpatialVector(0.0,0.0,1.1),0.,1.,0.,0.1);
	}

	/**************** Start Graphics ****************/

	// const int width = 2000, height = 1400;
	// const int width = 1600, height = 800;
	// const int width = 1000, height = 500;
	// const int width = 800, height = 600;
	const int width = 800, height = 400;

	window->setMinimumSize(width,height);
	if (window == NULL) exit(1);

	SoSelection *selectionRoot = new SoSelection;
	selectionRoot->policy = SoSelection::SINGLE;
	selectionRoot->ref();

	SoSeparator *root = new SoSeparator;
	//	root->ref(); // TODO Figure out ->ref();

	root->addChild(viz->makeRoot());
	// root->addChild(testText());

	if( true ) {
		selectionRoot->addChild(root);
	}

//	cout << 10000 << endl << flush;

	if( enableVizContainers ){
		SoSeparator *roots = new SoSeparator;
		for(std::vector<std::shared_ptr<VizHTM>>::iterator it = vizContainer.begin();
				it != vizContainer.end();
				++it ) {
			roots->addChild((*it)->makeRoot());
		}
		selectionRoot->addChild(roots);
	}

//	cout << 10100 << endl << flush;

	// Offscreen renderer for viz. TODO selectionRoot

	SpatialVector *cam_vizCenter = VectorFromLatLonDegrees(90.0, 0.0);
	float offscreenPerscamPositionScale = 3;

//	offscreen_viz = false;
//	if(offscreen_viz){
//		cout << "Offscreen rendering and saving to a file." << endl << flush;
//		cout << "width,height: " << width << ", " << height << endl << flush;
////		OffScreenViz *offscreen = new OffScreenViz(800,600);
//		OffScreenViz *offscreen = new OffScreenViz(width,height);
//		// offscreen->initImageDirectory("tmp/offscreen/"+formattedDateTime()+"/"+baseName+"/",4);
//		offscreen->initImageDirectory("/home/mrilee/workspace/VizHTM/tmp/offscreen/"+formattedDateTime()+"/"+baseName+"/",4);
//		offscreen->root = new SoSeparator;
//		offscreen->root->ref();
//		if( true ) {
//			if(!lookFrom_flag) {
//				cout << "Calling loadScene" << endl << flush;
//				loadScene(offscreen->root,selectionRoot,offscreen->vpRegion,cam_vizCenter,offscreenPerscamPositionScale);
//			} else {
//				cout << "Calling setLookFrom" << endl << flush;
//				setLookFrom(offscreen->root,selectionRoot,offscreen->vpRegion);
//			}
//		} else {
//			cout << "Loading test scene" << endl << flush;
//			loadTestScene(offscreen->root,offscreen->vpRegion);
//			cout << "Adding test scene to selectionRoot" << endl << flush;
//			selectionRoot->addChild(offscreen->root);
//		}
//		// offscreen->writeFileRGB(1);
//		offscreen->writeFile(1);
//		// offscreen->saveImage(1);
//
//	}
//	cout << 10200 << endl << flush;

	cout << "Checking for examiner visualization" << endl << flush;

	examiner_viz = true;
	// Examiner viewer
	if(examiner_viz){
		SoQtExaminerViewer *viewer = new SoQtExaminerViewer(window);

		{
			SbVec2f r_; float granularity;
			viewer->getPointSizeLimits(r_,granularity);
			cout << "pointsize range: " << r_[0] << "-" << r_[1] << ", granularity = " << granularity << endl << flush;
			viewer->getLineWidthLimits(r_,granularity);
			cout << "pointsize range: " << r_[0] << "-" << r_[1] << ", granularity = " << granularity << endl << flush;
		}

		//	SoQtFlyViewer *viewer = new SoQtFlyViewer(window); // Still fails with transparency sorted_*
		// Transparency on examiner viewer is broken on my Mac Pro.
		//	viewer->setTransparencyType(SoGLRenderAction::DELAYED_ADD); // Crash
		//	viewer->setTransparencyType(SoGLRenderAction::NONE); // No crash
		// viewer->setTransparencyType(SoGLRenderAction::DELAYED_BLEND); // Crash
		//	viewer->setTransparencyType(SoGLRenderAction::ADD); // Weird
		//	viewer->setTransparencyType(SoGLRenderAction::BLEND); // Weird. No transparency!
		// viewer->setTransparencyType(SoGLRenderAction::SORTED_OBJECT_BLEND);
		// viewer->setTransparencyType(SoGLRenderAction::BLEND);
		viewer->setTransparencyType(SoGLRenderAction::NONE);

		viewer->setSceneGraph(selectionRoot);
		viewer->setTitle(mainName);
		viewer->show();
		viewer->viewAll();

		if(lookFrom_flag) {
			SoCamera *camera = viewer->getCamera();
			if(camera) {

				float scale=1.00001,lon=0.0,lat=0.0;
				SpatialVector p;
				if(lookFromArgs != "") {
					cout << "lookFrom scanning " << lookFromArgs.c_str() << endl;
					sscanf(lookFromArgs.c_str(),"%f,%f,%f",&scale,&lon,&lat);
					cout << "Setting interactive lookFrom: (scale,lon,lat) = "
							<< scale << ", " << lon << ", " << lat << endl;
					SpatialVector *p_ = VectorFromLatLonDegrees(lat,lon);
					p = (*p_);
				} else {
					cout << "Using lookFrom = " << lookFrom << endl << flush;
					p = lookFrom;
				}
				p = p * scale;
				cout << "Setting camera to position p = " << p << endl << flush;

				SoSFVec3f position;
				position.setValue(p.x(),p.y(),p.z());
				//		((SoPerspectiveCamera*)camera)->position = position;
				camera->position = position;
				// camera->orientation.setValue(SbRotation::identity());
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

		offscreen_viz = true;
		if(offscreen_viz){
			cout << "Offscreen rendering and saving to a file." << endl << flush;
			cout << "width,height: " << width << ", " << height << endl << flush;
	//		OffScreenViz *offscreen = new OffScreenViz(800,600);
			OffScreenViz *offscreen = new OffScreenViz(width,height);
			// offscreen->initImageDirectory("tmp/offscreen/"+formattedDateTime()+"/"+baseName+"/",4);
			offscreen->initImageDirectory("/home/mrilee/workspace/VizHTM/tmp/offscreen/"+formattedDateTime()+"/"+baseName+"/",4);
			offscreen->root = new SoSeparator;
			offscreen->root->ref();
			SoCamera *camera = viewer->getCamera();
			SpatialVector *cam_vizCenter = new SpatialVector;
			SoSFVec3f position = camera->position;
			float x,y,z;
			position.getValue().getValue(x,y,z);
			cam_vizCenter->set(x, y, z);

			if( true ) {
				if(!lookFrom_flag) {
					cout << "Calling loadScene" << endl << flush;
					// loadScene(offscreen->root,selectionRoot,offscreen->vpRegion,cam_vizCenter,offscreenPerscamPositionScale,NULL);
					loadScene(offscreen->root,selectionRoot,offscreen->vpRegion,cam_vizCenter,offscreenPerscamPositionScale,camera);
				} else {
					cout << "Calling setLookFrom" << endl << flush;
					setLookFrom(offscreen->root,selectionRoot,offscreen->vpRegion);
				}
			} else {
				cout << "Loading test scene" << endl << flush;
				loadTestScene(offscreen->root,offscreen->vpRegion);
				cout << "Adding test scene to selectionRoot" << endl << flush;
				selectionRoot->addChild(offscreen->root);
			}
			// offscreen->writeFileRGB(1);
			offscreen->writeFile(1);
			// offscreen->saveImage(1);

		}

		delete viewer;
		//	root->unref();
	}

	cout << mainName << " done." << endl;

	return 0;
}


