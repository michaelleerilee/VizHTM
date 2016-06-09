/*
 * testViz.C
 *
 *  Created on: Mar 13, 2016
 *      Author: mrilee
 */

//#include <misc.h>
//#include <testOffScreenViz.h>

#include "misc.h"
#include "testOffScreenViz.h"

#include <stdlib.h>

#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>


#include <Inventor/SoDB.h>

#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoCube.h>
//#include <Inventor/nodes/>

#include <Inventor/SoOffscreenRenderer.h>
//#include <Inventor/>

namespace testOffScreenViz {

testOffScreenViz::testOffScreenViz() {
	loadCoinScene();

	std::string base     = "tmp/testOffScreenViz/"+formattedDateTime()+"/";
	std::string makePath = "mkdir -p " + base;
	int ok = system(makePath.c_str()); // TODO Check error conditions.

	for(int i=0; i<64; i++) {
		perscam->orientation.setValue(perscam->orientation.getValue()*SbRotation(SbVec3f(0, 1, 0), 0.01f));
		QImage image = getCoinCubeImgBuffer();
		image.save((formattedOutFileName(base, formattedZeroPaddedInteger(i), ".jpg")).c_str());
	}
}

testOffScreenViz::~testOffScreenViz() {
	// TODO Auto-generated destructor stub
}

void testOffScreenViz::loadCoinScene(){
    // Init Coin
    SoDB::init();
    // The root node
    root = new SoSeparator;
    root->ref();

    // Add the light _before_ you add the camera
    SoDirectionalLight * light = new SoDirectionalLight;
    root->addChild(light);

    vpRegion.setViewportPixels(0, 0, coinSceneWidth, coinSceneHeight);

//    SoPerspectiveCamera *
	perscam = new SoPerspectiveCamera();
    root->addChild(perscam);

//    SbRotation cameraRotation = SbRotation::identity();
    cameraRotation = new SbRotation(SbRotation::identity());
    (*cameraRotation) *= SbRotation(SbVec3f(0, 1, 0), 0.4f);
    perscam->orientation = *cameraRotation;

    SoCube * cube = new SoCube;
    root->addChild(cube);
    // make sure that the cube is visible
    perscam->viewAll(root, vpRegion);
    perscam->nearDistance = 0.1;
}

QImage testOffScreenViz::getCoinCubeImgBuffer(){
    SoOffscreenRenderer offscreenRenderer(vpRegion);
    offscreenRenderer.setComponents(
      SoOffscreenRenderer::Components::RGB_TRANSPARENCY
    );
    offscreenRenderer.render(root);

    QImage img(offscreenRenderer.getBuffer(), coinSceneWidth,
        coinSceneHeight, QImage::Format_ARGB32);

    // Important!
    return img.rgbSwapped();
}

} /* namespace testViz */
