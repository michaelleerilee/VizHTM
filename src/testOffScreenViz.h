/*
 * testViz.h
 *
 *  Created on: Mar 13, 2016
 *      Author: mrilee
 */

#ifndef TESTOFFSCREENVIZ_H_
#define TESTOFFSCREENVIZ_H_

#include <QtOpenGL/QGL>
#include <QtGui/QImage>
//#include <QtGui/QRegion>

#include <Inventor/SbRotation.h>
#include <Inventor/SbViewportRegion.h>

#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoSeparator.h>

namespace testOffScreenViz {

/*
 *
 * http://stackoverflow.com/questions/20126354/convert-coin3d-sooffscreenrenderer-to-qimage-and-render-with-opengl
 *
 */

class testOffScreenViz {
public:
	testOffScreenViz();
	virtual ~testOffScreenViz();

	void loadCoinScene();
	QImage getCoinCubeImgBuffer();

	SoSeparator *root;
	SbViewportRegion vpRegion;
    uint coinSceneWidth = 400;
	uint coinSceneHeight = 200;

	SoPerspectiveCamera *perscam;
	SbRotation *cameraRotation;

};

} /* namespace testViz */

#endif /* TESTOFFSCREENVIZ_H_ */
