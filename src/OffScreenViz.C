/*
 * OffScreenViz.C
 *
 *  Created on: Mar 14, 2016
 *      Author: mrilee
 */

//#include <OffScreenViz.h>
#include "OffScreenViz.h"

#include <Inventor/SoDB.h>
#include <Inventor/SoOffscreenRenderer.h>
#include <Inventor/nodes/SoTransparencyType.h>

OffScreenViz::OffScreenViz() {
	init();
}

void OffScreenViz::init() {
	// TODO Auto-generated constructor stub
//	if(!SoDB::isInitialized())SoDB::init(); // TODO Or should we fail?
	vpRegion = new SbViewportRegion(sceneWidth,sceneHeight);
//	vpRegion->setViewportPixels(0,0,sceneWidth,sceneHeight);
}

OffScreenViz::~OffScreenViz() {
	// TODO Auto-generated destructor stub
}

// TODO Add scene?

int OffScreenViz::initImageDirectory(const std::string basePath, const uint fieldWidth) {
	this->base       = basePath;
	this->fieldWidth = fieldWidth;
	std::string makePath = "mkdir -p " + base;
	return system(makePath.c_str());
}

void OffScreenViz::saveImage(const uint item) {
	QImage image = getImage();
	image.save((formattedOutFileName(base,formattedZeroPaddedInteger(item,fieldWidth),".png")).c_str()
			,"PNG");
//	image.save((formattedOutFileName(base,formattedZeroPaddedInteger(item,fieldWidth),".jpg")).c_str()
//			,"JPG");
}

QImage OffScreenViz::getImage() {
	SoOffscreenRenderer offscreenRenderer(*vpRegion);
	offscreenRenderer.setComponents(
//			SoOffscreenRenderer::Components::RGB);
			SoOffscreenRenderer::Components::RGB_TRANSPARENCY);
	// The following is okay, but won't reorder graphic elements.
	//	offscreenRenderer.getGLRenderAction()->setTransparencyType(SoGLRenderAction::DELAYED_BLEND);
	// Yay!!!  The following seems to work.
	offscreenRenderer.getGLRenderAction()->setTransparencyType(SoGLRenderAction::SORTED_OBJECT_SORTED_TRIANGLE_BLEND);
	offscreenRenderer.setBackgroundColor(SbColor(0.,0.,0.));
	offscreenRenderer.render(root);
	QImage imgFlipped(
			offscreenRenderer.getBuffer(),
			sceneWidth, sceneHeight,
			QImage::Format_RGB32);
//			QImage::Format_ARGB32);
	QImage img = imgFlipped.mirrored(false,true);
	return img.rgbSwapped();
}
