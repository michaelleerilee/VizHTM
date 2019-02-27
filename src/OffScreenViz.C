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
#include <Inventor/lists/SbPList.h>
#include <Inventor/SbString.h>

#include <QImageWriter>

#include <iostream>
#include <string>
#include <unistd.h>
#include <stdio.h>

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

void OffScreenViz::writeFileRGB(const uint item) {
	bool ok = false;
	std::cout << "Trying SoOffscreenRenderer writeToFile RGB output" << std::endl << std::flush;
	SoOffscreenRenderer offscreenRenderer(*vpRegion);
	offscreenRenderer.setComponents(
			SoOffscreenRenderer::Components::RGB);
	//SoOffscreenRenderer::Components::RGB_TRANSPARENCY);
	offscreenRenderer.render(root);
	std::string fileExtension = "jpg";
	std::cout << 100 << std::endl << std::flush;
	SbViewportRegion vpr = offscreenRenderer.getViewportRegion();
	SbVec2s size = vpr.getViewportSizePixels();
	std::cout << "size: " << size.toString().getString() << std::endl << std::flush;
	SoOffscreenRenderer::Components components =  offscreenRenderer.getComponents();
	std::cout << "comp: " << components << std::endl << std::flush;
	unsigned char *buffer = offscreenRenderer.getBuffer();
	printf("buffer: %p\n",buffer);

	if( true ) {
	// if( offscreenRenderer.isWriteSupported(fileExtension.c_str()) ) {
		std::cout << 110 << std::endl << std::flush;
		ok = offscreenRenderer.writeToRGB((formattedOutFileName(base,formattedZeroPaddedInteger(item,fieldWidth),"-wf.rgb")).c_str());
		std::cout << 120 << "result: ok? = " << ok << std::endl << std::flush;
	}
}

void OffScreenViz::writeFile(const uint item) {
	// bool ok = false;
	std::cout << "Trying SoOffscreenRenderer writeToFile output" << std::endl << std::flush;
	SoOffscreenRenderer offscreenRenderer(*vpRegion);
	offscreenRenderer.setComponents(
			// SoOffscreenRenderer::Components::RGB);
			SoOffscreenRenderer::Components::RGB_TRANSPARENCY);

	// The following is okay, but won't reorder graphic elements.
	//	offscreenRenderer.getGLRenderAction()->setTransparencyType(SoGLRenderAction::DELAYED_BLEND);
	// Yay!!!  The following seems to work.
	offscreenRenderer.getGLRenderAction()->setTransparencyType(SoGLRenderAction::SORTED_OBJECT_SORTED_TRIANGLE_BLEND);
	offscreenRenderer.setBackgroundColor(SbColor(0.,0.,0.));

	// Um, first render fails in a tile.
	offscreenRenderer.render(root);
	offscreenRenderer.render(root);
	std::string fileExtension = "png";
	std::cout << 100 << std::endl << std::flush;
	if( true ) {
	// if( offscreenRenderer.isWriteSupported(fileExtension.c_str()) ) {
		std::cout << 110 << std::endl << std::flush;
		offscreenRenderer.writeToFile(
				(formattedOutFileName(base,formattedZeroPaddedInteger(item,fieldWidth),"-wf.png")).c_str()
				,"png"
		);
		std::cout << 120 << std::endl << std::flush;
	} else {
		std::cout << 200 << std::endl << std::flush;
		int nFiletypes = offscreenRenderer.getNumWriteFiletypes();
		std::cout << 210 << std::endl << std::flush;
		for( int i=0; i<nFiletypes; ++i ) {
			SbPList extlist; SbString fullname, description;
			offscreenRenderer.getWriteFiletypeInfo(i,extlist,fullname,description);
			std::cout << fullname.getString() << std::endl << std::flush;
					// " " << extlist <<
		}
		std::cout << 220 << std::endl << std::flush;
	}
}

void OffScreenViz::saveImage(const uint item) {
	bool ok = false;

	QImage image = getImage();
	// bool ok = image.save((formattedOutFileName(base,formattedZeroPaddedInteger(item,fieldWidth),".png")).c_str(),"PNG");
	ok = image.save((formattedOutFileName(base,formattedZeroPaddedInteger(item,fieldWidth),".jpg")).c_str(),"JPG");

	if (!ok) {
		// std::cout << "OffScreenViz::saveImage failed to save '" << (formattedOutFileName(base,formattedZeroPaddedInteger(item,fieldWidth),".png")).c_str() << "'" << std::endl << std::flush;
		std::cout << "OffScreenViz::saveImage failed to save '" << (formattedOutFileName(base,formattedZeroPaddedInteger(item,fieldWidth),".jpg")).c_str() << "'" << std::endl << std::flush;

		char cwd[1024];
		getcwd(cwd, sizeof(cwd));
		std::cout << "Current working dir: " << cwd << std::endl << std::flush;
	}

	// QImage::save or QPixmap::save :

	std::cout << "Trying QImageWriter..." << std::endl << std::flush;
	QImageWriter writer((formattedOutFileName(base,formattedZeroPaddedInteger(item,fieldWidth),"-iw.jpg")).c_str(),"JPG");
	if(!writer.write(image))
	{
		std::cout << writer.errorString().toStdString() << std::endl << std::flush;
	}
}

/*
QImage OffScreenViz::getImage() {
	SoOffscreenRenderer offscreenRenderer(*vpRegion);
	offscreenRenderer.setComponents(
			// SoOffscreenRenderer::Components::RGB);
			SoOffscreenRenderer::Components::RGB_TRANSPARENCY);
	offscreenRenderer.render(root);
	offscreenRenderer.render(root);

	printf("osv::getImage::buffer: %p\n",offscreenRenderer.getBuffer());

	QImage img(
			offscreenRenderer.getBuffer(),
			sceneWidth, sceneHeight,
			// QImage::Format_RGB32);
			QImage::Format_ARGB32);
	return img.rgbSwapped();
}
*/

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
			//QImage::Format_ARGB32);
	QImage img = imgFlipped.mirrored(false,true);
	return img.rgbSwapped();
	// return imgFlipped.rgbSwapped();
}

