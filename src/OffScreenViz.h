/*
 * OffScreenViz.h
 *
 *  Created on: Mar 14, 2016
 *      Author: mrilee
 */

#ifndef OFFSCREENVIZ_H_
#define OFFSCREENVIZ_H_

#include <QtOpenGL/QGL>
#include <QtGui/QImage>

#include <Inventor/SbRotation.h>
#include <Inventor/SbViewportRegion.h>

#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoSeparator.h>

#include "misc.h"
#include <string>
#include <sstream>

#include <iostream>
#include <fstream>

class OffScreenViz {
public:
	OffScreenViz();
	OffScreenViz(uint w, uint h) {
		sceneWidth = w; sceneHeight = h;
		init();
	}
	virtual ~OffScreenViz();

	void init();
	int initImageDirectory(const std::string basePath="tmp/", const uint fieldWidth=3);

	QImage getImage();
	void writeFile(const uint item, bool enableTransparency = false);
	void writeFileRGB(const uint item);
	void saveImage(const uint item);

	std::string base;
	uint        fieldWidth;

	std::stringstream *notes = NULL;
	void addNoteStream(std::stringstream *ss) { this->notes = ss; } ;
	void writeNotes();

	SoSeparator         *root;
	SbViewportRegion    *vpRegion;

	uint sceneWidth  = 400; // 480;
	uint sceneHeight = 400; // 270;

};

#endif /* OFFSCREENVIZ_H_ */
