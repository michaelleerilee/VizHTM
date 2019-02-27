/*
 * TestWindow.cpp
 *
 *  Created on: Feb 21, 2019
 *      Author: mrilee
 *
 *  Copyright (C) 2019 Rilee Systems Technologies LLC
 */

#include "TestWindow.h"

#include <QGraphicsScene>
#include <QGraphicsView>
#include <QImage>
#include <QPixmap>



TestWindow::TestWindow()
{
	ui = this;

	this->setGeometry(10,10,200,200);

	QGraphicsScene *scene = new QGraphicsScene(this);
	QGraphicsView  *view  = new QGraphicsView();
	view->setScene(scene);
	// scene->addPixmap(...);

}

TestWindow::~TestWindow() {
	delete ui;
}

