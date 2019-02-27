/*
 * TestWindow.h
 *
 *  Created on: Feb 21, 2019
 *      Author: mrilee
 *
 *  Copyright (C) 2019 Rilee Systems Technologies LLC
 */

#ifndef TESTWINDOW_H_
#define TESTWINDOW_H_

#include <QMainWindow>
#include <QImage>

class TestWindow : public QMainWindow
{
	Q_OBJECT

public:
	TestWindow();
	virtual ~TestWindow();

private:
	TestWindow *ui;
	QImage image;
};

#endif /* TESTWINDOW_H_ */
