/*
 * VizHTM_main.C
 *
 *  Created on: Dec 16, 2015
 *      Author: mrilee
 */

#include "VizHTM.h"

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

// For offscreen rendering.
#include "misc.h"
#include "OffScreenViz.h"

#include "VizHTM_main.h"

using namespace std;

float lineWidth = -1;

/**
 *
 * @param htmIdLevel
 * @param range1
 * @param range2
 * @return
 */
HtmRange *HTMRangeAtLevelFromIntersection(int htmIdLevel, HtmRange *range1, HtmRange *range2){
	//	cout << "Comparing..." << endl << flush;
	if((!range1)||(!range2)) return 0;
	if((range1->nranges()<=0)||(range2->nranges()<=0)) return 0;
	HtmRange *resultRange = new HtmRange();
	resultRange->purge();
	Key lo1,hi1,lo2,hi2;
	range1->reset();
	uint64 indexp1 = range1->getNext(lo1,hi1);
	if (!indexp1) return 0;
//	cout << "indexp1: " << indexp1 << endl << flush;
//	cout << "l,lo,hi1: " << htmIdLevel << " " << lo1 << " " << hi1 << endl << flush;
//	cout << "a" << flush;
	do {
//		cout << "b" << endl << flush;
		KeyPair testRange1 = HTMRangeAtLevelFromHTMRange(htmIdLevel,lo1,hi1);
		range2->reset();
		uint64 indexp2 = range2->getNext(lo2,hi2);
		if(!indexp2) return 0;
//		cout << "indexp2: " << indexp2 << endl << flush;
//		cout << "l,lo,hi2: " << htmIdLevel << " " << lo2 << " " << hi2 << endl << flush;
		bool intersects = false;
		do {
//			cout << "c" << endl << flush;
			KeyPair testRange2 = HTMRangeAtLevelFromHTMRange(htmIdLevel,lo2,hi2);
			intersects = testRange2.lo <= testRange1.hi
					&& testRange2.hi >= testRange1.lo;
//						cout << "lh1,lh2: "
//								<< lo1 << " " << hi1 << ", "
//								<< lo2 << " " << hi2 << ", "
//								<< intersects << flush;
			if(intersects){
				Key lo_ = max(testRange1.lo,testRange2.lo);
				Key hi_ = min(testRange1.hi,testRange2.hi);
				resultRange->addRange(lo_,hi_);
//								cout << ", added "
//										<< lo_ << " "
//										<< hi_ << flush;
			}
//						cout << "." << endl << flush;
		} while (range2->getNext(lo2,hi2));
	} while (range1->getNext(lo1,hi1));
//	cout << "d" << flush;
	if(resultRange->nranges()>0)resultRange->defrag();
//	cout << "e" << flush;
	return resultRange;
}


void testAnEdge(VizHTM *viz) {
	int colorBase = viz->nEdgeColors;
	viz->addEdgeColor(0,1,1);
	viz->addEdgeColor(0,1,1);
	viz->addEdgeVertexColorIndices(colorBase,colorBase+1);
	int coordBase = viz->nCoordinates;
	viz->addEdgeIndices(coordBase,coordBase+1);
}

void testConstraint(VizHTM *viz, SpatialConstraint sc) {
	viz->addConstraint(sc,1.0,0.5,1.0);
}

void testConstraintAD(VizHTM *viz, SpatialVector a, float64 d) {
	SpatialConstraint sc = SpatialConstraint(a,d);
//	cout << " a in sc(a,d)? " << sc.contains(a) << endl << flush;
//	cout << " acos(a*a),s " << acos(a*a) << " " << sc.coneAngle() << endl << flush;
//	cout << " i in sc(a,d)? " << sc.contains(SpatialVector(1.0,0.0,0.0)) << endl << flush;
//	cout << " acos(i*a),s " << acos(SpatialVector(1.0,0.0,0.0)*a) << " " << sc.coneAngle() << endl << flush;
	testConstraint(viz,sc);
}

void testConvexHtmRangeIntersection(VizHTM * viz, RangeConvex convex, int htmIdLevel) {

	int saveLevel = 5;
	SpatialIndex *index = new SpatialIndex(htmIdLevel,saveLevel);
	SpatialDomain domain = SpatialDomain(index);
	domain.add(convex);
//	domain.setOlevel(htmIdLevel);

	cout << "intersection" << endl << flush;

	HtmRange *range = new HtmRange();
	HtmRange boundary, interior;

	range->purge();
	bool varlen_individualHTMIds = false; // true for individuals, false for ranges
	bool overlap = domain.intersect(index,range,varlen_individualHTMIds,&interior,&boundary);

	Key lo = 0; Key hi = 0;
	SpatialVector v1,v2,v3;
	if(true)
		cout
		<< " overlap: " << overlap
		<< " nRanges: " << range->nranges() << flush;
	range->defrag();
//	if(false)
//		cout
//		<< " nRangeDefrag: " << range->nranges()
//		<< " nConvexes: " << domain.numConvexes()
//		<< endl << flush;

	//		overlap = false;

	if(overlap) {
		range->reset();
		uint64 indexp = range->getNext(lo,hi);
		//		indexp = range->getNext(lo,hi);
		// cout << " indexp,lo,hi: " << indexp << " " << lo << " " << hi << endl << flush;
		int k = 0; int kMax = 10000;
		int triangles = 0; int tMax = 10000;
		do {
			k++;
			//			for(uint64 nodeIndex=lo; nodeIndex <= min(lo+10,hi); nodeIndex++){
			for(uint64 numericId=lo; numericId <= hi && triangles < tMax; numericId++){
				triangles ++;
				uint64 nodeIndex = index->nodeIndexFromId(numericId);

//				if(false){
//					cout
//					<< "working " << flush
//					<< " triangles: " << triangles << flush
//					<< " numericId: " << numericId << flush
//					<< " nodeIndex: " << nodeIndex << flush;
//					cout << " k=" << k << flush;
//				}

				index->nodeVertex(nodeIndex,v1,v2,v3);

//				if(false){
//					cout << "..." << flush;
//					cout
//					<< " v1: " << v1.x() << " " << v1.y() << " " << v1.z()
//					<< " v2: " << v2.x() << " " << v2.y() << " " << v2.z()
//					<< " v3: " << v1.x() << " " << v3.y() << " " << v3.z() << flush;
//				}

				if(true){
					float size = pow(0.5,htmIdLevel+3);
					float r = 0.4, g = 0.4, b = 0.6;
					SpatialVector x = 3.*v1+v2+v3; x.normalize(); x *= 1.0+1.0e-6;
					SpatialVector x_ = v1+v2+v3; x_.normalize();
//					if(false){
//						cout
//						<< " nI: " << nodeIndex
//						<< " x: " << x.x() << " " << x.y() << " " << x.z();
//					}
					char *str = new char[256];
					sprintf(str,"id: %llx\nix: %llu\n",numericId,nodeIndex);
					viz->addAnnotation((new SpatialVector(x)),str,size,r,g,b);
					viz->addEdge(x,x_,0.1,0.1,0.1);
				}

				int indexBase = viz->nCoordinates;
				int colorBase = viz->nFaceColors;
				if(true) {
					float r=0.; float g=0.; float b=0.;
					switch(numericId % 4) {
					case 0:
						r=1.;
						break;
					case 1:
						g=1.;
						break;
					case 2:
						b=1.;
						break;
					default:
						r=1.; g=1.; b=1.;
						break;
					}
					//	for(int i=0; i<3; i++) addEdgeColor(r,g,b);
					for(int i=0; i<3; i++) {
						viz->addFaceColor(r,g,b);
					}
					viz->addFaceVertexColorIndices3(colorBase,colorBase+1,colorBase+2);

					viz->addCoordinate64(v1.x(),v1.y(),v1.z());
					viz->addCoordinate64(v2.x(),v2.y(),v2.z());
					viz->addCoordinate64(v3.x(),v3.y(),v3.z());
					//	addEdgeIndicesTriangle(indexBase,indexBase+1,indexBase+2);
					viz->addFaceIndices3(indexBase,indexBase+1,indexBase+2);
					//	cout << " ( " << indexBase << " " << indexBase+1 << " " << indexBase+2 << ") ";
				}
//				if(false) cout << endl << flush;
			}
		} while (range->getNext(lo,hi) && k < kMax && triangles < tMax);
	}
	//	cout << "overlap done " << endl << flush;



}

void testTwoConstraints(VizHTM *viz, int htmIdLevel) {
	cout << 100 << endl << flush;
	SpatialVector o = SpatialVector(0.0,0.0,0.0);

	//	SpatialVector a = SpatialVector(1.0,0.0,0.0);
	SpatialVector a = randomVector();
	//	SpatialVector a = SpatialVector(1.0,1.0,0.0); a.normalize();

	SpatialVector as = 1.05*a;

	cout << 300 << endl << flush;

	//	float64 d = 0.9999;
	float64 d = 0.99;
	SpatialConstraint sc = SpatialConstraint(a,d);

	if(false) testConstraintAD(viz,a,d);

	RangeConvex convex = RangeConvex();
	convex.add(sc);

	//	SpatialVector a1 = a+randomVector()*0.15; a1.normalize();
	//	SpatialConstraint sc1 = SpatialConstraint(a1,0.99);

	cout << 400 << endl << flush;

	SpatialVector a1 = -1*a; a1.normalize();
	SpatialConstraint sc1 = SpatialConstraint(a1,0.995);
	sc1.invert();

	cout << 600 << endl << flush;

	// SpatialConstraint(randomVector(),0.99);
	viz->addConstraint(sc1,0.5,1.0,1.0);
	cout << 700 << endl << flush;
	convex.add(sc1);
	cout << 800 << endl << flush;

	testConvexHtmRangeIntersection(viz,convex,htmIdLevel);
	cout << 999 << endl << flush;
}

void testAddRectangle0(VizHTM *viz) {

	SpatialVector *v0 = VectorFromLatLonDegrees(5.0,0.0);
	SpatialVector *v1 = VectorFromLatLonDegrees(85.0,0.0);
	SpatialVector *v2 = VectorFromLatLonDegrees(85.0,85.0);
	SpatialVector *v3 = VectorFromLatLonDegrees(5.0,85.0);
	float r = 0.9;
	float g = 0.9;
	float b = 0.1;
	float a = -1.0;

	// viz->addArc(v1,v2,r,g,b,a,scale_,steps);
	float64 scale_ = 1.02;
	int steps = 100;
	viz->addArc(*v0,*v1,r,g,b,a,scale_,steps);
	viz->addArc(*v1,*v2,r,g,b,a,scale_,steps);
	viz->addArc(*v2,*v3,r,g,b,a,scale_,steps);
	viz->addArc(*v3,*v0,r,g,b,a,scale_,steps);

	r = 0.1;
	g = 0.9;
	b = 0.9;
	a = -1.0;

	viz->addRectangle(*v0,*v1,*v2,*v3,r,g,b);
	viz->addArcAtLatitudeDegrees(5.0,0.,85.0,1.0,0.,0.);
	viz->addArcAtLatitudeDegrees(85.0,0.,85.0,1.0,0.,0.);
}

void testAddRectangle(VizHTM *viz, int htmIdLevel) {
	int saveLevel = 5;

	SpatialIndex *index = new SpatialIndex(htmIdLevel,saveLevel);
	SpatialDomain domain = SpatialDomain(index);
//	domain.setOlevel(htmIdLevel); // Note this sets the olevel on the convexes.

	SpatialVector *v0 = VectorFromLatLonDegrees(10.0,0.0);
	SpatialVector *v1 = VectorFromLatLonDegrees(70.0,0.0);
	SpatialVector *v2 = VectorFromLatLonDegrees(70.0,80.0);
	SpatialVector *v3 = VectorFromLatLonDegrees(10.0,80.0);
	float r = 0.1;
	float g = 0.1;
	float b = 0.9;

	RangeConvex *rc = new RangeConvex(v0,v1,v2,v3);
	//	cout << "nConstraints: " << rc->numConstraints() << endl << flush;
	//	SpatialConstraint *sc = new SpatialConstraint(SpatialVector(0.,0.,1.),0.5);
	//	viz->addConstraint(*sc,1.0,1.0,1.0);
	//	rc->add(*sc);
//	rc->setOlevel(htmIdLevel); // Note this is supposed to be done when added to the domain. ///???
	domain.add(*rc);

	viz->addRectangle(*v0,*v1,*v2,*v3,r,g,b);
	viz->addArcAtLatitudeDegrees(10.0,0.,80.0,1.0,0.,0.);
	viz->addArcAtLatitudeDegrees(70.0,0.,80.0,1.0,0.,0.);

	//	return;

	HtmRange *range = new HtmRange();
	range->purge();
	bool varlen_individualHTMIds = false; // true for individuals, false for ranges
	bool overlap = domain.intersect(index,range,varlen_individualHTMIds);
	range->defrag();
	range->reset();
	Key lo = 0, hi = 0;
	uint64 indexp = range->getNext(lo,hi);
	SpatialVector x1,x2,x3;
	do {
		KeyPair adjustedRange = HTMRangeAtLevelFromHTMRange(htmIdLevel,lo,hi);

		lo = adjustedRange.lo;
		hi = adjustedRange.hi;

		for(uint64 numericId=lo; numericId<=hi;numericId++) {
			uint64 nodeIndex = index->nodeIndexFromId(numericId);
			if(nodeIndex!=0){
				index->nodeVertex(nodeIndex,x1,x2,x3);
				if(true) {
					float r=0.; float g=0.; float b=0.;
					switch(numericId % 4) {
					case 0:
						r=1.;
						break;
					case 1:
						g=1.;
						break;
					case 2:
						b=1.;
						break;
					default:
						r=1.; g=1.; b=1.;
						break;
					}
					//	for(int i=0; i<3; i++) addEdgeColor(r,g,b);
					int colorBase = viz->nFaceColors;
					for(int i=0; i<3; i++) {
						viz->addFaceColor(r,g,b);
					}
					viz->addFaceVertexColorIndices3(colorBase,colorBase+1,colorBase+2);

					int indexBase = viz->nCoordinates;
					viz->addCoordinate64(x1.x(),x1.y(),x1.z());
					viz->addCoordinate64(x2.x(),x2.y(),x2.z());
					viz->addCoordinate64(x3.x(),x3.y(),x3.z());
					viz->addFaceIndices3(indexBase,indexBase+1,indexBase+2);

					//					printf("id: %llx ix: %llu",numericId,nodeIndex);
					//					cout << endl << flush;

					if(true){
						float size = pow(0.5,htmIdLevel+3);
						float r = 0.4, g = 0.4, b = 0.6;
						SpatialVector x = 3.*x1+x2+x3; x.normalize(); x *= 1.0+1.0e-6;
						SpatialVector x_ = x1+x2+x3; x_.normalize();
//						if(false){
//							cout
//							<< " nI: " << nodeIndex
//							<< " x: " << x.x() << " " << x.y() << " " << x.z();
//						}
						char *str = new char[256];
						sprintf(str,"id: %llx\nix: %llu\n",numericId,nodeIndex);
						viz->addAnnotation((new SpatialVector(x)),str,size,r,g,b);
						viz->addEdge(x,x_,0.5,0.5,0.9);
					}
				}
			}
		}
	} while (range->getNext(lo,hi));
}

void testTenDegreeGridRGB(VizHTM *viz, float r, float g, float b) {
	float64 PI = atan2(0,-1);
	float64 k  = 2*PI/360.;

	float dLat = 10;
	float dLon = 10;
	for(float lat=-90; lat<=90; lat+= dLat) {
		for(float lon=0; lon<=360; lon+= dLon) {
			viz->addLatLonBoxEdgesDegrees(
					lat + 0.25, lon + 0.25,
					lat + dLat - 0.25, lon + dLon - 0.25,
					r,g,b);
		}
	}
}
/*
void plotHTMRange(VizHTM *viz, SpatialIndex index, HtmRange range) {
	range.purge();
	bool varlen_individualHTMIds = false; // true for individuals, false for ranges
	bool overlap = domain.intersect(index,range,varlen_individualHTMIds);
	range.defrag();
	range.reset();
	Key lo = 0, hi = 0;
	uint64 indexp = range.getNext(lo,hi);
	SpatialVector x1,x2,x3;
	do {
		KeyPair adjustedRange = HTMRangeAtLevelFromHTMRange(htmIdLevel,lo,hi);

		lo = adjustedRange.lo;
		hi = adjustedRange.hi;

		for(uint64 numericId=lo; numericId<=hi;numericId++) {
			uint64 nodeIndex = index.nodeIndexFromId(numericId);
			if(nodeIndex!=0){
				index.nodeVertex(nodeIndex,x1,x2,x3);
				if(true) {
					float r=0.; float g=0.; float b=0.;
					switch(numericId % 4) {
					case 0:
						r=1.;
						break;
					case 1:
						g=1.;
						break;
					case 2:
						b=1.;
						break;
					default:
						r=1.; g=1.; b=1.;
						break;
					}
					//	for(int i=0; i<3; i++) addEdgeColor(r,g,b);
					int colorBase = viz->nFaceColors;
					for(int i=0; i<3; i++) {
						viz->addFaceColor(r,g,b);
					}
					viz->addFaceVertexColorIndices3(colorBase,colorBase+1,colorBase+2);

					int indexBase = viz->nCoordinates;
					viz->addCoordinate64(x1.x(),x1.y(),x1.z());
					viz->addCoordinate64(x2.x(),x2.y(),x2.z());
					viz->addCoordinate64(x3.x(),x3.y(),x3.z());
					viz->addFaceIndices3(indexBase,indexBase+1,indexBase+2);

//					printf("id: %llx ix: %llu",numericId,nodeIndex);
//					cout << endl << flush;

					if(true){
						float size = pow(0.5,htmIdLevel+3);
						float r = 0.4, g = 0.4, b = 0.6;
						SpatialVector x = 3.*x1+x2+x3; x.normalize(); x *= 1.0+1.0e-6;
						SpatialVector x_ = x1+x2+x3; x_.normalize();
						if(false){
							cout
							<< " nI: " << nodeIndex
							<< " x: " << x.x() << " " << x.y() << " " << x.z();
						}
						char *str = new char[256];
						sprintf(str,"id: %llx\nix: %llu\n",numericId,nodeIndex);
						viz->addAnnotation((new SpatialVector(x)),str,size,r,g,b);
						viz->addEdge(x,x_,0.5,0.5,0.9);
					}
				}
			}
		}
	} while (range.getNext(lo,hi));
}
 */


void intersectTwoRectangles(
		VizHTM *const viz,
		const SpatialIndex  *index,
		const SpatialVector *u0,
		const SpatialVector *u1,
		const SpatialVector *u2,
		const SpatialVector *u3,
		const SpatialVector *v0,
		const SpatialVector *v1,
		const SpatialVector *v2,
		const SpatialVector *v3,
		HtmRange *rangeU,
		HtmRange *rangeV,
		HtmRange *rangeIntersection
		) {
	int htmIdLevel = index->getLeafLevel();

	SpatialDomain domain1 = SpatialDomain(index);
	SpatialDomain domain2 = SpatialDomain(index);
	cout << "1" << flush;
	if(true){
		RangeConvex rc = RangeConvex(u0,u1,u2,u3);
		domain1.add(rc);
	}
	cout << "2" << flush;
	if(true){
		RangeConvex rc = RangeConvex(v0,v1,v2,v3);
		domain2.add(rc);
	}
//		return;
	cout << "3" << flush;
	bool varlen_individualHTMIds = false; // true for individuals, false for ranges

	cout << "." << flush;
//	HtmRange range1 = HtmRange();
	HtmRange range1;
	range1.purge();
	cout << "." << flush;
	bool overlap1 = domain1.intersect(index,&range1,varlen_individualHTMIds);
	cout << "." << flush;
	range1.defrag();
	cout << "." << flush;
	range1.reset();
	cout << "." << flush;
	//	if(range1.nranges()==0)return;
	rangeU->addRange(&range1);
	cout << "." << flush;
	range1.reset();
	cout << "." << flush;
	cout << "4" << flush;
	cout << " overlap1: " << overlap1 << ";" << flush;
	cout << endl << flush;

	Key lo = 0, hi = 0;
	uint64 indexp = 0;

	cout << " range1.ranges(): " << range1.nranges() << endl << flush;
	range1.reset();
	indexp = range1.getNext(lo,hi);
	cout << " range1indexp " << indexp << endl << flush;
	cout << " range1.lo,hi " << lo << " " << hi << endl << flush;
	cout << "        level " << levelOfId(lo) << endl << flush;

//	cout << "." << flush;
	HtmRange range2 = HtmRange();
	range2.purge();
//	cout << "." << flush;
	bool overlap2 = domain2.intersect(index,&range2,varlen_individualHTMIds);
//	cout << "." << flush;
	range2.defrag();
	range2.reset();
//	cout << "." << flush;
//	if(range2.nranges()==0)return;
	rangeV->addRange(&range2);
//	cout << "." << flush;
//	rangeV->defrag();
//	cout << "." << flush;
	range2.reset();


	cout << "5" << flush;

//	cout << " overlap2: " << overlap2 << ";" << flush;
//	cout << endl << flush;
//
//	cout << " range2.ranges(): " << range2.nranges() << endl << flush;
	range2.reset();
	indexp = range2.getNext(lo,hi);
//	cout << " range2indexp " << indexp << endl << flush;
//	cout << " range2.lo,hi " << lo << " " << hi << endl << flush;
//	cout << "        level " << levelOfId(lo) << endl << flush;

	range1.reset(); range2.reset();

	if(range1.nranges()*range2.nranges()==0) return;
	cout << "6" << flush;
//	HtmRange *resultRange = HTMRangeAtLevelFromIntersection(htmIdLevel,&range1,&range2);
	HtmRange *resultRange = range1.HTMRangeAtLevelFromIntersection(&range2,htmIdLevel);
	if(!resultRange) return;
	cout << "7" << flush;
	rangeIntersection->addRange(resultRange); // Note: resultRange is copied piece by piece here.
	cout << "8" << flush;
	rangeIntersection->defrag();
	cout << "9" << flush;
//	cout << "+" << range1.nranges() << "," << range2.nranges() << flush;

	HtmRange *range = resultRange;
//	viz->addHTMRange(index,range,0.,0.7,0.7);

//	HtmRange *red = new HtmRange;
//	HtmRange *blue = new HtmRange;
//	HtmRange *green = new HtmRange;
//	range1.reset();
//	indexp = range1.getNext(lo,hi);
//	if(indexp) {
//		do {
////			cout << "r2.contains: " << lo << " " << hi << " " << range2.contains(lo,hi) << endl << flush;
//			int ret = range2.contains(lo,hi);
//			if(ret == -1) {
//				red->addRange(lo,hi);
//			} else if (ret == 0 ) {
//				green->addRange(lo,hi);
//			} else {
//				blue->addRange(lo,hi);
//			}
//		} while(range1.getNext(lo,hi));
//		viz->addHTMRange(index,red,0.7,0.,0.);
//		viz->addHTMRange(index,green,0.,0.7,0.);
//		viz->addHTMRange(index,blue,0.,0.,0.7);
//		cout << "red" << endl;
//		red->print(HtmRange::BOTH,cout);
//		cout << "green" << endl;
//		green->print(HtmRange::BOTH,cout);
//		cout << "blue" << endl;
//		blue->print(HtmRange::BOTH, cout);
//	}
//  delete red, blue, green;
//
//	range1.reset();
////	HtmRange *range1_ = &range1;
////	viz->addHTMRange(index,&range1,0.7,0.,0.);

	if(false) {
		range->reset();
		//	cout << "6" << flush;
		lo=0; hi=0;
		indexp = range->getNext(lo,hi);
		SpatialVector x1,x2,x3;

		//	cout << "7" << flush;
		//		if(indexp) //?
		//	cout << "lo,hi: " << lo << " " << hi << endl << flush;
		//	cout << "indexp: " << indexp << endl << flush;
		if(indexp) //?
			do {
				for(uint64 numericId=lo; numericId<=hi;numericId++) {
					uint64 nodeIndex = index->nodeIndexFromId(numericId);
					if(nodeIndex!=0){
						viz->addEdgesFromIndexAndId(index,numericId,0.2,1.0,0.2);
						if(false){
							index->nodeVertex(nodeIndex,x1,x2,x3);
							if(true) {
								float r=0.; float g=0.; float b=0.;
								switch(numericId % 4) {
								case 0:
									r=1.;
									break;
								case 1:
									g=1.;
									break;
								case 2:
									b=1.;
									break;
								default:
									r=1.; g=1.; b=1.;
									break;
								}
								//					for(int i=0; i<3; i++) viz->addEdgeColor(r,g,b);
								int colorBase = viz->nFaceColors;
								for(int i=0; i<3; i++) {
									viz->addFaceColor(r,g,b);
								}
								viz->addFaceVertexColorIndices3(colorBase,colorBase+1,colorBase+2);

								int indexBase = viz->nCoordinates;
								viz->addCoordinate64(x1.x(),x1.y(),x1.z());
								viz->addCoordinate64(x2.x(),x2.y(),x2.z());
								viz->addCoordinate64(x3.x(),x3.y(),x3.z());
								viz->addFaceIndices3(indexBase,indexBase+1,indexBase+2);
							}
							//					printf("id: %llx ix: %llu",numericId,nodeIndex);
							//					cout << endl << flush;

							if(true){
								float size = pow(0.5,htmIdLevel+3);
								float r = 0.4, g = 0.4, b = 0.6;
								SpatialVector x = 3.*x1+x2+x3; x.normalize(); x *= 1.0+1.0e-6;
								SpatialVector x_ = x1+x2+x3; x_.normalize();
								if(false){
									cout
									<< " nI: " << nodeIndex
									<< " x: " << x.x() << " " << x.y() << " " << x.z();
								}
								char *str = new char[256];
								sprintf(str,"id: %llx\nix: %llu\n",numericId,nodeIndex);
								viz->addAnnotation((new SpatialVector(x)),str,size,r,g,b);
								viz->addEdge(x,x_,0.5,0.5,0.9);
							}
						}
					}
				}
			} while (range->getNext(lo,hi));
	}

	delete resultRange; // Should replace with a local var.

}


void intersectTwoRectangles_TRMMAndNMQ_1(
		VizHTM *const viz,
		const SpatialIndex  *index,
		const SpatialVector *u0,
		const SpatialVector *u1,
		const SpatialVector *u2,
		const SpatialVector *u3,
		const SpatialVector *v0,
		const SpatialVector *v1,
		const SpatialVector *v2,
		const SpatialVector *v3,
		HtmRange *rangeU,
		HtmRange *rangeV,
		HtmRange *rangeIntersection
		) {
	int htmIdLevel = index->getLeafLevel();

	SpatialDomain domain1 = SpatialDomain(index);
	SpatialDomain domain2 = SpatialDomain(index);
	if(true){
		RangeConvex rc = RangeConvex(u0,u1,u2,u3);
		domain1.add(rc);
	}
	if(true){
		RangeConvex rc = RangeConvex(v0,v1,v2,v3);
		domain2.add(rc);
	}
	bool varlen_individualHTMIds = false; // true for individuals, false for ranges

	HtmRange range1;
	range1.purge();
	bool overlap1 = domain1.intersect(index,&range1,varlen_individualHTMIds);
	range1.defrag();
	range1.reset();
	rangeU->addRange(&range1);
	range1.reset();

	Key lo = 0, hi = 0;
	uint64 indexp = 0;

	range1.reset();
	indexp = range1.getNext(lo,hi);

	HtmRange range2 = HtmRange();
	range2.purge();
	bool overlap2 = domain2.intersect(index,&range2,varlen_individualHTMIds);
	range2.defrag();
	range2.reset();
	rangeV->addRange(&range2);
	range2.reset();

	range2.reset();
	indexp = range2.getNext(lo,hi);

	range1.reset(); range2.reset();

	if(range1.nranges()*range2.nranges()==0) return;

	HtmRange *resultRange = range1.HTMRangeAtLevelFromIntersection(&range2,htmIdLevel);
	if(!resultRange) return;

	rangeIntersection->addRange(resultRange); // Note: resultRange is copied piece by piece here.

	rangeIntersection->defrag();

	HtmRange *range = resultRange;

	if(false) {
		range->reset();
		//	cout << "6" << flush;
		lo=0; hi=0;
		indexp = range->getNext(lo,hi);
		SpatialVector x1,x2,x3;

		//	cout << "7" << flush;
		//		if(indexp) //?
		//	cout << "lo,hi: " << lo << " " << hi << endl << flush;
		//	cout << "indexp: " << indexp << endl << flush;
		if(indexp) //?
			do {
				for(uint64 numericId=lo; numericId<=hi;numericId++) {
					uint64 nodeIndex = index->nodeIndexFromId(numericId);
					if(nodeIndex!=0){
						viz->addEdgesFromIndexAndId(index,numericId,0.2,1.0,0.2);
						if(false){
							index->nodeVertex(nodeIndex,x1,x2,x3);
							if(true) {
								float r=0.; float g=0.; float b=0.;
								switch(numericId % 4) {
								case 0:
									r=1.;
									break;
								case 1:
									g=1.;
									break;
								case 2:
									b=1.;
									break;
								default:
									r=1.; g=1.; b=1.;
									break;
								}
								//					for(int i=0; i<3; i++) viz->addEdgeColor(r,g,b);
								int colorBase = viz->nFaceColors;
								for(int i=0; i<3; i++) {
									viz->addFaceColor(r,g,b);
								}
								viz->addFaceVertexColorIndices3(colorBase,colorBase+1,colorBase+2);

								int indexBase = viz->nCoordinates;
								viz->addCoordinate64(x1.x(),x1.y(),x1.z());
								viz->addCoordinate64(x2.x(),x2.y(),x2.z());
								viz->addCoordinate64(x3.x(),x3.y(),x3.z());
								viz->addFaceIndices3(indexBase,indexBase+1,indexBase+2);
							}
							//					printf("id: %llx ix: %llu",numericId,nodeIndex);
							//					cout << endl << flush;

							if(true){
								float size = pow(0.5,htmIdLevel+3);
								float r = 0.4, g = 0.4, b = 0.6;
								SpatialVector x = 3.*x1+x2+x3; x.normalize(); x *= 1.0+1.0e-6;
								SpatialVector x_ = x1+x2+x3; x_.normalize();
								if(false){
									cout
									<< " nI: " << nodeIndex
									<< " x: " << x.x() << " " << x.y() << " " << x.z();
								}
								char *str = new char[256];
								sprintf(str,"id: %llx\nix: %llu\n",numericId,nodeIndex);
								viz->addAnnotation((new SpatialVector(x)),str,size,r,g,b);
								viz->addEdge(x,x_,0.5,0.5,0.9);
							}
						}
					}
				}
			} while (range->getNext(lo,hi));
	}

	delete resultRange; // Should replace with a local var.

}


typedef float (*colorMap) (uint64 id, uint64 loId, uint64 hiId, float lo, float hi );

float linearInterp(uint64 id, uint64 loId, uint64 hiId, float lo, float hi ) {
	if (hiId == loId ) return 0.5*(lo+hi);
	return lo + (hi-lo)*((float) (id-loId))/((float) (hiId-loId));
}
float hat(uint64 id, uint64 loId, uint64 hiId, float lo, float hi) {
	if (hiId == loId ) return 0.5*(lo+hi);
	float x = (float) (id-loId)/ (float) (hiId-loId); // map to [0,1]
	float d = hi-lo;
	return lo + 2.0 * d * (0.5-abs( x - 0.5 ));
}

void plotEdgesOrArcFromHTMNameInterval(
		VizHTM *viz, const char* loName, const char* hiName,
		float r0, float g0, float b0,
		float r1, float g1, float b1,
		colorMap cmr, colorMap cmg, colorMap cmb,
		bool edges = true
		) {

	SpatialIndex *index = new SpatialIndex(loName);
	int indexLevel = index->getMaxlevel();
	if (indexLevel != levelOfName(hiName)) {
		cout
		<< "levelOfName(loName) " << levelOfName(loName)
		<< " != levelOfName(hiName) " << levelOfName(hiName)
		<< endl << flush;
		return;
	}
	uint64 loId = index->idByName(loName);
	uint64 hiId = index->idByName(hiName);
//	cout << "loId,hiId: " << loId << " " << hiId << endl << flush;
	for(uint64 id = loId; id <= hiId; id++) {
		if(edges) {
			viz->addEdgesFromIndexAndId(
					index,id,
					cmr(id,loId,hiId,r0,r1),
					cmg(id,loId,hiId,g0,g1),
					cmb(id,loId,hiId,b0,b1),
					0.1);
		} else {
			viz->addArcFromIndexAndId(
					index,id,
					cmr(id,loId,hiId,r0,r1),
					cmg(id,loId,hiId,g0,g1),
					cmb(id,loId,hiId,b0,b1),
					0.1);
		}
	}

//	viz->addEdgesFromIndexAndName(index,loName,r0,g0,b0);
//	viz->addEdgesFromIndexAndName(index,hiName,r1,g1,b1);

	delete index;
}

void plotBlockingSphere(VizHTM* viz, float r, float g, float b, float radius) {
	viz->addSphere(SpatialVector(0.,0.,0.),r,g,b,radius);
}

void testPlotEdgesFromHTMNameInterval(
		VizHTM *viz
	) {

	viz->lineWidth = 1.5;
	viz->sphereComplexity = 1.0;
	plotBlockingSphere(viz,0.3,0.3,0.3,0.999999);
	plotEdgesOrArcFromHTMNameInterval(
			viz,"S000","N333",
			1.,0.,0.,
			1.,0.,0.,
			linearInterp,hat,linearInterp,
			false
			);
//	plotEdgesOrArcFromHTMNameInterval(
//			viz,"N030","N030",
//			1.,1.,1.,
//			1.,1.,1.,
//			linearInterp,hat,linearInterp,
//			false
//			);
	plotEdgesOrArcFromHTMNameInterval(
			viz,"N0300","N0333",
			1.,1.,0.,
			1.,1.,0.,
			linearInterp,hat,linearInterp,
			false
			);
//	plotEdgesOrArcFromHTMNameInterval(
//			viz,"N0303","N0303",
//			0.,1.,1.,
//			0.,1.,1.,
//			linearInterp,hat,linearInterp,
//			false
//			);
	plotEdgesOrArcFromHTMNameInterval(
			viz,"N03330","N03333",
			1.,0.,1.,
			1.,0.,1.,
			linearInterp,hat,linearInterp,
			false
			);
	plotEdgesOrArcFromHTMNameInterval(
			viz,"N0333330","N0333333",
			.5,5.,1.,
			.5,5.,1.,
			linearInterp,hat,linearInterp,
			true
			);
	plotEdgesOrArcFromHTMNameInterval(
			viz,"N033333330","N033333333",
			1.,0.,1.,
			1.,0.,1.,
			linearInterp,hat,linearInterp,
			true
			);
}


void testPlotEdgesFromHTMNameInterval0(
		VizHTM *viz
	) {
	plotBlockingSphere(viz,0.1,0.1,0.1,0.999);
	plotEdgesOrArcFromHTMNameInterval(
			viz,"N3333000000","N3333333333",
			1.,0.,0.,
			0.,1.,1.,
			linearInterp,hat,linearInterp
			);
}

SpatialVector *pointFromReferenceAndDeltas(
		const SpatialVector target,
		const SpatialVector aheadOfTarget,
		const float64 deltaAlongDegrees,
		const float64 deltaCrossDegrees) {
	float64 k = 2.*M_PI/360.0;
	SpatialVector L = target^aheadOfTarget;
	return new SpatialVector(
			sin(k*deltaCrossDegrees)*L
			+cos(k*deltaCrossDegrees)*
			(sin(k*deltaAlongDegrees)*aheadOfTarget
					+ cos(k*deltaAlongDegrees)*target));
}

vector<SpatialVector> plotCircularTrack(
		VizHTM *viz,
		float64 latDegrees, float64 lonDegrees,
		float64 groundTrackDegreesFromEast,
		float r, float g, float b, float alpha,
		float64 thetaAlongTrackDegrees0, float64 thetaAlongTrackDegrees1,
		float64 phiCrossTrackDegrees0, float64 phiCrossTrackDegrees1,
		float64 deltaAlongTrackDegrees, float64 deltaCrossTrackDegrees,
		float64 EarthRotationDegreePerAlong360Degrees=0.0
		) {

	if( EarthRotationDegreePerAlong360Degrees != 0.0 ) {
		cout << "plotCircularTrack::Warning! EarthRotationAndPrecessionNotImplemented!!!" << endl << flush;
	}

	vector<SpatialVector> granules;

	SpatialVector target; target.setLatLonDegrees(latDegrees,lonDegrees);
	const SpatialVector north = SpatialVector(0.,0.,1.);
	SpatialVector dTargetHat, L;
	if(abs(north*target)!=1) {
		SpatialVector east  = north^target; east.normalize();
		dTargetHat =
				cos(groundTrackDegreesFromEast)*east +
				sin(groundTrackDegreesFromEast)*north;
		L = target^dTargetHat;
	} else {
		cout << "plotCircularTrack::UnimplementedCase: target parallel to north. Please fix." << endl << flush;
		cout << "target: " << target << endl << flush;
		cout << "north*target: " << north*target << endl << flush;
		exit(1);
	}

	SpatialVector aheadOfTarget = L^target; // Should be a unit vector.

	SpatialVector origin(0.,0.,0.);
//	viz->triaxis();
	testTenDegreeGridRGB(viz,0.3,0.3,0.3);
//	viz->addEdge(origin,target,0.5,1.,0.5);
//	viz->addEdge(origin,aheadOfTarget,0.5,1.,0.5);
//	viz->addEdge(origin,L,0.5,1.,1.0);

//	float64 k = 2.*M_PI/360.0;

	const float64 RadiansPerDegree = 2.*M_PI/360.0;
	const float64 shiftDuringOrbit0 = (EarthRotationDegreePerAlong360Degrees/360.0)*thetaAlongTrackDegrees0*RadiansPerDegree;
	SpatialVector target0 = target.rotatedAbout(north,shiftDuringOrbit0); target0.normalize();
	SpatialVector aheadOfTarget0 = L^target0; aheadOfTarget0.normalize();
	SpatialVector* x0 =
			pointFromReferenceAndDeltas(
					target0,
					aheadOfTarget0,
					thetaAlongTrackDegrees0,
					0.);

	const float64 shiftDuringOrbit1 = (EarthRotationDegreePerAlong360Degrees/360.0)*thetaAlongTrackDegrees1*RadiansPerDegree;
	SpatialVector target1 = target.rotatedAbout(north,shiftDuringOrbit1); target1.normalize();
	SpatialVector aheadOfTarget1 = L^target1; aheadOfTarget1.normalize();
	SpatialVector* x1 =
			pointFromReferenceAndDeltas(
					target1,
					aheadOfTarget1,
					thetaAlongTrackDegrees1,
					0.);

//	cout << "s0: " << shiftDuringOrbit0 << endl << flush;
//	cout << "s1: " << shiftDuringOrbit1 << endl << flush;
//	cout << "t0: " << target0 << endl << flush;
//	cout << "t:  " << target << endl << flush;
//	cout << "t1: " << target1 << endl << flush;
//	cout << "x0: " << *x0 << endl << flush;
//	cout << "x1: " << *x1 << endl << flush;
	viz->addArc(*x0,target,1,1,1,-1.0,-1.0,200);
	viz->addArc(target,*x1,1,1,1,-1.0,-1.0,200);

	int steps = (thetaAlongTrackDegrees1-thetaAlongTrackDegrees0)/deltaAlongTrackDegrees;
	int crossSteps = (phiCrossTrackDegrees1-phiCrossTrackDegrees0)/deltaCrossTrackDegrees;

	if (steps < 1) steps = 1;
	if (crossSteps < 1) crossSteps = 1;
//	steps=1;
//	crossSteps=1;

	SpatialVector *v0, *v1, *v2, *v3;
	for(int isteps=0;isteps<steps;isteps++) {
		for(int csteps=0;csteps<crossSteps;csteps++) {
			v0 = pointFromReferenceAndDeltas(target,aheadOfTarget,thetaAlongTrackDegrees0+isteps*deltaAlongTrackDegrees,phiCrossTrackDegrees0+csteps*deltaCrossTrackDegrees);
			v1 = pointFromReferenceAndDeltas(target,aheadOfTarget,thetaAlongTrackDegrees0+(isteps+1)*deltaAlongTrackDegrees,phiCrossTrackDegrees0+csteps*deltaCrossTrackDegrees);
			v2 = pointFromReferenceAndDeltas(target,aheadOfTarget,thetaAlongTrackDegrees0+(isteps+1)*deltaAlongTrackDegrees,phiCrossTrackDegrees0+(csteps+1)*deltaCrossTrackDegrees);
			v3 = pointFromReferenceAndDeltas(target,aheadOfTarget,thetaAlongTrackDegrees0+isteps*deltaAlongTrackDegrees,phiCrossTrackDegrees0+(csteps+1)*deltaCrossTrackDegrees);
			viz->addArc(*v0,*v1,r,g,b,alpha);
			viz->addArc(*v1,*v2,r,g,b,alpha);
			viz->addArc(*v2,*v3,r,g,b,alpha);
			viz->addArc(*v3,*v0,r,g,b,alpha);
			granules.push_back(*v0);
			granules.push_back(*v1);
			granules.push_back(*v2);
			granules.push_back(*v3);
		}
	}
	return granules;
}


void testPlotDataSetIntersection(VizHTM *viz) {
	SpatialIndex index(5,5);
	SpatialVector u0, u1, u2, u3, v0, v1, v2, v3;

	HtmRange *rangeU = new HtmRange;
	HtmRange *rangeV = new HtmRange;
	HtmRange *rangeIntersect = new HtmRange;

	rangeU->purge(); rangeV->purge(); rangeIntersect->purge();

	if(false){
		u0.setLatLonDegrees(10.0,0.0);
		u1.setLatLonDegrees(10.0,20.0);
		u2.setLatLonDegrees(0.0,20.0);
		u3.setLatLonDegrees(0.0,0.0);

		v0.setLatLonDegrees(5.0,10.0);
		v1.setLatLonDegrees(5.0,40.0);
		v2.setLatLonDegrees(18.0,40.0);
		v3.setLatLonDegrees(18.0,0.0);

		intersectTwoRectangles(
				viz,
				&index,
				&u0,&u1,&u2,&u3,
				&v0,&v1,&v2,&v3,
				rangeU,
				rangeV,
				rangeIntersect
				);
		viz->addHTMRange(&index,rangeU,1.0,0.0,0.0,0.7);
		viz->addHTMRange(&index,rangeV,0.0,1.0,0.0,0.7);
		viz->addHTMRange(&index,rangeIntersect,0.0,0.5,1.0,0.2);

		delete rangeU, rangeV, rangeIntersect;
		return;
	}

	const SpatialVector zHat = SpatialVector(0.,0.,1.);
	for(int i=0; i<10; i++){
		cout << "i: " << i << endl << flush;

		double delta = 10.;
		double lat = -80. + 160.0*uniformDouble();
		double lon = 360.0 * uniformDouble();
		cout << " udll: " << delta << " " << lat << " " << lon << endl << flush;
		u0.setLatLonDegrees(lat,lon);
		u1.setLatLonDegrees(lat+delta,lon);
		u2.setLatLonDegrees(lat+delta,lon+delta);
		u3.setLatLonDegrees(lat,lon+delta);

		double dlat = 12.0*(uniformDouble()-0.5);
		double dlon = 12.0*(uniformDouble()-0.5);
		cout << "dlat,dlon: " << dlat << " " << dlon << endl << flush;
		double theta = 90.0*uniformDouble();
//		lat = uniformDouble(-80.,80.);
//		lon = uniformDouble(10.,350.);
		lat+=dlat; lon+=dlon;
		cout << " vdll: " << delta << " " << lat << " " << lon << endl << flush;
		v0.setLatLonDegrees(lat,lon+sin(theta)*delta);
		v1.setLatLonDegrees(lat+cos(theta)*delta,lon+sin(theta)*delta);
		v2.setLatLonDegrees(lat+sin(theta)*delta,lon+cos(theta)*delta);
		v3.setLatLonDegrees(lat+sin(theta)*delta,lon);

		if(i==i){
//		if(i==2) {
//		if(false){
//		if(i==14){
			cout << "intersect " << endl << flush;
			intersectTwoRectangles(
					viz,
					&index,
					&u0,&u1,&u2,&u3,
					&v0,&v1,&v2,&v3,
					rangeU,
					rangeV,
					rangeIntersect
			);
		}
//		if(i!=i){
		if(i==i){
//		if(i==2){
			if(true){
				float r = 0.0;
				float g = 0.0;
				float b = 1.0;
				cout << "blue" << endl << flush;
				viz->addRectangle(u0,u1,u2,u3,r,g,b);
			}
			if(true){
				float r = 1.0;
				float g = 0.0;
				float b = 0.0;
				cout << "red" << endl << flush;
				viz->addRectangle(v0,v1,v2,v3,r,g,b);
			}
		}
		cout << "+" << endl << flush;
	}

	cout << "add ranges" << endl << flush;

	viz->addHTMRange(&index,rangeU,1.0,0.0,0.0,0.7);
	viz->addHTMRange(&index,rangeV,0.0,1.0,0.0,0.7);
	viz->addHTMRange(&index,rangeIntersect,0.0,0.5,1.0,0.2);

	delete rangeU, rangeV, rangeIntersect;
}


array<vector<SpatialVector>,2> DataIntersectionDriver1(
		VizHTM *viz) {
	vector<SpatialVector> granules0, granules1;

	{
		float64 latDegrees = 0.3;
		float64 lonDegrees = 45.7;
		float64 groundTrackDegreesFromEast = 45.0;
		float r = 0.8;
		float g = 0.3;
		float b = 1.0;
		float a = -1;
		float64 thetaAlongTrackDegrees0 = -30;
		float64 thetaAlongTrackDegrees1 =  -5;
		float64 phiCrossTrackDegrees0 =   -5;
		float64 phiCrossTrackDegrees1 =    0;
		float64 deltaAlongTrackDegrees = 3.5;
		float64 deltaCrossTrackDegrees = 8;

		granules0 = plotCircularTrack(
				viz,
				latDegrees,  lonDegrees,
				groundTrackDegreesFromEast,
				r,  g,  b,  a,
				thetaAlongTrackDegrees0,  thetaAlongTrackDegrees1,
				phiCrossTrackDegrees0,  phiCrossTrackDegrees1,
				deltaAlongTrackDegrees,  deltaCrossTrackDegrees
		);
	}

	{
		float64 latDegrees = 0.7;
		float64 lonDegrees = 45.3;
		float64 groundTrackDegreesFromEast = 5.0;
		float r = 0.0;
		float g = 0.8;
		float b = 0.8;
		float a = -1;
		float64 thetaAlongTrackDegrees0 = 7.5;
		float64 thetaAlongTrackDegrees1 =  10;
		float64 phiCrossTrackDegrees0 =   -5;
		float64 phiCrossTrackDegrees1 =    -2.5;
		float64 deltaAlongTrackDegrees = 2.5;
		float64 deltaCrossTrackDegrees = 2.5;

		granules1 = plotCircularTrack(
				viz,
				latDegrees,  lonDegrees,
				groundTrackDegreesFromEast,
				r,  g,  b,  a,
				thetaAlongTrackDegrees0,  thetaAlongTrackDegrees1,
				phiCrossTrackDegrees0,  phiCrossTrackDegrees1,
				deltaAlongTrackDegrees,  deltaCrossTrackDegrees
		);
	}
	return {granules0,granules1};
}

array<vector<SpatialVector>,2> DataIntersectionDriver(
		VizHTM *viz) {
	vector<SpatialVector> granules0, granules1;

	{
		float64 latDegrees = 0.3;
		float64 lonDegrees = 45.7;
		float64 groundTrackDegreesFromEast = 45.0;
		float r = 0.8;
		float g = 0.3;
		float b = 1.0;
		float a = -1;
		float64 thetaAlongTrackDegrees0 = -30;
		float64 thetaAlongTrackDegrees1 =  30;
		float64 phiCrossTrackDegrees0 =   -4;
		float64 phiCrossTrackDegrees1 =    4;
		float64 deltaAlongTrackDegrees = 4;
		float64 deltaCrossTrackDegrees = 4;

		granules0 = plotCircularTrack(
				viz,
				latDegrees,  lonDegrees,
				groundTrackDegreesFromEast,
				r,  g,  b,  a,
				thetaAlongTrackDegrees0,  thetaAlongTrackDegrees1,
				phiCrossTrackDegrees0,  phiCrossTrackDegrees1,
				deltaAlongTrackDegrees,  deltaCrossTrackDegrees
		);
	}

	{
		float64 latDegrees = 0.7;
		float64 lonDegrees = 45.3;
		float64 groundTrackDegreesFromEast = 5.0;
		float r = 0.0;
		float g = 0.8;
		float b = 0.8;
		float a = -1;
		float64 thetaAlongTrackDegrees0 = -40;
		float64 thetaAlongTrackDegrees1 =  40;
		float64 phiCrossTrackDegrees0 =   -5;
		float64 phiCrossTrackDegrees1 =    5;
		float64 deltaAlongTrackDegrees = 2.5;
		float64 deltaCrossTrackDegrees = 2.5;

		granules1 = plotCircularTrack(
				viz,
				latDegrees,  lonDegrees,
				groundTrackDegreesFromEast,
				r,  g,  b,  a,
				thetaAlongTrackDegrees0,  thetaAlongTrackDegrees1,
				phiCrossTrackDegrees0,  phiCrossTrackDegrees1,
				deltaAlongTrackDegrees,  deltaCrossTrackDegrees
		);
	}
	return {granules0,granules1};
}



void testPlotDataSetIntersection0_PlotHtmRangeContains(VizHTM *viz, uint htmIdLevel=4, uint saveLevel=5) {
	SpatialIndex index(htmIdLevel,saveLevel);
	viz->lineWidth = 1.5;
	viz->sphereComplexity = 1.0;
	if(htmIdLevel>2) {
		plotBlockingSphere(viz,0.3,0.3,0.3,0.99);
	}

	array<vector<SpatialVector>,2> granuleSets;
	vector<SpatialVector> *granules0, *granules1;

	// Get the data
	granuleSets = DataIntersectionDriver(viz);
	granules0 = &(granuleSets[0]);
	granules1 = &(granuleSets[1]);
	// Granules now contain sequences of 4-vertex lists corresponding to geometries to be compared.

//	cout << "100" << endl << flush;
//	vector<SpatialVector> *tmp;
//	tmp = granules0;
//	granules0 = granules1;
//	granules1 = tmp;

	HtmRange *rangeU = new HtmRange;
	HtmRange *rangeV = new HtmRange;
	HtmRange *rangeIntersect = new HtmRange;
	rangeU->purge();
	rangeV->purge();
	rangeIntersect->purge();

//	cout << "200" << endl << flush;

	bool focus = false; // Focus means restrict attention to a particular box pair.
	int iFocus = 0;
	int jFocus = 0;
//	SpatialIndex index(3,5);

//	cout << "300" << endl << flush;

//	cout << "intersecting" << endl << flush;
	int count=0;
	int i=0, j=0;
//	cout << "uv<" << (granules0->size()/4)
//			<< "," << (granules1->size()/4)
//			<< ">" << endl << flush;
	for(
			vector<SpatialVector>::iterator iterU = granules0->begin();
			iterU != granules0->end();
			iterU++, i++) {
//		cout<< "u" << flush;
		// Extract a box from the u-list.
		SpatialVector u0 = *iterU++;
		SpatialVector u1 = *iterU++;
		SpatialVector u2 = *iterU++;
		SpatialVector u3 = *iterU;
		j=0;
		for(
				vector<SpatialVector>::iterator iterV = granules1->begin();
				iterV != granules1->end();
				iterV++, j++) {
//			cout << "v" << flush;
//			cout << "<ij= " << i << " " << j << " >" << flush;

			// Extract a box from the v-list.
			SpatialVector v0 = *iterV++;
			SpatialVector v1 = *iterV++;
			SpatialVector v2 = *iterV++;
			SpatialVector v3 = *iterV;

//			cout << "." << flush;
			if(!focus){ // The general case.
				intersectTwoRectangles(
						viz,
						&index,
						&u0,&u1,&u2,&u3,
						&v0,&v1,&v2,&v3,
						rangeU,
						rangeV,
						rangeIntersect
				);
			} else { // The special case for debugging, restricted to a particular pair.
				if((iFocus==i)&&(jFocus==j)){
					intersectTwoRectangles(  // Visualizes the HtmRange-based comparison (on triangles).
							viz,
							&index,
							&u0,&u1,&u2,&u3,
							&v0,&v1,&v2,&v3,
							rangeU,
							rangeV,
							rangeIntersect
					);
					viz->addRectangle(u0,u1,u2,u3,1.0,0.,0.); // Visualize the raw boxes.
					viz->addRectangle(v0,v1,v2,v3,0.0,0.,1.);
				}
			}
			count++;

		}
	}

	rangeU->defrag(); rangeV->defrag(); // TODO Defrag led to a little strangeness. May require double-checking.

	Key lo=-1, hi=-1;
	HtmRange red, blue, green;
	HtmRange *range1, *range2;
	range1=rangeU;
	range2=rangeV;
	// range2=rangeIntersect;
	range1->reset();
	int indexp = range1->getNext(lo,hi);
	if(indexp) {
		do {
			//			cout << "r2.contains: " << lo << " " << hi << " " << range2.contains(lo,hi) << endl << flush;
			int ret = range2->contains(lo,hi);
//			cout << "r2.cont: lo,hi,ret: " << lo << " " << " " << hi << " " << ret << endl << flush;
			if(ret == -1) {
				red.addRange(lo,hi);
			} else if (ret == 0 ) {
				green.addRange(lo,hi);
			} else {
				blue.addRange(lo,hi);
			}
		} while(range1->getNext(lo,hi));
		viz->addHTMRange(&index,&red,0.7,0.,0.);
		viz->addHTMRange(&index,&green,0.,0.7,0.);
		viz->addHTMRange(&index,&blue,0.,0.,0.7);
		//			cout << "red" << endl;
		//			red->print(HtmRange::BOTH,cout);
		//			cout << "green" << endl;
		//			green->print(HtmRange::BOTH,cout);
		//			cout << "blue" << endl;
		//			blue->print(HtmRange::BOTH, cout);
	}

//	cout << "r1: " << endl << flush;
//	range1->print(HtmRange::BOTH,cout,false);
//	cout << "r2: " << endl << flush;
//	range2->print(HtmRange::BOTH,cout,false);
//
//	cout << " r1.contains(141,141) == true : " << range1->contains(141,141) << endl << flush;
//	cout << " r1.contains(143,143) == true : " << range1->contains(143,143) << endl << flush;
//	cout << " r1.contains(130,131) == true : " << range1->contains(130,131) << endl << flush;

	delete rangeU, rangeV, rangeIntersect;

}



void testPlotDataSetIntersection0_PlotIntersectTwoRectanglesOutput(VizHTM *viz, uint htmIdLevel=4, uint saveLevel=5) {
	SpatialIndex index(htmIdLevel,saveLevel);

	viz->lineWidth = 1.5;
	viz->sphereComplexity = 1.0;
	if(htmIdLevel>2) {
		plotBlockingSphere(viz,0.3,0.3,0.3,0.99);
	}

	array<vector<SpatialVector>,2> granuleSets;
	vector<SpatialVector> *granules0, *granules1;

	granuleSets = DataIntersectionDriver(viz);
	granules0 = &(granuleSets[0]);
	granules1 = &(granuleSets[1]);

//	vector<SpatialVector> *tmp;
//	tmp = granules0;
//	granules0 = granules1;
//	granules1 = tmp;

	HtmRange *rangeU = new HtmRange;
	HtmRange *rangeV = new HtmRange;
	HtmRange *rangeIntersect = new HtmRange;
	rangeU->purge();
	rangeV->purge();
	rangeIntersect->purge();

	bool focus = false;
	int iFocus = 0;
	int jFocus = 0;
//	cout << "intersecting" << endl << flush;
	int count=0;
	int i=0, j=0;
//	cout << "uv<" << (granules0->size()/4)
//			<< "," << (granules1->size()/4)
//			<< ">" << flush;
	for(
			vector<SpatialVector>::iterator iterU = granules0->begin();
			iterU != granules0->end();
			iterU++, i++) {
//		cout<< "u" << flush;
		SpatialVector u0 = *iterU++;
		SpatialVector u1 = *iterU++;
		SpatialVector u2 = *iterU++;
		SpatialVector u3 = *iterU;
		j=0;
		for(
				vector<SpatialVector>::iterator iterV = granules1->begin();
				iterV != granules1->end();
				iterV++, j++) {
//			cout << "v" << flush;
//			cout << "<ij= " << i << " " << j << " >" << flush;

			SpatialVector v0 = *iterV++;
			SpatialVector v1 = *iterV++;
			SpatialVector v2 = *iterV++;
			SpatialVector v3 = *iterV;

//			cout << "." << flush;
			if(!focus){
				intersectTwoRectangles(
						viz,
						&index,
						&u0,&u1,&u2,&u3,
						&v0,&v1,&v2,&v3,
						rangeU,
						rangeV,
						rangeIntersect
				);
			} else {
				if((iFocus==i)&&(jFocus==j)){
					intersectTwoRectangles(
							viz,
							&index,
							&u0,&u1,&u2,&u3,
							&v0,&v1,&v2,&v3,
							rangeU,
							rangeV,
							rangeIntersect
					);
					viz->addRectangle(u0,u1,u2,u3,1.0,0.,0.);
					viz->addRectangle(v0,v1,v2,v3,0.0,0.,1.);
				}
			}
			count++;
		}
	}
	viz->addHTMRange(&index,rangeIntersect,0.9,0.1,0.9,0.);
	viz->addHTMRange(&index,rangeU,0.5,0.9,0.1,0.);
	viz->addHTMRange(&index,rangeV,0.1,0.9,0.5,0.);

	delete rangeU, rangeV, rangeIntersect;
}

void testMultiResolutionHtmRange(VizHTM *viz) {

	SpatialIndex index;
	uint64 htmId;
	HtmRange *r = new HtmRange; // Use bitshifted format

	htmId = index.idByName("N01");
	r->addRange(htmId,htmId);
	viz->addHTMRange(r,0.1,0.8,0.8,0.0);

	htmId = index.idByName("N0133");
	r->purge();
	r->addRange(htmId,htmId);
	viz->addHTMRange(r,0.8,0.1,0.1,0.0);

	htmId = index.idByName("N013333");
	r->purge();
	r->addRange(htmId,htmId);
	viz->addHTMRange(r,0.8,0.1,0.8,0.0);

	r->purge();
	htmId = index.idByName("N21");
	r->addRange(htmId,htmId);
	htmId = index.idByName("N2133");
	r->addRange(htmId-3,htmId);
	htmId = index.idByName("N213333");
	r->addRange(htmId-3,htmId);
	viz->addHTMRange(r,0.1,0.6,0.9,0.0);

}


void printHtmIdInfo(uint64 id){
	BitShiftNameEncoding rightJustified;
	EmbeddedLevelNameEncoding leftJustified;
	string name;
	name = rightJustified.nameById(id);
	cout << "HtmIdInfo " << id << " ";
	cout << hex << id << dec << " " ;
	cout << rightJustified.nameById(id) << " ";
	cout << leftJustified.idByName(rightJustified.nameById(id)) << " ";
	cout << hex << leftJustified.idByName(rightJustified.nameById(id)) << dec << " ";

	leftJustified.setName(rightJustified.nameById(id));
	cout << leftJustified.getId_NoEmbeddedLevel() << " ";
	cout << hex << leftJustified.getId_NoEmbeddedLevel() << dec << " ";
	cout << endl << flush;
}

void testChunks(VizHTM *viz) {
	cout << "testChunks start" << endl << flush;
	BitShiftNameEncoding rightJustified;
	EmbeddedLevelNameEncoding leftJustified;
	int64 idDelta = 100000000000000000;
	//              100000000000000000
	int64 idMax  = 9000000000000000000;
	int level = 1;
	SpatialIndex *index = new SpatialIndex(level);
	HtmRange *r = new HtmRange;
	for(int64 id = 0; id < idMax; id += idDelta) {
		leftJustified.setIdFromSciDBLeftJustifiedFormat(id+level);
		cout << id+level << " zero-lj: " << flush;
		cout << leftJustified.getName() << endl << flush;
		rightJustified.setName(leftJustified.getName());
		r->purge();
		int64 id_ = rightJustified.getId();
//		r->addRange(id_,id_);
//		r->reset();
//		viz->addHTMRange(r,0.8,0.2,0.5,0.0);
		viz->addArcFromIndexAndId(index,id_,0.8,0.5,0.5,0.0);
//		r->reset();
	}
	cout << "testChunks done" << endl << flush;
}

void testLevelChunk(VizHTM *viz) {

	bool referenceGrid = false;
	bool level3 = false;
	bool level4 = false;
	bool faces  = true;
	bool level5 = true;
	bool printInfo = true;
	viz->lineWidth = 1.5;

	if(referenceGrid) {
		testTenDegreeGridRGB(viz,0.6,0.6,0.6);
		viz->sphereComplexity = 1.0;
//		if(htmIdLevel>2) {
			plotBlockingSphere(viz,0.3,0.3,0.3,0.99);
//		}
	}

	SpatialIndex *index3 = new SpatialIndex(3,5);
	SpatialIndex *index4 = new SpatialIndex(4,5);
	SpatialIndex *index5 = new SpatialIndex(5,5);
	uint64 htmId;
	HtmRange *r = new HtmRange;

	if(level3){
		r->purge();
		r->addRange(1023,1023);
		viz->addHTMRange(r,0.1,0.8,0.1,0.0);
		//	viz->addArcFromIndexAndId(index3,1023,0.1,0.8,0.1,0.0);
	}

	if(level4){
		r->purge();
		r->addRange(4092,4095);
		//	r->addRange(4092,4092);
		viz->addHTMRange(r,0.8,0.1,0.8,0.0);
		//	viz->addArcFromIndexAndId(index4,4092,0.1,0.8,0.8,0.0);
	}

	if(faces){
		float intensity = 0.5;
		viz->addFaceFromIndexAndId(
				index4,4092,
				intensity, 0, 0,
				intensity, 0, 0,
				intensity, 0, 0
		);
		viz->addFaceFromIndexAndId(
				index4,4093,
				0, intensity, 0,
				0, intensity, 0,
				0, intensity, 0
		);
		viz->addFaceFromIndexAndId(
				index4,4094,
				0, 0, intensity,
				0, 0, intensity,
				0, 0, intensity
		);
		viz->addFaceFromIndexAndId(
				index4,4095,
				intensity, intensity, intensity,
				intensity, intensity, intensity,
				intensity, intensity, intensity
		);
	}

	if(level5){
		//	r->purge();
		////	r->addRange(16368,16383);
		//	r->addRange(16368,16368);
		//	viz->addHTMRange(r,0.1,0.1,0.8,0.0);
		for(int i=0;i<16;i++){
			viz->addArcFromIndexAndId(index5,16368+i,0.1,0.8,0.8,0.0);
		}
	}

	if(printInfo) {
		printHtmIdInfo(1023);
		printHtmIdInfo(4092);
		printHtmIdInfo(4093);
		printHtmIdInfo(4094);
		printHtmIdInfo(4095);
		printHtmIdInfo(16368);
		printHtmIdInfo(16368+3);
		printHtmIdInfo(16372);
		printHtmIdInfo(16372+3);
		printHtmIdInfo(16376);
		printHtmIdInfo(16376+3);
		printHtmIdInfo(16380);
		printHtmIdInfo(16383);
	}


	delete index3, index4, index5, r;
}


void plotHTMInterval(VizHTM *viz, SpatialIndex index, htmRange interval) {
	size_t htmIdLevel = index.getMaxlevel();
	Key lo = interval.lo;
	Key hi = interval.hi;
	SpatialVector x1,x2,x3;
	KeyPair adjustedRange = HTMRangeAtLevelFromHTMRange(htmIdLevel,lo,hi);
	lo = adjustedRange.lo;
	hi = adjustedRange.hi;

	for(uint64 numericId=lo; numericId<=hi;numericId++) {
		uint64 nodeIndex = index.nodeIndexFromId(numericId);
		if(nodeIndex!=0){
			index.nodeVertex(nodeIndex,x1,x2,x3);
			if(true) {
				float r=0.; float g=0.; float b=0.;
				switch(numericId % 4) {
				case 0:
					r=1.;
					break;
				case 1:
					g=1.;
					break;
				case 2:
					b=1.;
					break;
				default:
					r=1.; g=1.; b=1.;
					break;
				}
				//	for(int i=0; i<3; i++) addEdgeColor(r,g,b);
				int colorBase = viz->nFaceColors;
				for(int i=0; i<3; i++) {
					viz->addFaceColor(r,g,b);
				}
				viz->addFaceVertexColorIndices3(colorBase,colorBase+1,colorBase+2);

				int indexBase = viz->nCoordinates;
				viz->addCoordinate64(x1.x(),x1.y(),x1.z());
				viz->addCoordinate64(x2.x(),x2.y(),x2.z());
				viz->addCoordinate64(x3.x(),x3.y(),x3.z());
				viz->addFaceIndices3(indexBase,indexBase+1,indexBase+2);

				//					printf("id: %llx ix: %llu",numericId,nodeIndex);
				//					cout << endl << flush;

				if(true){
					float size = pow(0.5,htmIdLevel+3);
					float r = 0.4, g = 0.4, b = 0.6;
					SpatialVector x = 3.*x1+x2+x3; x.normalize(); x *= 1.0+1.0e-6;
					SpatialVector x_ = x1+x2+x3; x_.normalize();
//					if(false){
//						cout
//						<< " nI: " << nodeIndex
//						<< " x: " << x.x() << " " << x.y() << " " << x.z();
//					}
					char *str = new char[256];
					sprintf(str,"id: %llx\nix: %llu\n",numericId,nodeIndex);
					viz->addAnnotation((new SpatialVector(x)),str,size,r,g,b);
					viz->addEdge(x,x_,0.5,0.5,0.9);
				}
			}
		}
	}
}

void testHTMRange(VizHTM *viz, int htmIdLevel, const char *n0, const char *n1) {
	int buildLevel = 5;
	SpatialIndex *index = new SpatialIndex(htmIdLevel,buildLevel);
	BitShiftNameEncoding name0 = BitShiftNameEncoding(n0);
	BitShiftNameEncoding name1 = BitShiftNameEncoding(n1);

	htmRange interval;
	interval.lo = name0.getId();
	interval.hi = name1.getId();

	plotHTMInterval(viz,*index,interval);
}

void testTwoRectangle(VizHTM *viz, int htmIdLevel) {
	int saveLevel = 5;

	SpatialIndex *index = new SpatialIndex(htmIdLevel,saveLevel);

	SpatialDomain domain1 = SpatialDomain(index);
//	domain1.setOlevel(htmIdLevel); // Note this sets the olevel on the convexes.

	SpatialDomain domain2 = SpatialDomain(index);
//	domain2.setOlevel(htmIdLevel); // Note this sets the olevel on the convexes.

	if(true){
		SpatialVector *v0 = VectorFromLatLonDegrees(10.0,0.0);
		SpatialVector *v1 = VectorFromLatLonDegrees(30.0,-10.0);
		SpatialVector *v2 = VectorFromLatLonDegrees(60.0,10.0);
		SpatialVector *v3 = VectorFromLatLonDegrees(40.0,20.0);
		float r = 1.0;
		float g = 0.0;
		float b = 1.0;

		RangeConvex *rc = new RangeConvex(v0,v1,v2,v3);
		//	cout << "nConstraints: " << rc->numConstraints() << endl << flush;
		//	SpatialConstraint *sc = new SpatialConstraint(SpatialVector(0.,0.,1.),0.5);
		//	viz->addConstraint(*sc,1.0,1.0,1.0);
		//	rc->add(*sc);
//		rc->setOlevel(htmIdLevel); // Note this is supposed to be done when added to the domain.  ///???
		viz->addRectangle(*v0,*v1,*v2,*v3,r,g,b);
	}

	if(true){
		SpatialVector *v0 = VectorFromLatLonDegrees(30.0,20.0);
		SpatialVector *v1 = VectorFromLatLonDegrees(30.0,-10.0);
		SpatialVector *v2 = VectorFromLatLonDegrees(20.0,0.0);
		SpatialVector *v3 = VectorFromLatLonDegrees(20.0,30.0);
		float r = 1.0;
		float g = 1.0;
		float b = 0.0;

		RangeConvex *rc = new RangeConvex(v0,v1,v2,v3);
		//	cout << "nConstraints: " << rc->numConstraints() << endl << flush;
		//	SpatialConstraint *sc = new SpatialConstraint(SpatialVector(0.,0.,1.),0.5);
		//	viz->addConstraint(*sc,1.0,1.0,1.0);
		//	rc->add(*sc);
//		rc->setOlevel(htmIdLevel); // Note this is supposed to be done when added to the domain. ///???
		domain2.add(*rc);
		viz->addRectangle(*v0,*v1,*v2,*v3,r,g,b);
	}

	//	return;

	bool varlen_individualHTMIds = false; // true for individuals, false for ranges

	HtmRange *range1 = new HtmRange();
	range1->purge();
	bool overlap1 = domain1.intersect(index,range1,varlen_individualHTMIds);
	range1->defrag();
	range1->reset();

	HtmRange *range2 = new HtmRange();
	range2->purge();
	bool overlap2 = domain2.intersect(index,range2,varlen_individualHTMIds);
	range2->defrag();
	range2->reset();

	Key lo = 0, hi = 0;
	uint64 indexp = 0;

	HtmRange *resultRange = HTMRangeAtLevelFromIntersection(htmIdLevel,range1,range2);

	HtmRange *range = resultRange;
	range->reset();
	indexp = range->getNext(lo,hi);
	SpatialVector x1,x2,x3;
	//	if(indexp)
	do {
		for(uint64 numericId=lo; numericId<=hi;numericId++) {
			uint64 nodeIndex = index->nodeIndexFromId(numericId);
			if(nodeIndex!=0){
				index->nodeVertex(nodeIndex,x1,x2,x3);
				if(true) {
					float r=0.; float g=0.; float b=0.;
					switch(numericId % 4) {
					case 0:
						r=1.;
						break;
					case 1:
						g=1.;
						break;
					case 2:
						b=1.;
						break;
					default:
						r=1.; g=1.; b=1.;
						break;
					}
					//	for(int i=0; i<3; i++) addEdgeColor(r,g,b);
					int colorBase = viz->nFaceColors;
					for(int i=0; i<3; i++) {
						viz->addFaceColor(r,g,b);
					}
					viz->addFaceVertexColorIndices3(colorBase,colorBase+1,colorBase+2);

					int indexBase = viz->nCoordinates;
					viz->addCoordinate64(x1.x(),x1.y(),x1.z());
					viz->addCoordinate64(x2.x(),x2.y(),x2.z());
					viz->addCoordinate64(x3.x(),x3.y(),x3.z());
					viz->addFaceIndices3(indexBase,indexBase+1,indexBase+2);

					//					printf("id: %llx ix: %llu",numericId,nodeIndex);
					//					cout << endl << flush;

					if(true){
						float size = pow(0.5,htmIdLevel+3);
						float r = 0.4, g = 0.4, b = 0.6;
						SpatialVector x = 3.*x1+x2+x3; x.normalize(); x *= 1.0+1.0e-6;
						SpatialVector x_ = x1+x2+x3; x_.normalize();
						if(false){
							cout
							<< " nI: " << nodeIndex
							<< " x: " << x.x() << " " << x.y() << " " << x.z();
						}
						char *str = new char[256];
						sprintf(str,"id: %llx\nix: %llu\n",numericId,nodeIndex);
						viz->addAnnotation((new SpatialVector(x)),str,size,r,g,b);
						viz->addEdge(x,x_,0.5,0.5,0.9);
					}
				}
			}
		}
	} while (range->getNext(lo,hi));
}


void testTenDegreeGrid(VizHTM *viz) {
	int saveLevel = 5;
	float64 PI = atan2(0,-1);
	float64 k  = 2*PI/360.;

	float dLat = 10;
	float dLon = 10;
	for(float lat=-90; lat<=90; lat+= dLat) {
		for(float lon=0; lon<=360; lon+= dLon) {
			float r = 0.125+0.125*(1+cos(k*lon)); // r=0;
			float g = 0.125+0.125*(1+sin(k*lon)); // g=0;
			float b = 0.125+0.125*(1+sin(k*lat)); // b=0;
			viz->addLatLonBoxEdgesDegrees(
					lat + 0.25, lon + 0.25,
					lat + dLat - 0.25, lon + dLon - 0.25,
					r,g,b);
		}
	}
}


void testText1(VizHTM* viz, SpatialVector a) {
	viz->triaxis();
	if(true){
		cout << "Starting text test" << flush;
		SpatialVector o = SpatialVector(0.,0.,0.);
		viz->addEdge(o,a,1.,1.,1.);
		float size = 0.25;
		float r = 0.8, g = 0.8, b = 0.9;
		SpatialVector *x = new SpatialVector(a); x->normalize(); (*x) *= 1.01;
		char *str = new char[256];
		//		strcpy(str,"N0xx\n");
		sprintf(str,"%s\n%s","NOxy",a.toString());
		//		strcat(str,a.toString());
		viz->addAnnotation(x,str,size,r,g,b);
		cout << "'" << str << "'" << flush;
		cout << "...done." << endl << flush;
	}
	if(true){
		cout << "Starting text test 2" << flush;
		SpatialVector o = SpatialVector(0.,0.,0.);
		SpatialVector a = SpatialVector(0.,0.,1.);
		viz->addEdge(o,a,1.,1.,1.);
		float size = 0.25;
		float r = 0.8, g = 0.8, b = 0.9;
		SpatialVector *x = new SpatialVector(a); x->normalize(); (*x) *= 1.01;
		char *str = new char[256];
		//		strcpy(str,"N0xx\n");
		sprintf(str,"%s\n%s","NOxy",a.toString());
		//		strcat(str,a.toString());
		viz->addAnnotation(x,str,size,r,g,b);
		cout << "'" << str << "'" << flush;
		cout << "...done." << endl << flush;
	}
}

void testShowEdgeProjections(VizHTM *viz, SpatialVector a, float64 d) {
	SpatialVector o = SpatialVector(0.,0.,0.);
	if(true) {
		viz->triaxis();
		viz->addEdge(o,a,1.,1.,1.);
		viz->addConstraint(a,d,0.5,0.5,1.0);
		viz->addEdgeProjections(a*d);
	}
}

float lat_(int i, int nPoints) { return -90.0 + (i*180.)/(nPoints-1); }
float lon_(int i, int nPoints) { return -(i*16*4*90.)/(nPoints-1); }
float r_  (int i, int nPoints) { return (i/(1.0*nPoints)); }

void testLatLonSpiral(VizHTM* viz) { // Try to draw latlon lines.
	int indexBase = viz->nCoordinates;
	int edgeColorBase = viz->nEdgeColors;
	int nPoints = 1001;
	int nEdges = nPoints-1;
	float last_lat = lat_(0,nPoints);
	float last_lon = lon_(0,nPoints);
	//	float last_r   = r_  (0,nPoints);
	float64 *last_xyz= xyzFromLatLonDegrees(last_lat,last_lon);

	for(int i=1; i <= nEdges; i++) {
		float lat = lat_(i,nPoints);
		float lon = lon_(i,nPoints);
		float r   = r_  (i,nPoints);
		float g   = 1.0;
		float b   = 1.0;
		float64 *xyz= xyzFromLatLonDegrees(lat,lon);

		viz->addEdge(
				SpatialVector(last_xyz[0],last_xyz[1],last_xyz[2]),
				SpatialVector(xyz[0],xyz[1],xyz[2]),
				r, g, b);

		last_xyz = xyz;
		last_lat = lat;
		last_lon = lon;
		//		last_r   = r;

	}

	viz->addEdge(SpatialVector(0.0,0.0,1.0),SpatialVector(0.0,0.0,-1.0),0.9,0.9,1.0);

	//	for(int i = 0;i < nPoints; i++) {
	//		float lat = -90.0 + (i*180.)/(nPoints-1);	float lon = -(i*16*4*90.)/(nPoints-1);
	//		float* xyz = xyzFromLatLonDegrees(lat,lon);
	//		//		cout << "1000: " << i << " ( " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " ) " << endl;
	//		viz->addCoordinate(xyz[0],xyz[1],xyz[2]);
	//	}
	//	for(int i=0;i<nPoints;i++) {
	//		edgeColorBase = viz->nEdgeColors;
	//		float r = (i/(1.0*nPoints));
	//		float g = 1.0;
	//		float b = 1.0;
	//		viz->addEdgeColor(r,g,b);
	//		viz->addEdgeVertexColorIndices(edgeColorBase,edgeColorBase);
	//	}
	//	int nEdges = nPoints-1;
	//	for(int iPoint = 0; iPoint < nEdges; iPoint++) {
	//		viz->edgeIndices[viz->nEdgeIndices] = indexBase + iPoint; viz->nEdgeIndices++;
	//		viz->edgeIndices[viz->nEdgeIndices] = indexBase + iPoint + 1; viz->nEdgeIndices++;
	//		viz->edgeIndices[viz->nEdgeIndices] = SO_END_LINE_INDEX; viz->nEdgeIndices++;
	//	}

	//		cout << "1000: " << nCoordinates << " " << nEdgeColors << " " << nFaceColors << " " << nEdgeIndices << " " << nFaceIndices << endl << endl;
}

void testTriaxis(VizHTM* viz) {
	cout << "viz triaxis..." << flush;
	viz->triaxis();
	cout << "created." << endl << flush;
}
void testIJKRGBFace(VizHTM* viz) {
	cout << "viz test ijk-rgb-face..." << flush;
	viz->triaxis();
	viz->addFace3(
			1.0,0.0,0.0,
			0.0,1.0,0.0,
			0.0,0.0,1.0,
			1.0,0.0,0.0,
			0.0,1.0,0.0,
			0.0,0.0,1.0
	);
	cout << "done." << endl << flush;
}

//TODO Note that the convex hull is implemented in the HTM interface.
void testGeorgia(VizHTM* viz, int  level, int hullSteps,
		float hullR, float hullG, float hullB) {
//	cout << system("pwd") << endl << flush;
	string line;
	string fileName = "data/Georgia.Specification";
	ifstream specFile(fileName);
	float64 lat1, lon1;
	float64 lat0, lon0, latStart, lonStart;
	if(specFile.is_open()){
		cout << "Found file '" << fileName << "'. Reading..." << endl << flush;
		//		int count = 8;
		//		if(getline(specFile,line))
		//			while(specFile >> lat >> lon && count > 0) {
		//				cout << 8-count << " " << lat << " " << lon << endl << flush;
		//				count--;
		//			}
		LatLonDegrees64ValueVector latlon, latlon1;
		int count=0;
		int loadEvery = 1;
		if(getline(specFile,line)) {
			if(specFile >> latStart >> lonStart) {
				lat0 = latStart; lon0 = lonStart;
				latlon.push_back(LatLonDegrees64(lat0,lon0));
				SpatialVector *xStart = VectorFromLatLonDegrees(lat0,lon0);
				SpatialVector *x0 = xStart;
				count++;
				// Memory leak bait.
				while(specFile >> lat1 >> lon1) {
					//			cout << 8-count << " " << lat << " " << lon << endl << flush;
					if (!(count % loadEvery)) {
						latlon.push_back(LatLonDegrees64(lat1,lon1));
					}
					count++;
				}
//				for(int l=0; l<count; l++) latlon1.push_back(latlon[l]);
//				std::reverse(latlon.begin(),latlon.end());
//				int k=0;
				for(int l=latlon.size()-1; l>0; l--) {
					lat1 = latlon[l].lat;
					lon1 = latlon[l].lon;
					latlon1.push_back(LatLonDegrees64(lat1,lon1));
					SpatialVector *x1 = VectorFromLatLonDegrees(lat1,lon1);
					viz->addEdge((*x0),(*x1),0.2,0.5,0.9);
					x0 = x1;
				}
				latlon1.push_back(LatLonDegrees64(latStart,lonStart));
				viz->addEdge((*x0),(*xStart),0.2,0.5,0.9);
			}
		}
		specFile.close();
		cout << " latlon count: " << count << endl << flush;
		int buildLevel = 5;
//		cout << 1000 << endl << flush;
		htmInterface *htm = new htmInterface(level,buildLevel);
		cout << 1100 << endl << flush;
//		HTMRangeValueVector htmRangeVector = htm->convexHull(latlon1,hullSteps);
		HTMRangeValueVector htmRangeVector = htm->convexHull(latlon,hullSteps);

		SpatialVector *xStart = &(htm->polyCorners_[0].c_);
		SpatialVector *x0 = xStart;
		cout << " polyCorners: 0 " << flush;
		for(int j=1; j<htm->polyCorners_.size();j++){
			cout << j << " " << flush;
			SpatialVector *x1 = &(htm->polyCorners_[j].c_);
//			cout << *x0 << ", " << *x1 << "; " << flush;
			viz->addEdge((*x0)*1.001,(*x1)*1.001,hullR,hullG,hullB);
			x0=x1;
		}
		cout << endl << flush;
		viz->addEdge((*x0)*1.001,(*xStart)*1.001,hullR,hullG,hullB);

		cout << 1200 << endl << flush;
		cout << " htmRangeVector-size: " << htmRangeVector.size() << endl << flush;
		cout << " htmRangeVector-j:    ";
		for( int j=0; j<htmRangeVector.size(); j++) {
			htmRange hr = htmRangeVector[j];
			cout << " (j=" << j << " " << hr.lo << " " << hr.hi << " ) ";
			plotHTMInterval(viz,htm->index(),hr);
		}
		cout << endl << flush;
	} else {
		cout << "Couldn't open '" << fileName << "'." << endl << flush;
	}
	cout << "Leaving testGeorgia" << endl << flush;
}

void testAddEdgesFromIndexAndName(VizHTM* viz, const char* htmIdName, int saveLevel=5) {
	SpatialIndex *index = new SpatialIndex(htmIdName,saveLevel);
	viz->addEdgesFromIndexAndName(index,htmIdName,0.6,0.6,0.9);
}


string shapeFile = "data/ne_50m_coastline/ne_50m_coastline.shp";
float shapeFilesBlockingSphereScale = 0.999;

void testShapeFiles(VizHTM* viz) {
	bool verbose = false;
	if(true) {
		plotBlockingSphere(viz,0.2,0.2,0.2,shapeFilesBlockingSphereScale);
		testTenDegreeGridRGB(viz,0.6,0.6,0.6);
		//	string shapeFile = "data/ne_110m_coastline/ne_110m_coastline.shp"; // okay
		// string shapeFile = "data/ne_50m_coastline/ne_50m_coastline.shp"; // okay
		// string shapeFile = "data/ne_10m_coastline/ne_10m_coastline.shp"; // Doesn't work yet.
		float r = 0.5;
		float g = 0.9;
		float b = 0.9;
		viz->addShapeFile(shapeFile,r,g,b);
		cout << "testShape done" << endl << flush;
	}
}

void testTIFF(VizHTM *viz) {
	// http://www.libtiff.org/libtiff.html
	string tifFilename = "data/NE1_50M_SR_W/NE1_50M_SR_W.tif";
	QImage imageIn = QImage(tifFilename.c_str(),"tiff");
	QImage image = imageIn.mirrored(false,true);

	if (!image.isNull()) {
//	if(true) {
		uint32 w = image.width();  // lon
		uint32 h = image.height(); // lat

//		int w = 10800;
//		int h =  5400;

		cout << "w,h= " << w << ", " << h << endl;
		size_t npixels = w * h;

		float64 lat0 = -90.0;
		float64 dlat = 180.0/(h-1);
		float64 lon0 = -180.0;
		float64 dlon = 360.0/(w-1);

		// The whole planet
		int step = 2; // About the max memory...
		// step = 4; // for speed & fun
		// step = 100; // for speed & fun
		int i0 = 0; // w
		int j0 = 0; // h
		int i1 = w-step;
		int j1 = h-step;

		// Part of North America
		step = 1;
//		step = 2;
		i0 = 2000;
		j0 = 3000;
		int d1 = 2000;
//		d1 = 1000;
		i1 = i0 + d1;
		j1 = j0 + d1;

//		// English Channel
//		int i0 = 5400; // w
//		int j0 = 4200; // h
//		int d1 = 50;
//		int i1 = i0+d1;
//		int j1 = j0+d1;

		// TODO fix last pixel column

		viz->addLatLonBoxEdgesDegrees(
				lat0+j0*dlat,lon0+i0*dlon,
				lat0+j1*dlat,lon0+i1*dlon,
				0.0,0.8,0.0
				);

		// ...process raster data...
		for(int j=j0; j<j1; j += step) { // lat
			for(int i=i0; i<i1; i += step) { // lon
				QColor qc0, qc1, qc2, qc3;
				qc0 = QColor(image.pixel(i,j)); // x,y
				qc1 = QColor(image.pixel(i,j+step));
				qc2 = QColor(image.pixel(i+step,j+step));
				qc3 = QColor(image.pixel(i+step,j));

//				int ic = (i-i0)*255/d1;
//				int jc = (j-j0)*255/d1;
//				int dc = step*255/d1;
//				qc0 = QColor(ic,0,jc);
//				qc1 = QColor(ic,0,jc+dc);
//				qc2 = QColor(ic+dc,0,jc+dc);
//				qc3 = QColor(ic+dc,0,jc);

				viz->addFace4FromLatLonDegrees(
						lat0+j*dlat,        lon0+i*dlon,
						lat0+(j+step)*dlat, lon0+i*dlon,
						lat0+(j+step)*dlat, lon0+(i+step)*dlon,
						lat0+(j)*dlat,      lon0+(i+step)*dlon,
						qc0.redF(), qc0.greenF(), qc0.blueF(),
						qc1.redF(), qc1.greenF(), qc1.blueF(),
						qc2.redF(), qc2.greenF(), qc2.blueF(),
						qc3.redF(), qc3.greenF(), qc3.blueF()
						);



			}
		}
	}
}

void testDelaunay0(VizHTM *viz) {

	// Input
	int npt     = 5;
	int hashSize = 3*npt/2;
	cout << "Calling prime_c" << flush;
	int sizht   = prime_c(&hashSize);
	cout << ". " << sizht << endl << flush;
	int maxbf   = 1000;
	int maxfc   = 1000;
	int maxPt   = 100; // npt <= maxPt
	float64 vcl[3*maxPt];

	int i,j;
	i=0;j=0; vcl[j*3 + i] = 1.0;
	i=1;j=0; vcl[j*3 + i] = 0.5;
	i=2;j=0; vcl[j*3 + i] = 0.0;

	i=0;j=1; vcl[j*3 + i] = -1.0;
	i=1;j=1; vcl[j*3 + i] = 1.0;
	i=2;j=1; vcl[j*3 + i] = 0.0;

	i=0;j=2; vcl[j*3 + i] = 0.0;
	i=1;j=2; vcl[j*3 + i] = 0.0;
	i=2;j=2; vcl[j*3 + i] = 1.0;

	i=0;j=3; vcl[j*3 + i] = 0.0;
	i=1;j=3; vcl[j*3 + i] = 0.0;
	i=2;j=3; vcl[j*3 + i] = -1.0;

	i=0;j=4; vcl[j*3 + i] = 0.0;
	i=1;j=4; vcl[j*3 + i] = 1.0;
	i=2;j=4; vcl[j*3 + i] = 0.0;

	// InOut
	int vm[npt];

	// Remember Fortran is option base 1.
	vm[0] = 1;
	vm[1] = 2;
	vm[2] = 3;
	vm[3] = 4;
	vm[4] = 5;

	// Output
	int nbf;
	int nfc;
	int nface;
	int ntetra;
	int bf[3*maxbf];
	int fc[7*maxfc];
	int ht[sizht];
	int ierror;

	cout << "dtris3_c start" << endl << flush;
	dtris3_c(
			&npt,
			&sizht,
			&maxbf,
			&maxfc,
			vcl, //
			vm,  //
			&nbf,
			&nfc,
			&nface,
			&ntetra,
			bf, //
			fc, //
			ht,  //
			&ierror
			);
	cout << "ierror " << ierror << endl << flush;
	cout << "dtris3_c end" << endl << flush;

	int nt; // number of tetrahedra
	int maxNt = 1000;
	int tetra[4*maxNt];
	tetlst_c(
			&nfc,
			vm,
			fc,
			&nt,
			tetra
			);

	viz->triaxis();

	cout << "ntetra,nt: " << ntetra << ", " << nt << endl << flush;

	for(int iTet=0; iTet < nt; iTet++) {
		cout << "iTet: " << iTet << endl << flush;
		// vcl[j*3+i]
		// tetra[iTet*4 + iEdge]
		int iVerts[4];
		for(int iVert=0; iVert < 4; iVert++) {
			int iV = tetra[iTet*4+iVert];
			iVerts[iVert] = iV - 1; // Fortran is option base 1.
		}
		SpatialVector verts[4];
		for(int iVert=0; iVert<4; iVert++) {
			verts[iVert]
				  = SpatialVector(vcl[iVerts[iVert]*3],vcl[iVerts[iVert]*3+1],vcl[iVerts[iVert]*3+2]);
		}
		for(int i=0; i < 4; i++) {
			for(int j=i+1; j < 4; j++ ) {
				cout << "adding i,j: " << i << ", " << j << " : " << iVerts[i] << ", " << iVerts[j] << endl << flush;
				if(iTet == 0){
					viz->addEdge(verts[i],verts[j],0.0,0.5,0.5);
				} else if( iTet == 1) {
					viz->addEdge(verts[i],verts[j],0.5,0.5,0.0);
				} else {
					viz->addEdge(verts[i],verts[j],0.0,1.0,0.0);
				}

			}
		}
		cout << endl << flush;
	}
}



void testDelaunay1(VizHTM *viz) {


	// plotBlockingSphere(viz,0.2,0.2,0.2,0.999);
	testTenDegreeGridRGB(viz,0.6,0.6,0.6);


	// Input
	int npt     = 10;
	int hashSize = 3*npt/2;
	cout << "Calling prime_c" << flush;
	int sizht   = prime_c(&hashSize);
	cout << ". " << sizht << endl << flush;
	int maxbf   = 10000;
	int maxfc   = 10000;
	int maxPt   = 10000; // npt <= maxPt
	float64 vcl[3*maxPt];

	int i,j;
	for(j=0;j<npt;j++) {
		SpatialVector a = randomVector(); a.normalize();
		i=0; vcl[j*3 + i] = a.x();
		i=1; vcl[j*3 + i] = a.y();
		i=2; vcl[j*3 + i] = a.z();
	}

	// InOut
	int vm[npt];

	// Remember Fortran is option base 1.
	for(int ipt=0; ipt<npt;ipt++) {
		vm[ipt] = ipt + 1;
	}

	// Output
	int nbf;
	int nfc;
	int nface;
	int ntetra;
	int bf[3*maxbf];
	int fc[7*maxfc];
	int ht[sizht];
	int ierror;

	cout << "dtris3_c start" << endl << flush;
	dtris3_c(
			&npt,
			&sizht,
			&maxbf,
			&maxfc,
			vcl, //
			vm,  //
			&nbf,
			&nfc,
			&nface,
			&ntetra,
			bf, //
			fc, //
			ht,  //
			&ierror
			);
	cout << "ierror " << ierror << endl << flush;
	cout << "nbf    " << nbf    << endl << flush;
	cout << "nfc    " << nfc    << endl << flush;
	cout << "dtris3_c end" << endl << flush;

	if(true){
		cout << " boundary start " << endl << flush;
		 plotBlockingSphere(viz,0.2,0.2,0.2,0.999);

		for(int iFC=0; iFC<nfc; iFC++) {
			cout << " iFC: " << iFC << flush;
			int A = fc[7*iFC+0];
			cout << " A: " << A << flush;
			int B = fc[7*iFC+1];
			cout << " B: " << B << flush;
			int g7 = fc[7*iFC+6];
			cout << " g7: " << g7 << flush;
			if(A>0) {
				int E = fc[7*iFC+4]; // aka boundaryp
				cout << " E: " << E;
				if(E<0){
					SpatialVector verts[4]; // 1..3 the face, 4 the interior vert
					for(int iV=0; iV<4; iV++) {
						int icl = vm[fc[iFC*7+iV]-1] - 1; // Fortran Fortran!!!
						cout << " icl: " << icl << flush;
						float64 x = vcl[icl*3+0];
						float64 y = vcl[icl*3+1];
						float64 z = vcl[icl*3+2];
						verts[iV] = SpatialVector(x,y,z);
					}
					cout << "." << flush;
					// Cross a couple
					SpatialVector v01 = verts[1] - verts[0];
					SpatialVector v02 = verts[2] - verts[0];
					SpatialVector v01c02 = v01 ^ v02;
					SpatialVector v012   = (verts[0] + verts[1] + verts[2])*(1.0/3.0);
					SpatialVector v012to3 = verts[3] - v012;
					if(v012to3 * v01c02 > 0.0) {
						v01c02 = v01c02 * -1.0;
					}
					SpatialVector vcross = v012 + v01c02;
					if(v012 * v01c02 > 0.0) {
						// These triangles are on the surface of the sphere. Final step: check that they're interior.
						// This is a triangulation of the convex hull of an arbitrary set of points on the sphere.
						SpatialVector tmp = v012; tmp.normalize();
						viz->addEdge(tmp,vcross,1.0,1.0,1.0); // draw the normal starting from the centroid
						viz->addEdge(SpatialVector(0,0,0),tmp,0.5,0.5,0.8); // draw to face centroid
						for(int i=0; i<3; i++) {
							for(int j=i+1; j<3; j++) {
								viz->addArc(verts[i],verts[j],0.9,0.0,0.0); // arc are better reps of the region
							}
						}
					}
				}
			}
			cout << endl << flush;
		}
		cout << " boundary done " << endl << flush;
	}

	if(false) {
		int nt; // number of tetrahedra
		int maxNt = 1000;
		int tetra[4*maxNt];
		tetlst_c(
				&nfc,
				vm,
				fc,
				&nt,
				tetra
		);

		viz->triaxis();

		cout << "ntetra,nt: " << ntetra << ", " << nt << endl << flush;

		for(int iTet=0; iTet < nt; iTet++) {
			cout << "iTet: " << iTet << endl << flush;
			// vcl[j*3+i]
			// tetra[iTet*4 + iEdge]
			int iVerts[4];
			for(int iVert=0; iVert < 4; iVert++) {
				int iV = tetra[iTet*4+iVert];
				iVerts[iVert] = iV - 1; // Fortran is option base 1.
			}
			SpatialVector verts[4];
			for(int iVert=0; iVert<4; iVert++) {
				verts[iVert]
					  = SpatialVector(vcl[iVerts[iVert]*3],vcl[iVerts[iVert]*3+1],vcl[iVerts[iVert]*3+2]);
			}
			for(int i=0; i < 4; i++) {
				for(int j=i+1; j < 4; j++ ) {
					cout << "adding i,j: " << i << ", " << j << " : " << iVerts[i] << ", " << iVerts[j] << endl << flush;
					if(iTet == 0){
						viz->addEdge(verts[i],verts[j],0.0,0.5,0.5);
					} else if( iTet == 1) {
						viz->addEdge(verts[i],verts[j],0.5,0.5,0.0);
					} else {
						viz->addEdge(verts[i],verts[j],0.0,1.0,0.0);
					}

				}
			}
			cout << endl << flush;
		}
	}
}


void testDelaunay(VizHTM *viz) {


	plotBlockingSphere(viz,0.2,0.2,0.2,0.999);
	testTenDegreeGridRGB(viz,0.6,0.6,0.6);


	// Input
	int npt     = 12;
	int hashSize = 3*npt/2;
	cout << "Calling prime_c" << flush;
	int sizht   = prime_c(&hashSize);
	cout << ". " << sizht << endl << flush;
	int maxbf   = 10000;
	int maxfc   = 10000;
	int maxPt   = 10000; // npt <= maxPt
	float64 vcl[3*maxPt];

	{int i,j;
	j = 0;
	SpatialVector a;
	a.setLatLonDegrees(30.0,30.0);
	i=0; vcl[j*3 + i] = a.x();
	i=1; vcl[j*3 + i] = a.y();
	i=2; vcl[j*3 + i] = a.z();}

	{int i,j;
	j = 1;
	SpatialVector a;
	a.setLatLonDegrees(30.0,60.0);
	i=0; vcl[j*3 + i] = a.x();
	i=1; vcl[j*3 + i] = a.y();
	i=2; vcl[j*3 + i] = a.z();}

	{int i,j;
	j = 2;
	SpatialVector a;
	a.setLatLonDegrees(60.0,45.0);
	i=0; vcl[j*3 + i] = a.x();
	i=1; vcl[j*3 + i] = a.y();
	i=2; vcl[j*3 + i] = a.z();}

	{int i,j;
	j = 3;
	SpatialVector a;
	a.setLatLonDegrees(45.0,30.0);
	i=0; vcl[j*3 + i] = a.x();
	i=1; vcl[j*3 + i] = a.y();
	i=2; vcl[j*3 + i] = a.z();}

	{int i,j;
	j = 4;
	SpatialVector a;
	a.setLatLonDegrees(45.0,60.0);
	i=0; vcl[j*3 + i] = a.x();
	i=1; vcl[j*3 + i] = a.y();
	i=2; vcl[j*3 + i] = a.z();}

	{int i,j;
	j = 5;
	SpatialVector a;
	a.setLatLonDegrees(15.0,45);
	i=0; vcl[j*3 + i] = a.x();
	i=1; vcl[j*3 + i] = a.y();
	i=2; vcl[j*3 + i] = a.z();}

	{int i,j;
	j = 6;
	SpatialVector a;
	a.setLatLonDegrees(45.0,37.5);
	i=0; vcl[j*3 + i] = a.x();
	i=1; vcl[j*3 + i] = a.y();
	i=2; vcl[j*3 + i] = a.z();}

	{int i,j;
	j = 7;
	SpatialVector a;
	a.setLatLonDegrees(30.0,37.5);
	i=0; vcl[j*3 + i] = a.x();
	i=1; vcl[j*3 + i] = a.y();
	i=2; vcl[j*3 + i] = a.z();}

	{int i,j;
	j = 8;
	SpatialVector a;
	a.setLatLonDegrees(45.0,52.5);
	i=0; vcl[j*3 + i] = a.x();
	i=1; vcl[j*3 + i] = a.y();
	i=2; vcl[j*3 + i] = a.z();}

	{int i,j;
	j = 9;
	SpatialVector a;
	a.setLatLonDegrees(30.0,52.5);
	i=0; vcl[j*3 + i] = a.x();
	i=1; vcl[j*3 + i] = a.y();
	i=2; vcl[j*3 + i] = a.z();}

	{int i,j;
	j = 10;
	SpatialVector a;
	a.setLatLonDegrees(37.5,33.75);
	i=0; vcl[j*3 + i] = a.x();
	i=1; vcl[j*3 + i] = a.y();
	i=2; vcl[j*3 + i] = a.z();}

	{int i,j;
	j = 11;
	SpatialVector a;
	a.setLatLonDegrees(37.5,56.25);
	i=0; vcl[j*3 + i] = a.x();
	i=1; vcl[j*3 + i] = a.y();
	i=2; vcl[j*3 + i] = a.z();}

	// InOut
	int vm[npt];

	// Remember Fortran is option base 1.
	for(int ipt=0; ipt<npt;ipt++) {
		vm[ipt] = ipt + 1;
	}

	// Output
	int nbf;
	int nfc;
	int nface;
	int ntetra;
	int bf[3*maxbf];
	int fc[7*maxfc];
	int ht[sizht];
	int ierror;

	cout << "dtris3_c start" << endl << flush;
	dtris3_c(
			&npt,
			&sizht,
			&maxbf,
			&maxfc,
			vcl, //
			vm,  //
			&nbf,
			&nfc,
			&nface,
			&ntetra,
			bf, //
			fc, //
			ht,  //
			&ierror
			);
	cout << "ierror " << ierror << endl << flush;
	cout << "nbf    " << nbf    << endl << flush;
	cout << "nfc    " << nfc    << endl << flush;
	cout << "dtris3_c end" << endl << flush;

	if(true){
		cout << " boundary start " << endl << flush;
//		 plotBlockingSphere(viz,0.2,0.2,0.2,0.999);

		for(int iFC=0; iFC<nfc; iFC++) {
			cout << " iFC: " << iFC << flush;
			int A = fc[7*iFC+0];
			cout << " A: " << A << flush;
			int B = fc[7*iFC+1];
			cout << " B: " << B << flush;
			int g7 = fc[7*iFC+6];
			cout << " g7: " << g7 << flush;
			if(A>0) {
				int E = fc[7*iFC+4]; // aka boundaryp
				cout << " E: " << E;
				if(E<0){
					SpatialVector verts[4]; // 1..3 the face, 4 the interior vert
					for(int iV=0; iV<4; iV++) {
						int icl = vm[fc[iFC*7+iV]-1] - 1; // Fortran Fortran!!!
						cout << " icl: " << icl << flush;
						float64 x = vcl[icl*3+0];
						float64 y = vcl[icl*3+1];
						float64 z = vcl[icl*3+2];
						verts[iV] = SpatialVector(x,y,z);
					}
					cout << "." << flush;
					// Cross a couple
					SpatialVector v01 = verts[1] - verts[0];
					SpatialVector v02 = verts[2] - verts[0];
					SpatialVector v01c02 = v01 ^ v02;
					SpatialVector v012   = (verts[0] + verts[1] + verts[2])*(1.0/3.0);
					SpatialVector v012to3 = verts[3] - v012;
					if(v012to3 * v01c02 > 0.0) {
						v01c02 = v01c02 * -1.0;
					}
					SpatialVector vcross = v012 + v01c02;
					if(v012 * v01c02 > 0.0) {
						// These triangles are on the surface of the sphere. Final step: check that they're interior.
						// This is a triangulation of the convex hull of an arbitrary set of points on the sphere.
						SpatialVector tmp = v012; tmp.normalize();
						viz->addEdge(tmp,vcross,1.0,1.0,1.0); // draw the normal starting from the centroid
						viz->addEdge(SpatialVector(0,0,0),tmp,0.5,0.5,0.8); // draw to face centroid
						for(int i=0; i<3; i++) {
							for(int j=i+1; j<3; j++) {
								viz->addArc(verts[i],verts[j],0.9,0.0,0.0); // arc are better reps of the region
							}
						}
					}
				}
			}
			cout << endl << flush;
		}
		cout << " boundary done " << endl << flush;
	}

	if(false) {
		int nt; // number of tetrahedra
		int maxNt = 1000;
		int tetra[4*maxNt];
		tetlst_c(
				&nfc,
				vm,
				fc,
				&nt,
				tetra
		);

		viz->triaxis();

		cout << "ntetra,nt: " << ntetra << ", " << nt << endl << flush;

		for(int iTet=0; iTet < nt; iTet++) {
			cout << "iTet: " << iTet << endl << flush;
			// vcl[j*3+i]
			// tetra[iTet*4 + iEdge]
			int iVerts[4];
			for(int iVert=0; iVert < 4; iVert++) {
				int iV = tetra[iTet*4+iVert];
				iVerts[iVert] = iV - 1; // Fortran is option base 1.
			}
			SpatialVector verts[4];
			for(int iVert=0; iVert<4; iVert++) {
				verts[iVert]
					  = SpatialVector(vcl[iVerts[iVert]*3],vcl[iVerts[iVert]*3+1],vcl[iVerts[iVert]*3+2]);
			}
			for(int i=0; i < 4; i++) {
				for(int j=i+1; j < 4; j++ ) {
					cout << "adding i,j: " << i << ", " << j << " : " << iVerts[i] << ", " << iVerts[j] << endl << flush;
					if(iTet == 0){
						viz->addEdge(verts[i],verts[j],0.0,0.5,0.5);
					} else if( iTet == 1) {
						viz->addEdge(verts[i],verts[j],0.5,0.5,0.0);
					} else {
						viz->addEdge(verts[i],verts[j],0.0,1.0,0.0);
					}

				}
			}
			cout << endl << flush;
		}
	}
}

void testNearestNeighbors(VizHTM *viz) {
	int level = 5;
	SpatialIndex index(level,5);

	{
		string name="N012301";

		testHTMRange(viz,level,name.c_str(),name.c_str());
		uint64 htmId = index.idByName(name.c_str());

		uint64 neighbors[3];
		index.NeighborsAcrossEdgesFromHtmId(neighbors,htmId);
		for(int j=0; j<3; j++) {
			char tmpName[64];
			index.nameById(neighbors[j],tmpName);
			testHTMRange(viz,level,tmpName,tmpName);
		}

		uint64 neighborsV[9];
		index.NeighborsAcrossVerticesFromHtmId(neighborsV,htmId);
		for(int j=0; j<9; j++) {
			char tmpName[64];
			index.nameById(neighborsV[j],tmpName);
			testHTMRange(viz,level,tmpName,tmpName);
		}
	}

	{
		string name="N011301";

		testHTMRange(viz,level,name.c_str(),name.c_str());
		uint64 htmId = index.idByName(name.c_str());

		uint64 neighbors[3];
		index.NeighborsAcrossEdgesFromHtmId(neighbors,htmId);
		for(int j=0; j<3; j++) {
			char tmpName[64];
			index.nameById(neighbors[j],tmpName);
			testHTMRange(viz,level,tmpName,tmpName);
		}

	}

	{
		string name="N010301";

		testHTMRange(viz,level,name.c_str(),name.c_str());
		uint64 htmId = index.idByName(name.c_str());

		uint64 neighborsV[9];
		index.NeighborsAcrossVerticesFromHtmId(neighborsV,htmId);
		for(int j=0; j<9; j++) {
			char tmpName[64];
			index.nameById(neighborsV[j],tmpName);
			testHTMRange(viz,level,tmpName,tmpName);
		}
	}

}

// TODO Calculate the triangle-covering of a circle.

struct colorRGBAS {
	colorRGBAS(
			float ri = 0,
			float gi = 0,
			float bi = 0,
			float ai = -1, // Not generally implemented.
			float scalei = 1.0) {
		r=ri; g=gi; b=bi; a=ai; scale=scalei;
	}; // Scale is used to scale the unit vector for layering graphics.
	float r,g,b,a,scale;
};

colorRGBAS swathStratiformColor            (0.8,0.8,0.0,0.5,1.0007);
colorRGBAS swathConvectiveColor            (0.0,0.0,0.8,0.5,1.0007);
colorRGBAS swathOtherColor                 ();
colorRGBAS swathConvectiveHTMGeometryColor (1.0,0.8,0.8,-1.0,1.0012);

colorRGBAS swathGeometryColor        (0.4,0.6,0.0,-1.0,1.0008);
colorRGBAS swathHTMRGeometryColor    (0.8,0.9,0.0,-1.0,1.001);

colorRGBAS nmqColor                  (0.7,0.0,0.0,-1.0,1.0005);
colorRGBAS nmqGridColor              (0.2,0.8,0.2,-1.0,1.0006); // Larger step scale
colorRGBAS nmqLatLonBoxStepScaleColor(0.3,1.0,0.3,-1.0,1.0007); // Smaller step scale
colorRGBAS nmqHTMRGeometryColor      (0.6,0.6,1.0,-1.0,1.0009);

colorRGBAS nmq_trmm_intersectionColor        (1.0,0.1,0.1,-1.0,1.00119);
colorRGBAS nmq_trmm_PartialIntersectionColor (1.0,0.5,0.8,-1.0,1.00111);
colorRGBAS trmm_nmq_PartialIntersectionColor (1.0,0.5,0.7,-1.0,1.00112);

bool swathGeometry_viz               = false;
bool swathStratiform_viz             = false;
bool swathConvective_viz             = false;
bool swathOtherRainType_viz          = false;
bool swathConvectiveHTM_viz          = false;
bool swathGeometryHTM_viz            = false;

int swathIndexStart = 190000;
int swathIndexEnd   = 260000;

bool nmqHTM                          = false; ///< Calculate nmqRange at high resolution;

bool nmqLatLonBox_viz                = false; // At the larger geomStepScale
bool nmqLatLonBoxStepScale_viz       = false; // At the smaller stepScale
bool nmqPrecipRateHSR_viz            = false;
bool nmqPrecipRateHSRGeometry_viz    = false;
bool nmqNoSignalGeometry_viz         = false;

bool nmqHTMRGeometry_viz             = false;

bool nmq_trmm_HTMRPartial_viz        = false;
bool nmq_trmm_HTMRDoesNotContain_viz = false;
bool nmq_trmm_HTMRContains_viz       = false;

bool trmm_nmq_HTMRPartial_viz        = false;
bool trmm_nmq_HTMRDoesNotContain_viz = false;
bool trmm_nmq_HTMRContains_viz       = false;

int nmqLatStart = 1900; int nmqLatEnd = 2450;
int nmqLonStart = 4450; int nmqLonEnd = 5550;

bool offscreen_viz                   = false;
bool examiner_viz                    = false;

bool addFiduciaryTriangles           = false;

SpatialVector *nmq_trmm_vizCenter = VectorFromLatLonDegrees(33.63,-76.575);
float offscreenPerscamPositionScale = 1.1;
float examinerCameraPositionScale   = 1.17;

//uint64          swathLevel = 6; // Good for testing
// uint64           swathLevel = 7;
//uint64           swathLevel = 9;
uint64           swathLevel = 7;
//uint64           swathLevel = 11;
uint64          buildLevel = 5;
// uint64          buildLevel = 11; // default is 5;
SpatialIndex  swathIndex(swathLevel,buildLevel);

string baseName = "TrmmNmq";

void initializeFocusOnNMQSquare() {
	swathIndexStart = swathIndexEnd - 35000;
	swathIndexEnd   = swathIndexEnd - 4000;
	offscreenPerscamPositionScale = 1.02;
	examinerCameraPositionScale   = 1.01;
}


void initializeFocusOnNMQCorner() {

	bool cornerFocus = true; // 34.245,-74.51

	examinerCameraPositionScale   = 1.01;
	offscreenPerscamPositionScale = 1.25;

	bool narrow = true;
	string slide = "8a1";

	if( slide == "1") {
		swathLevel = 9;
		offscreenPerscamPositionScale = 1.1;
		nmqHTM = true; // was false...
	} else if( slide == "2") {
		swathLevel = 10;
		offscreenPerscamPositionScale = 1.02;
	} else if ( slide == "3" ) {
		swathLevel = 11;
		offscreenPerscamPositionScale = 1.02;
		nmqHTM = false;  // Use full NMQ resolution?
	} else if ( slide == "4" ) {
		swathLevel = 11;
		offscreenPerscamPositionScale = 1.02;
		nmqHTM = true;  // Use full NMQ resolution?
	} else if ( slide == "4a" ) {
		swathLevel = 9;
		offscreenPerscamPositionScale = 1.02;
		nmqHTM = true;  // Use full NMQ resolution?
	} else if ( slide == "4b" ) {
	// Slide 4b Same as 4 with nmq grid
		swathLevel = 11;
		offscreenPerscamPositionScale = 1.02;
		nmqHTM = true;  // Use full NMQ resolution?
//		nmqLatLonBox_viz = true; // blocks of ten pixels
		nmqLatLonBoxStepScale_viz = true;
	} else if ( slide == "4c" ) {
		narrow = false;
		// Make even narrower
		swathIndexStart = swathIndexEnd - 12000;
		swathIndexEnd   = swathIndexEnd - 8000;
		nmqLatStart = 2000;
		nmqLatEnd   = 2150;
		nmqLonStart = 5250;
		nmqLonEnd   = 5550;
//		// While aligning grid.
//		swathLevel = 7;
//		offscreenPerscamPositionScale = 1.25;
//		nmqHTM = false;
	} else if ( slide == "5" ) {
		// Slide 5
		swathLevel = 12;
		offscreenPerscamPositionScale = 1.02;
		nmqHTM = true;  // Use full NMQ resolution?
	} else if ( slide == "6" ) {
	// Slide 6
		swathLevel = 8;
		offscreenPerscamPositionScale = 1.03;
		nmqHTM = false;  // Use full NMQ resolution?
	} else if ( slide == "7" ) {
		swathLevel = 7;
		offscreenPerscamPositionScale = 1.25;
		nmqHTM = false;
		narrow = false;
		nmq_trmm_HTMRContains_viz        = true;
		trmm_nmq_HTMRContains_viz        = true;
		nmq_trmm_HTMRDoesNotContain_viz  = true;
		trmm_nmq_HTMRDoesNotContain_viz  = true;
		nmq_trmm_HTMRPartial_viz         = true;
		trmm_nmq_HTMRPartial_viz         = true;
	} else if ( slide == "7a" ) {
		swathLevel = 7;
		offscreenPerscamPositionScale = 1.25;
		nmqHTM = false;
		narrow = false;
		nmq_trmm_HTMRContains_viz        = false;
		trmm_nmq_HTMRContains_viz        = false;
		nmq_trmm_HTMRDoesNotContain_viz  = false;
		trmm_nmq_HTMRDoesNotContain_viz  = false;
		nmq_trmm_HTMRPartial_viz         = false;
		trmm_nmq_HTMRPartial_viz         = false;
	} else if ( slide == "8" ) {
		swathLevel = 5;
		offscreenPerscamPositionScale = 1.25;
		nmqHTM = false;
		narrow = false;
		nmq_trmm_HTMRContains_viz        = true;
		trmm_nmq_HTMRContains_viz        = true;
		nmq_trmm_HTMRDoesNotContain_viz  = true;
		trmm_nmq_HTMRDoesNotContain_viz  = true;
		nmq_trmm_HTMRPartial_viz         = true;
		trmm_nmq_HTMRPartial_viz         = true;
	} else if ( slide == "8a" ) {
		swathLevel = 5;
		offscreenPerscamPositionScale = 1.25;
		examinerCameraPositionScale = offscreenPerscamPositionScale;
		nmqHTM = false;
		narrow = false;
		cornerFocus = false;
		nmq_trmm_HTMRContains_viz        = true;
		trmm_nmq_HTMRContains_viz        = true;
		nmq_trmm_HTMRDoesNotContain_viz  = true;
		trmm_nmq_HTMRDoesNotContain_viz  = true;
		nmq_trmm_HTMRPartial_viz         = true;
		trmm_nmq_HTMRPartial_viz         = true;
	} else if ( slide == "8a1" ) {
		swathLevel = 6;
		offscreenPerscamPositionScale = 1.25;
		examinerCameraPositionScale = offscreenPerscamPositionScale;
		nmqHTM = false;
		narrow = false;
		cornerFocus = false;
		nmq_trmm_HTMRContains_viz        = true;
		trmm_nmq_HTMRContains_viz        = true;
		nmq_trmm_HTMRDoesNotContain_viz  = true;
		trmm_nmq_HTMRDoesNotContain_viz  = true;
		nmq_trmm_HTMRPartial_viz         = true;
		trmm_nmq_HTMRPartial_viz         = true;
	} else if ( slide == "8b" ) {
		swathLevel = 7;
		offscreenPerscamPositionScale = 1.25;
		examinerCameraPositionScale = offscreenPerscamPositionScale;
		nmqHTM = false;
		narrow = false;
		cornerFocus = false;
		nmq_trmm_HTMRContains_viz        = true;
		trmm_nmq_HTMRContains_viz        = true;
		nmq_trmm_HTMRDoesNotContain_viz  = true;
		trmm_nmq_HTMRDoesNotContain_viz  = true;
		nmq_trmm_HTMRPartial_viz         = true;
		trmm_nmq_HTMRPartial_viz         = true;
	} else if ( slide == "8b1" ) {
		swathLevel = 8;
		offscreenPerscamPositionScale = 1.25;
		examinerCameraPositionScale = offscreenPerscamPositionScale;
		nmqHTM = false;
		narrow = false;
		cornerFocus = false;
		nmq_trmm_HTMRContains_viz        = true;
		trmm_nmq_HTMRContains_viz        = true;
		nmq_trmm_HTMRDoesNotContain_viz  = true;
		trmm_nmq_HTMRDoesNotContain_viz  = true;
		nmq_trmm_HTMRPartial_viz         = true;
		trmm_nmq_HTMRPartial_viz         = true;
	} else if( slide == "8c1" ) {
		swathLevel = 8;
		offscreenPerscamPositionScale = 2.5;
		examinerCameraPositionScale = offscreenPerscamPositionScale;
		nmqHTM = false;
		narrow = false;
		cornerFocus = false;
		nmq_trmm_HTMRContains_viz        = true;
		trmm_nmq_HTMRContains_viz        = true;
		nmq_trmm_HTMRDoesNotContain_viz  = true;
		trmm_nmq_HTMRDoesNotContain_viz  = true;
		nmq_trmm_HTMRPartial_viz         = true;
		trmm_nmq_HTMRPartial_viz         = true;

		addFiduciaryTriangles            = true;

	} else if( slide == "8c" ) {
		swathLevel = 9;
		offscreenPerscamPositionScale = 1.25;
		examinerCameraPositionScale = offscreenPerscamPositionScale;
		nmqHTM = false;
		narrow = false;
		cornerFocus = false;
		nmq_trmm_HTMRContains_viz        = true;
		trmm_nmq_HTMRContains_viz        = true;
		nmq_trmm_HTMRDoesNotContain_viz  = true;
		trmm_nmq_HTMRDoesNotContain_viz  = true;
		nmq_trmm_HTMRPartial_viz         = true;
		trmm_nmq_HTMRPartial_viz         = true;
	} else if( slide == "8d" ) {
		swathLevel = 9;
		offscreenPerscamPositionScale = 1.25;
		examinerCameraPositionScale = offscreenPerscamPositionScale;
		nmqHTM = false;
		narrow = false;
		cornerFocus = false;
		nmq_trmm_HTMRContains_viz        = true;
		trmm_nmq_HTMRContains_viz        = true;
		nmq_trmm_HTMRDoesNotContain_viz  = true;
		trmm_nmq_HTMRDoesNotContain_viz  = true;
		nmq_trmm_HTMRPartial_viz         = true;
		trmm_nmq_HTMRPartial_viz         = true;

		nmqLatLonBox_viz                 = true;

	} else {
		cout << "slide " << slide << " not found, exiting. " << endl << flush;
		exit( 1 );
	}

	if(narrow){
		// Narrow to corner
		swathIndexStart = swathIndexEnd - 15000;
		swathIndexEnd   = swathIndexEnd - 4500;
		//	nmqLatStart = 2000;
		//	nmqLatEnd   = 2350;
		nmqLonStart = 4800;
		nmqLonEnd   = 5550;
	}

	cout << "writing slide " << slide << endl << flush;

	if(cornerFocus) {
		nmq_trmm_vizCenter = VectorFromLatLonDegrees(34.245,-74.51);
	}

	baseName +=
			"-Level="+std::to_string(swathLevel)
	+"-PositionScale="+std::to_string(offscreenPerscamPositionScale)
	+"-nmqHtm="+std::to_string(nmqHTM)
	+"-Slide="+slide;

}

void initializeZoomCornerVisualization() {

	offscreen_viz = false; examiner_viz = !offscreen_viz;

	//	nmqHTMRGeometry_viz = true;

	nmqPrecipRateHSR_viz = true;

	swathStratiform_viz  = true;
	swathConvective_viz  = true;
	swathOtherRainType_viz = true;

	nmq_trmm_HTMRContains_viz        = true;
	trmm_nmq_HTMRContains_viz        = true;

	nmq_trmm_HTMRDoesNotContain_viz  = true;
	trmm_nmq_HTMRDoesNotContain_viz  = true;

	nmq_trmm_HTMRPartial_viz = true;
	trmm_nmq_HTMRPartial_viz = true;

	initializeFocusOnNMQCorner();
	//	initializeFocusOnNMQSquare();

	swathIndex = SpatialIndex(swathLevel,buildLevel);

}

void initialize_nmq_trmm_viz() {

	initializeZoomCornerVisualization();

}


HtmRange      swathRange;
bool          swath_varlenHTMIDs = false;

HtmRange      convectiveRainRange;

void testSwath(VizHTM *viz) {
	ios::fmtflags coutInitialState(std::cout.flags());

	std::ifstream lat_ifs( "/Users/mrilee/data/DERECHOS/TRMM/2A23.20091231.69087.7.Latitude",
			std::ios::binary );

	std::ifstream lon_ifs( "/Users/mrilee/data/DERECHOS/TRMM/2A23.20091231.69087.7.Longitude",
			std::ios::binary );

	std::ifstream rainType_ifs( "/Users/mrilee/data/DERECHOS/TRMM/2A23.20091231.69087.7.rainType",
			std::ios::binary );

	float lat = -999; float lon = -999; short int rainType = -999;

	SpatialVector v;
	float64 const km = 1.0;
	float64 const RE = 6371.0; // mean radius km
	float64 gFOVRadians = 5.0 * km / RE; // Could use asin, but the accuracy does not deserve it.
//	float64 PIo2 = 0.5*atan2(0,-1);

	int iStart = swathIndexStart;
	int iEnd   = swathIndexEnd;

	// Swath inside the nmq rectangle.
//	iStart = iEnd - 21000;
//	iEnd = iEnd - 12000;

	cout << "TRMM swath: iStart, iEnd: " << iStart << ", " << iEnd << endl << flush;

	// DEBUGGING
//	iEnd = iStart + 100;

	swathRange.purge();
	convectiveRainRange.purge();

	bool tapFileFlag = true;
	ofstream tapFile;
	if(tapFileFlag) {
		tapFile.open("/Users/mrilee/Downloads-1/swathTapFile.csv",ios::out|ios::trunc);
	}

	for(int i=0; i< 9247*49; i++) {
		lat_ifs.read(reinterpret_cast<char*>(&lat),sizeof(float));
		lon_ifs.read(reinterpret_cast<char*>(&lon),sizeof(float));
		rainType_ifs.read(reinterpret_cast<char*>(&rainType),sizeof(short int));
		if(iStart <= i && i < iEnd) {

			if(tapFileFlag) {
				tapFile << "( "<< lon << "," << lat << "," << rainType << " ),"<< endl << flush;
			}

			v.setLatLonDegrees(lat,lon);

			if(true){
				SpatialConstraint c(v,cos(0.5*gFOVRadians));
				RangeConvex rc; rc.add(c);
				SpatialDomain d; d.add(rc);
				HtmRange r;
				r.purge();
				bool overlap = d.intersect(&swathIndex,&r,swath_varlenHTMIDs);
				r.reset();
				Key lo=-1, hi=-1;
				int indexp = r.getNext(lo,hi);
				if(indexp) {
					do {
						swathRange.addRange(lo,hi);
					} while(r.getNext(lo,hi));
				}
			}

//			viz->addConstraint(c,1.0,1.0,1.0);  // Verified constraint has correct geometry.
			if(true) {
				//			viz->addCircleFacet(v,0.5*gFOVRadians,1.0,1.0,1.0);
				if(100 <= rainType && rainType < 200 ) { // stratiform tags
					//					viz->addCircleFacet(v,0.5*gFOVRadians,0.9,0.9,0.0,0.5);
					if(swathStratiform_viz){
						viz->addCircleFacet(v,0.5*gFOVRadians,
								swathStratiformColor.r,
								swathStratiformColor.g,
								swathStratiformColor.b,
								swathStratiformColor.a,
								swathStratiformColor.scale
						);}
				} else if (200 <= rainType && rainType < 300 ) { // convective
					//					viz->addCircleFacet(v,0.5*gFOVRadians,0.0,0.0,0.9,0.5);
					if(swathConvective_viz) {
						viz->addCircleFacet(v,0.5*gFOVRadians,
								swathConvectiveColor.r,
								swathConvectiveColor.g,
								swathConvectiveColor.b,
								swathConvectiveColor.a,
								swathConvectiveColor.scale
						);}
					if(true){
						SpatialConstraint c(v,cos(0.5*gFOVRadians));
						RangeConvex rc; rc.add(c);
						SpatialDomain d; d.add(rc);
						HtmRange r;
						r.purge();
						bool overlap = d.intersect(&swathIndex,&r,swath_varlenHTMIDs);
						r.reset();
						Key lo=-1, hi=-1;
						int indexp = r.getNext(lo,hi);
						if(indexp) {
							do {
								convectiveRainRange.addRange(lo,hi);
							} while(r.getNext(lo,hi));
						}
					}

				} else if (300 <= rainType && rainType < 400 ) { // others
					if(swathOtherRainType_viz) {
						viz->addCircleFacet(v,0.5*gFOVRadians,0.0,0.5,0.5,0.5);
					}
				}
			}
			if(swathGeometry_viz) { // Geometry
				viz->addCircleEdges(v,0.5*gFOVRadians,
						swathGeometryColor.r,
						swathGeometryColor.g,
						swathGeometryColor.b
						);
			}
		}
	}

	if(tapFileFlag) {
		tapFile.close();
	}
	cout << "TRMM swath: HTM viz" << endl << flush;

//	swathDomain.add(swathConvex);
//	swathRange.purge();
//	bool overlap = swathDomain.intersect(&swathIndex,&swathRange,swath_varlenHTMIDs);
//	cout << "overlap: " << overlap << endl << flush;
//	cout << "nranges: " << swathRange.nranges() << endl << flush;

	if(swathGeometryHTM_viz) {
		swathRange.reset();
		//	// Iterate over components of swathRange.
		////	viz->addHTMRange(&swathIndex,&swathRange,0.0,1.0,0.0,-1,1.01);
		viz->addHTMRange(&swathIndex,&swathRange,
				swathHTMRGeometryColor.r,
				swathHTMRGeometryColor.g,
				swathHTMRGeometryColor.b,
				swathHTMRGeometryColor.a,
				swathHTMRGeometryColor.scale
		);
		//	cout << " swathGeometryColor.r "<< swathGeometryColor.r << endl << flush;
		//	cout << " swathGeometryColor.g "<< swathGeometryColor.g << endl << flush;
		//	cout << " swathGeometryColor.b "<< swathGeometryColor.b << endl << flush;
	}

	if(swathConvectiveHTM_viz) {
		convectiveRainRange.reset();
		viz->addHTMRange(&swathIndex,&convectiveRainRange,
				swathConvectiveHTMGeometryColor.r,
				swathConvectiveHTMGeometryColor.g,
				swathConvectiveHTMGeometryColor.b,
				swathConvectiveHTMGeometryColor.a,
				swathConvectiveHTMGeometryColor.scale
		);
	}

	lat_ifs.close();
	lon_ifs.close();
	rainType_ifs.close();

	std::cout.flags(coutInitialState);
	cout << "TRMM swath: done" << endl << flush;
}

HtmRange nmqRangeCoarse, nmqRange;

void testnmq(VizHTM *viz) {

	std::cout << " test nmq start " << endl << flush;

	const int nLat = 3501;
	const int nLon = 7001;

	std::cout << "allocating lat_ " << flush;
	float lat_[nLat];
	std::cout << " lon_ " << flush;
	float lon_[nLon];
//	std::cout << " preciprate_hsr_ " << flush;
//	float preciprate_hsr_[nLat*nLon];

	std::cout << "." << endl << flush;

	string baseName("/Users/mrilee/data/DERECHOS/nmq/netCDF/PRECIPRATE_HSR/20091231-132500.netcdf");

	std::cout << " latitude " << endl << flush;
	std::ifstream lat_ifs( baseName+".Latitude", std::ios::binary );
	for(int i=0; i< nLat; i++){
		lat_ifs.read(reinterpret_cast<char*>(&(lat_[i])),sizeof(float));
	}
	lat_ifs.close();

	std::cout << " longitude " << endl << flush;
	std::ifstream lon_ifs( baseName+".Longitude", std::ios::binary );
	for(int i=0; i< nLon; i++){
		lon_ifs.read(reinterpret_cast<char*>(&(lon_[i])),sizeof(float));
	}
	cout << "100" << endl << flush;
	lon_ifs.close();
	cout << "101" << endl << flush;
//
	std::cout << " radar " << endl << flush;

//	SpatialVector v;
	float64 const km = 1.0;
	float64 const RE = 6371.0; // mean radius km
	float64 gFOVRadians = 1.0 * km / RE; // Could use asin, but the accuracy does not deserve it.

	int latStart = nmqLatStart;
	int latEnd   = nmqLatEnd;

	int lonStart = nmqLonStart;
	int lonEnd   = nmqLonEnd;

//	// focus on east coast
//	latStart = 1900; latEnd = 2450;
//	lonStart = 4450; lonEnd = 5550;
//
//	// across us
////	latStart = 1750; latEnd = 3000;
////	lonStart = 2500; lonEnd = 5500;
//
////	latEnd=latStart+100;
////	lonEnd=lonStart+100;
//
//	// Smaller scale for testing failure at high resolution.
////	latEnd=latStart+40;
////	lonEnd=lonStart+40;
//
////		latEnd=latStart+41;
////		lonEnd=lonStart+41;
//
//	//	latEnd = 2005;
//	//	lonEnd = 5005;

	int   stepScale = 1;
	// stepScale = 10;
	float dLat = 0.005*stepScale; // TODO Note hardwired resolution in degrees.
	float dLon = 0.005*stepScale;

	int   geomStepScale = 10;
	// stepScale = 10;
	float geomDLat = 0.005*geomStepScale;
	float geomDLon = 0.005*geomStepScale;

	float mn = 1.0e9;
	float mx = -1.0e9;

	viz->lineWidth = 1.5;

	nmqRange.purge();
//	SpatialIndex index0(swathLevel,11);
//	swathIndex.setMaxlevel(swathLevel);

	uint64 nmqRangeCount = 0;
	cout << 200 << endl << flush;

	bool tapFileFlag = true;
	ofstream tapFile;
	if(tapFileFlag) {
		tapFile.open("/Users/mrilee/Downloads-1/nmqTapFile.csv",ios::out|ios::trunc);
	}

	std::ifstream preciprate_hsr_ifs( baseName+".PrecipitationRate_HSR", std::ios::binary );
	for(int jLat=0; jLat < nLat; jLat++) {
		for(int iLon=0; iLon < nLon; iLon++) {
//			v.setLatLonDegrees(lat_[jLat],lon_[iLon]);
			float preciprate_hsr;
			preciprate_hsr_ifs.read(reinterpret_cast<char*>(&preciprate_hsr),sizeof(float));

			// Overlay geometry
			if(true) {
				if( (jLat % geomStepScale == 0) && (iLon % geomStepScale == 0)) {
					if((latStart <= jLat) && (jLat <= latEnd) &&
							(lonStart <= iLon) && (iLon <= lonEnd)) {
						if(true) { // geometry
							if(nmqLatLonBox_viz){
								viz->addLatLonBoxEdgesDegrees(
										lat_[jLat]-geomDLat,lon_[iLon]-geomDLon,
										lat_[jLat]+geomDLat,lon_[iLon]+geomDLon,
										nmqGridColor.r,nmqGridColor.g,nmqGridColor.b
								);
							}
//							SpatialVector *v0 = VectorFromLatLonDegrees(lat_[jLat]-geomDLat,lon_[iLon]-geomDLon);
//							SpatialVector *v1 = VectorFromLatLonDegrees(lat_[jLat]-geomDLat,lon_[iLon]+geomDLon);
//							SpatialVector *v2 = VectorFromLatLonDegrees(lat_[jLat]+geomDLat,lon_[iLon]+geomDLon);
//							SpatialVector *v3 = VectorFromLatLonDegrees(lat_[jLat]+geomDLat,lon_[iLon]-geomDLon);
//							RangeConvex rc(v0,v1,v2,v3);
							SpatialVector v0, v1, v2, v3;
							 v0.setLatLonDegrees(lat_[jLat]-geomDLat,lon_[iLon]-geomDLon);
							 v1.setLatLonDegrees(lat_[jLat]-geomDLat,lon_[iLon]+geomDLon);
							 v2.setLatLonDegrees(lat_[jLat]+geomDLat,lon_[iLon]+geomDLon);
							 v3.setLatLonDegrees(lat_[jLat]+geomDLat,lon_[iLon]-geomDLon);
//							 viz->addRectangle(v0,v1,v2,v3,1.0,0.0,0.0);
							 RangeConvex rc = RangeConvex(&v0,&v1,&v2,&v3);
							SpatialDomain d; d.add(rc);
							HtmRange r;
							r.purge();
//							cout << jLat << " " << iLon << " : " << lat_[jLat] << ", " << lon_[iLon] << flush;
//							bool overlap = d.intersect(&index0,&r,swath_varlenHTMIDs);
							bool overlap = d.intersect(&swathIndex,&r,swath_varlenHTMIDs);
							Key lo=-999, hi=-999;
							r.reset();
							int indexp = r.getNext(lo,hi);
							if(indexp) {
								do {
									nmqRangeCoarse.addRange(lo,hi);
//									long losCount = nmqRange.my_los->getCount();
//									long hisCount = nmqRange.my_his->getCount();
//									cout << "( " << nmqRangeCount << ", "
//											<< losCount << ", "
//											<< hisCount << ", "
//											<< (losCount == hisCount)
//											<< " )" << endl << flush;
//									cout << " " << lo << ", " << hi << "; " << flush;
								} while(r.getNext(lo,hi));
							}
//							cout << endl << flush;
						}
					}
				}
			}


			if(true){
				// % mod the stepScale...
				if( (jLat % stepScale == 0) && (iLon % stepScale == 0)) {

					//			if(preciprate_hsr < mn) mn = preciprate_hsr;
					//			if(preciprate_hsr > mx) mx = preciprate_hsr;
					if((latStart <= jLat) && (jLat <= latEnd) &&
							(lonStart <= iLon) && (iLon <= lonEnd)) {
						//				cout << jLat << " " << iLon << " " << preciprate_hsr << endl << flush;

						tapFile
						<< "( "
						<< lon_[iLon] << ","
						<< lat_[jLat] << ","
						<< preciprate_hsr
						<< " ),"
						<< endl << flush;

						if(preciprate_hsr < mn) mn = preciprate_hsr;
						if(preciprate_hsr > mx) mx = preciprate_hsr;
						float scale = 150.0;
						scale = 55;
						if(nmqLatLonBoxStepScale_viz) { // geometry
							viz->addLatLonBoxEdgesDegrees(
									lat_[jLat]-dLat,lon_[iLon]-dLon,
									lat_[jLat]+dLat,lon_[iLon]+dLon,
									0.0,1.0,0.0
							);
						}
						if((0 < preciprate_hsr) && (preciprate_hsr <= scale)) {
							if(nmqPrecipRateHSR_viz) { // data
								float r = sqrt(sqrt(preciprate_hsr/scale));
								float g = sqrt(preciprate_hsr/scale);
								float b = 0.3;
								viz->addFace4FromLatLonDegrees(
										lat_[jLat]-dLat,lon_[iLon]-dLon, // fix with dLat, dLon
										lat_[jLat]+dLat,lon_[iLon]-dLon,
										lat_[jLat]+dLat,lon_[iLon]+dLon,
										lat_[jLat]-dLat,lon_[iLon]+dLon,
										r,g,b,
										r,g,b,
										r,g,b,
										r,g,b
								);
							}
							if(nmqPrecipRateHSRGeometry_viz) { // geometry
								viz->addLatLonBoxEdgesDegrees(
										lat_[jLat]-dLat,lon_[iLon]-dLon,
										lat_[jLat]+dLat,lon_[iLon]+dLon,
										0.0,1.0,0.0
								);
							}
							//					viz->addCircleFacet(
							//							v,
							//							0.5*gFOVRadians,
							//							preciprate_hsr/scale,
							//							preciprate_hsr/scale,
							//							0.1
							//							);
						} else { // no signal
							if(nmqNoSignalGeometry_viz) {
								float r = 0.15;
								float g = 0.15;
								float b = 0.20;
								viz->addFace4FromLatLonDegrees(
										lat_[jLat]-dLat,lon_[iLon]-dLon,
										lat_[jLat]+dLat,lon_[iLon]-dLon,
										lat_[jLat]+dLat,lon_[iLon]+dLon,
										lat_[jLat]-dLat,lon_[iLon]+dLon,
										r,g,b,
										r,g,b,
										r,g,b,
										r,g,b
								);
							}
						}

						if(nmqHTM) {
							SpatialVector v0, v1, v2, v3;
							v0.setLatLonDegrees(lat_[jLat]-dLat,lon_[iLon]-dLon);
							v1.setLatLonDegrees(lat_[jLat]-dLat,lon_[iLon]+dLon);
							v2.setLatLonDegrees(lat_[jLat]+dLat,lon_[iLon]+dLon);
							v3.setLatLonDegrees(lat_[jLat]+dLat,lon_[iLon]-dLon);
							RangeConvex rc = RangeConvex(&v0,&v1,&v2,&v3);
							SpatialDomain d; d.add(rc);
							HtmRange r;
							r.purge();
							bool overlap = d.intersect(&swathIndex,&r,swath_varlenHTMIDs);
							Key lo=-999, hi=-999;
							r.reset();
							int indexp = r.getNext(lo,hi);
							if(indexp) {
								do {
									nmqRangeCount++;
									nmqRange.addRange(lo,hi);
								} while(r.getNext(lo,hi));
							}
						}
					}
				}
			}
		}
	}
	preciprate_hsr_ifs.close();

	if(tapFileFlag) {
		tapFile.close();
	}

//	cout << 1900 << " nmqRangeCount: " << nmqRangeCount << endl << flush;
//	cout << 1901 << " los-stat:      " << endl << flush;
//	nmqRange.my_los->stat();
//	cout << 1902 << " his-stat:      " << endl << flush;
//	nmqRange.my_his->stat();
//	cout << 1904 << " los:           " << nmqRange.getLosLength() << endl << flush;
//	cout << 1905 << " his:           " << nmqRange.getHisLength() << endl << flush;
//
////	nmqRange.print(HtmRange::BOTH,cout);
//	cout << "lows" << endl;
//	nmqRange.print(HtmRange::LOWS,cout); cout << endl;
//	cout << "highs" << endl;
//	nmqRange.print(HtmRange::HIGHS,cout); cout << endl;
//
//	cout << 1910 << " nmqRange.nranges(): " << nmqRange.nranges() << endl << flush;

	// Dies in the following!

	// Plot the nmq triangles.
	if(nmqHTMRGeometry_viz){
		nmqRange.reset();
		////	viz->addHTMRange(&swathIndex,&nmqRange,1.0,0.5,0.0,-1,1.009);
		viz->addHTMRange(
				&swathIndex,
				&nmqRange,
				nmqHTMRGeometryColor.r,
				nmqHTMRGeometryColor.g,
				nmqHTMRGeometryColor.b,
				nmqHTMRGeometryColor.a,
				nmqHTMRGeometryColor.scale);
	}

//	cout << 2000 << endl << flush;

//	return;

	Key lo=-1, hi=-1;
	HtmRange red, blue, green;
	HtmRange *range1, *range2;
	SpatialIndex *index = &swathIndex;
//	SpatialIndex index(swathLevel,buildLevel);

//	nmqRange.defrag();
//	swathRange.defrag();

	// Note in the following we play games with the colors.
	red.purge(); blue.purge(); green.purge();
	if(nmqHTM) {
		range1 = &nmqRange;
	} else {
		range1 = &nmqRangeCoarse;
	}
	range2 = &swathRange;
	range1->reset();
	range2->reset();
	int indexp = range1->getNext(lo,hi);
	if(indexp) {
		do {
			int ret = range2->contains(lo,hi);
			if(ret == -1) { // Partial?
				red.addRange(lo,hi);
			} else if (ret == 0 ) { // Does not contain
				green.addRange(lo,hi);
			} else { // Contains
				blue.addRange(lo,hi);
			}
		} while(range1->getNext(lo,hi));

		if(nmq_trmm_HTMRPartial_viz){
			//		viz->addHTMRange(index,&red,0.7,0.,0.,-1,1.008);
//			viz->addHTMRange(index,&red,1.0,1.0,0.,-1,1.008); // error?
			viz->addHTMRange(index,&red,
					nmq_trmm_PartialIntersectionColor.r,
					nmq_trmm_PartialIntersectionColor.g,
					nmq_trmm_PartialIntersectionColor.b,
					nmq_trmm_PartialIntersectionColor.a,
					nmq_trmm_PartialIntersectionColor.scale
					);
		}
		if(nmq_trmm_HTMRDoesNotContain_viz){
			//		viz->addHTMRange(index,&green,0.,0.7,0.,-1,1.008);
			viz->addHTMRange(index,&green,
					nmqHTMRGeometryColor.r,
					nmqHTMRGeometryColor.g,
					nmqHTMRGeometryColor.b,
					nmqHTMRGeometryColor.a,
					nmqHTMRGeometryColor.scale
			);
		}

		if(nmq_trmm_HTMRContains_viz){
			//		viz->addHTMRange(index,&blue,0.,0.,0.7,-1,1.008);
			viz->addHTMRange(index,&blue, // intersection
					nmq_trmm_intersectionColor.r,
					nmq_trmm_intersectionColor.g,
					nmq_trmm_intersectionColor.b,
					nmq_trmm_intersectionColor.a,
					nmq_trmm_intersectionColor.scale
			);
		}
	}

//	return;

	red.purge(); blue.purge(); green.purge();
	range1 = &swathRange;
	if(nmqHTM) {
		range2 = &nmqRange;
	} else {
		range2 = &nmqRangeCoarse;
	}
	range1->reset();
	range2->reset();
	indexp = range1->getNext(lo,hi);
	if(indexp) {
		do {
			int ret = range2->contains(lo,hi);
			if(ret == -1) { // Partial?
				red.addRange(lo,hi);
			} else if (ret == 0 ) { // Does not contain
				green.addRange(lo,hi);
			} else { // Contains
				blue.addRange(lo,hi);
			}
		} while(range1->getNext(lo,hi));

		if(trmm_nmq_HTMRPartial_viz){
			//		viz->addHTMRange(index,&red,0.7,0.,0.,-1,1.008);
//			viz->addHTMRange(index,&red,1.0,1.0,0.,-1,1.008); // error?
			viz->addHTMRange(index,&red,
					trmm_nmq_PartialIntersectionColor.r,
					trmm_nmq_PartialIntersectionColor.g,
					trmm_nmq_PartialIntersectionColor.b,
					trmm_nmq_PartialIntersectionColor.a,
					trmm_nmq_PartialIntersectionColor.scale
					);
		}

		if(trmm_nmq_HTMRDoesNotContain_viz){
			//		viz->addHTMRange(index,&green,0.,0.7,0.,-1,1.008);
			viz->addHTMRange(index,&green,
					swathHTMRGeometryColor.r,
					swathHTMRGeometryColor.g,
					swathHTMRGeometryColor.b,
					swathHTMRGeometryColor.a,
					swathHTMRGeometryColor.scale
			);
		}

		if(trmm_nmq_HTMRContains_viz){
			// The following is redundant with the comparison above. Use as a check.
			//		viz->addHTMRange(index,&blue,0.,0.,0.7,-1,1.011);
			viz->addHTMRange(index,&blue, // intersection
					nmq_trmm_intersectionColor.r,
					nmq_trmm_intersectionColor.g,
					nmq_trmm_intersectionColor.b,
					nmq_trmm_intersectionColor.a,
					nmq_trmm_intersectionColor.scale
			);
		}
	}


	cout << "mn,mx: " << mn << " " << mx << endl << flush;
	std::cout << " nmq done " << endl << flush;

}


void testHstmRange(VizHTM *viz) {
	examiner_viz = true;
	offscreen_viz = false;

	viz->lineWidth = 1.5;
	viz->sphereComplexity = 1.0;
	//	plotBlockingSphere(viz,0.2,0.2,0.2,0.999);
	plotBlockingSphere(viz,0.2,0.2,0.2,0.99);
	testTenDegreeGridRGB(viz,0.6,0.6,0.6);

	EmbeddedLevelNameEncoding leftJustified;
	BitShiftNameEncoding rightJustified;
	HstmRange hstmRange;

	int indexLevel = 2, level;
	SpatialIndex *index = new SpatialIndex(indexLevel);

	uint64 id_, id0_, id1_;
	string name, name0, name1;
	HtmRange range;

	if(true){
		name = "N011";
		leftJustified.setName(name.c_str());
		id_ = leftJustified.getId();
		hstmRange.addRange(id_,id_);
		if(false) {
			level = leftJustified.getLevel();
			if( level != indexLevel ) {
				delete index;
				indexLevel = level;
				index = new SpatialIndex(indexLevel);
			}
			viz->addArcFromIndexAndName(index,name.c_str(),0.8,0.8,0.0,-1.0);
		}

		name0 = "N0120";
		leftJustified.setName(name0.c_str());
		id0_ = leftJustified.getId();
		name1 = "N0210";
		leftJustified.setName(name1.c_str());
		id1_ = leftJustified.getId();
		hstmRange.addRange(id0_,id1_);
		if(false){
			range.purge();
			range.addRange(rightJustified.idByName(name0.c_str()),rightJustified.idByName(name1.c_str()));
			level = leftJustified.getLevel();
			if( level != indexLevel ) {
				delete index;
				indexLevel = level;
				index = new SpatialIndex(indexLevel);
			}
			viz->addHTMRange(index,&range,0.0,0.8,0.8,-1.0);
		}
	}
	if(true) {
		name = "N01";
		leftJustified.setName(name.c_str());
		id_ = leftJustified.getId();
		hstmRange.addRange(id_,id_);
		if(false) {
			level = leftJustified.getLevel();
			if( level != indexLevel ) {
				delete index;
				indexLevel = level;
				index = new SpatialIndex(indexLevel);
			}
			viz->addArcFromIndexAndName(index,name.c_str(),0.5,0.9,0.5,-1.0);
		}
	}

	name0 = "S333300";
	leftJustified.setName(name0.c_str());
	id0_ = leftJustified.getId();
	name1 = "N0303333";
	leftJustified.setName(name1.c_str());
	id1_ = leftJustified.getId();
	hstmRange.addRange(id0_,id1_);
	if(false) {
		range.purge();
		range.addRange(rightJustified.idByName(name0.c_str()),rightJustified.idByName(name1.c_str()));
		level = leftJustified.getLevel();
		if( level != indexLevel ) {
			delete index;
			indexLevel = level;
			index = new SpatialIndex(indexLevel);
		}
		viz->addHTMRange(index,&range,0.0,0.8,0.8,-1.0);
	}

	viz->addHstmRange(&hstmRange,0.2,1.0,0.2,-1.0);
}

void testMultiResolution0(VizHTM *viz) {
	examiner_viz = true;
	offscreen_viz = false;

	viz->lineWidth = 1.5;
	viz->sphereComplexity = 1.0;
	//	plotBlockingSphere(viz,0.2,0.2,0.2,0.999);
	plotBlockingSphere(viz,0.2,0.2,0.2,0.99);
	testTenDegreeGridRGB(viz,0.6,0.6,0.6);

	testShapeFiles(viz);

	EmbeddedLevelNameEncoding leftJustified;
	BitShiftNameEncoding rightJustified;
	HstmRange hstmRange;

	int indexLevel = 4, level;
	int buildLevel = 8;
	SpatialIndex *index = new SpatialIndex(indexLevel,buildLevel);

	uint64 id_, id0_, id1_;
	string name, name0, name1;
	HtmRange range, interior, boundary;
	KeyPair kp;
	int indexp, k;

	SpatialVector *v0 = VectorFromLatLonDegrees(10.0,0.0);
	SpatialVector *v1 = VectorFromLatLonDegrees(70.0,0.0);
	SpatialVector *v2 = VectorFromLatLonDegrees(70.0,80.0);
	SpatialVector *v3 = VectorFromLatLonDegrees(10.0,80.0);

	viz->addRectangle(*v0,*v1,*v2,*v3,1.0,0.,0.);
//	viz->addArcAtLatitudeDegrees(10.0,0.,80.0,1.0,0.,0.);
//	viz->addArcAtLatitudeDegrees(70.0,0.,80.0,1.0,0.,0.);

	RangeConvex *rc = new RangeConvex(v0,v1,v2,v3);
	SpatialDomain domain = SpatialDomain(index);
	domain.add(*rc);

	bool varlen_individualHTMIds = false; // true for individuals, false for ranges

	range.purge(); interior.purge(); boundary.purge();
	bool overlap = domain.intersect(index,&range,varlen_individualHTMIds,&interior,&boundary);

//	while((indexp = range.getNext(kp)) >= 0) {
//		hstmRange.addRange(kp.lo,kp.hi);
//	}
//	viz->addHstmRange(&hstmRange,0.2,1.0,0.2,-1.0);

//	viz->addHTMRange(&interior,0.0,1.0,0.0);
//	viz->addHTMRange(&boundary,1.0,0.0,0.0);

////	BitShiftNameEncoding right;
//	cout << 1000 << endl << flush;
	hstmRange.purge();
//	interior.reset();
//	k = 0;
//	while((indexp = interior.getNext(kp)) > 0) {
//		Key l = rightJustified.leftJustifiedId(kp.lo);
//		Key h = rightJustified.leftJustifiedId(kp.hi);
////		cout << 200 << " lh " << hex << l << " " << h << dec << endl << flush;
////		cout << 201 << " kp " << hex << kp.lo << " " << kp.hi << dec << endl << flush;
//		hstmRange.addRange(l,h);
//		k++;
////		if(k>0)break;
//	}

//	viz->addHstmRange(&hstmRange,0.0,1.0,1.0);
//	return;
//	cout << 2000 << endl << flush;

//	for( indexLevel = 6; indexLevel >= 5; indexLevel-- ) {
	//	for( indexLevel = 3; indexLevel >= 2; indexLevel-- ) {
	for( indexLevel = 8; indexLevel >= 2; indexLevel-- ) {
//	for( indexLevel = 5; indexLevel <= 6; indexLevel++ ) {
		cout << "Constructing level - indexLevel: " << indexLevel << endl << flush;
		range.purge(); interior.purge(); boundary.purge();
		delete index;
		index = new SpatialIndex(indexLevel,buildLevel);
		domain = SpatialDomain(index);
		domain.add(*rc);
		overlap = domain.intersect(index,&range,varlen_individualHTMIds,&interior,&boundary);

		interior.reset();
		k = 0;
		while((indexp = interior.getNext(kp)) > 0) {
			Key l = rightJustified.leftJustifiedId(kp.lo);
			Key h = rightJustified.leftJustifiedId(kp.hi);
//#define hexOut(a,b,c) cout << a << " 0x" << hex << setfill('0') << setw(16) << b << ".." << c << dec << endl << flush;
//			hexOut("lh  ",l,h);
//#undef hexOut
			hstmRange.addRange(l,h);
			k++;
			if(false) if( hstmRange.range->my_los->getLength() != hstmRange.range->my_his->getLength() ) {
				cout << "addRange error at iteration " << k << endl << flush;
				cout << "nranges "
						<< hstmRange.range->my_los->getLength() << " "
						<< hstmRange.range->my_his->getLength() << endl << flush;
				cout << " lh " << hex << l << " " << h << dec << endl << flush;
				cout << " lh1 " << hex << leftJustified.maskOffLevelBit(l) << " "
						<< leftJustified.maskOffLevelBit(h)  << dec << endl << flush;
				cout << " kp " << hex << kp.lo << " " << kp.hi << dec << endl << flush;
				cout << " lhl " << leftJustified.nameById(l) << " "
					           << leftJustified.nameById(h) << endl << flush;
				cout << " lhr " << rightJustified.nameById(kp.lo) << " "
						        << rightJustified.nameById(kp.hi) << endl << flush;
				int kount = 0;
				hstmRange.reset();
				while( (indexp = hstmRange.getNext(kp)) > 0 ) {
					cout << kount << " kp 0x" << hex << kp.lo << ", 0x" << kp.hi << ",  " << dec
							<< endl << flush;
					kount++;
				}
				cout << "lows" << endl << flush;
				hstmRange.range->print( HtmRangeMultiLevel_NameSpace::HtmRangeMultiLevel::LOWS,cout,true);
				cout << "highs" << endl << flush;
				hstmRange.range->print( HtmRangeMultiLevel_NameSpace::HtmRangeMultiLevel::HIGHS,cout,true);
				cout << " exiting " << endl << flush;
				exit(1);
			}
			//		if(k>0)break;
		}
		cout << "nranges "
				<< hstmRange.range->my_los->getLength() << " "
				<< hstmRange.range->my_his->getLength() << endl << flush;
		if(false) {
			hstmRange.reset();
			k = 0;
			cout << "Printing hstmRange" << endl << flush;
			while( (indexp = hstmRange.getNext(kp)) > 0 ) {
				cout << k << " kp 0x" << hex << kp.lo << ", 0x" << kp.hi << ",  " << dec
						<< endl << flush;
				k++;
			}
		}
	}

	cout << "constructing hstm visualization" << endl << flush;
//	viz->addHstmRange(&hstmRange,0.0,1.0,1.0);
//	viz->addHTMRange(&interior,0.0,1.0,0.0);
//	viz->addHTMRange(&boundary,1.0,0.0,0.0);
//	return;

	range.purge(); interior.purge(); boundary.purge();
	indexLevel = 8;
	delete index;
	index = new SpatialIndex(indexLevel,buildLevel);
	domain = SpatialDomain(index);
	domain.add(*rc);
	overlap = domain.intersect(index,&range,varlen_individualHTMIds,&interior,&boundary);

	boundary.reset();
	k = 0;
	while((indexp = boundary.getNext(kp)) > 0) {
		Key l = rightJustified.leftJustifiedId(kp.lo);
		Key h = rightJustified.leftJustifiedId(kp.hi);
//		cout << 200 << " lh " << hex << l << " " << h << dec << endl << flush;
//		cout << 201 << " kp " << hex << kp.lo << " " << kp.hi << dec << endl << flush;
		hstmRange.addRange(l,h);
		k++;
//		if(k>0)break;
	}

	viz->addHstmRange(&hstmRange,0.0,1.0,1.0);

}

void testMultiResolution(VizHTM *viz) {
	examiner_viz = true;
	offscreen_viz = false;

	viz->lineWidth = 1.5;
	viz->sphereComplexity = 1.0;
	//	plotBlockingSphere(viz,0.2,0.2,0.2,0.999);
	plotBlockingSphere(viz,0.2,0.2,0.2,0.99);
	testTenDegreeGridRGB(viz,0.6,0.6,0.6);

	testShapeFiles(viz);

	EmbeddedLevelNameEncoding leftJustified;
	BitShiftNameEncoding rightJustified;
	HstmRange hstmRange;

	int indexLevel = 4, level;
	int buildLevel = 8;
	SpatialIndex *index = new SpatialIndex(indexLevel,buildLevel);

	uint64 id_, id0_, id1_;
	string name, name0, name1;
	HtmRange range, interior, boundary;
	KeyPair kp;
	int indexp, k;

	SpatialVector *v0 = VectorFromLatLonDegrees(10.0,0.0);
	SpatialVector *v1 = VectorFromLatLonDegrees(70.0,0.0);
	SpatialVector *v2 = VectorFromLatLonDegrees(70.0,80.0);
	SpatialVector *v3 = VectorFromLatLonDegrees(10.0,80.0);

	viz->addRectangle(*v0,*v1,*v2,*v3,1.0,0.,0.);
	//	viz->addArcAtLatitudeDegrees(10.0,0.,80.0,1.0,0.,0.);
	//	viz->addArcAtLatitudeDegrees(70.0,0.,80.0,1.0,0.,0.);

	RangeConvex *rc = new RangeConvex(v0,v1,v2,v3);
	SpatialDomain domain = SpatialDomain(index);
	domain.add(*rc);

	SpatialVector *v_ = VectorFromLatLonDegrees(10.0,0.0);
	SpatialVector v = *v_;
	v = v*(-1);
	SpatialConstraint sc(v,0.95);
    sc.invert();
    SpatialDomain include = SpatialDomain(index);
	RangeConvex *rc1 = new RangeConvex; rc1->add(sc);
	include.add(*rc1);
	HtmRange includeIntersectRange, includeInteriorRange, includeBoundaryRange;

	bool varlen_individualHTMIds = false; // true for individuals, false for ranges

	includeIntersectRange.purge(); includeInteriorRange.purge(); includeBoundaryRange.purge();
	bool overlapInclude = domain.intersect(
			index,
			&includeIntersectRange,varlen_individualHTMIds,
			&includeInteriorRange, &includeBoundaryRange );

	range.purge(); interior.purge(); boundary.purge();
	bool overlap = domain.intersect(index,&range,varlen_individualHTMIds,&interior,&boundary);

	//	while((indexp = range.getNext(kp)) >= 0) {
	//		hstmRange.addRange(kp.lo,kp.hi);
	//	}
	//	viz->addHstmRange(&hstmRange,0.2,1.0,0.2,-1.0);

	//	viz->addHTMRange(&interior,0.0,1.0,0.0);
	//	viz->addHTMRange(&boundary,1.0,0.0,0.0);

	////	BitShiftNameEncoding right;
	//	cout << 1000 << endl << flush;
	hstmRange.purge();
	//	interior.reset();
	//	k = 0;
	//	while((indexp = interior.getNext(kp)) > 0) {
	//		Key l = rightJustified.leftJustifiedId(kp.lo);
	//		Key h = rightJustified.leftJustifiedId(kp.hi);
	////		cout << 200 << " lh " << hex << l << " " << h << dec << endl << flush;
	////		cout << 201 << " kp " << hex << kp.lo << " " << kp.hi << dec << endl << flush;
	//		hstmRange.addRange(l,h);
	//		k++;
	////		if(k>0)break;
	//	}

	//	viz->addHstmRange(&hstmRange,0.0,1.0,1.0);
	//	return;
	//	cout << 2000 << endl << flush;


	int totalTrianglesAdded = 0;

//	vector<RangeConvex> geometry;
//	geometry.push_back(*rc);
//	geometry.push_back(*rc1);

	//	for( indexLevel = 6; indexLevel >= 5; indexLevel-- ) {
	//	for( indexLevel = 3; indexLevel >= 2; indexLevel-- ) {
//	for( indexLevel = 6; indexLevel <= 6; indexLevel++ ) { // Forward -- Try interior subtraction
	for( indexLevel = 2; indexLevel <= 8; indexLevel++ ) { // Forward -- Try interior subtraction
		//		for( indexLevel = 2; indexLevel <= 8; indexLevel++ ) { // Forward -- Try interior subtraction
		//	 for( indexLevel = 8; indexLevel >= 4; indexLevel-- ) {  // Good w/o interior subtraction.
		//	for( indexLevel = 5; indexLevel <= 6; indexLevel++ ) {
		cout << "Constructing level - indexLevel: " << indexLevel << endl << flush;
		range.purge(); interior.purge(); boundary.purge();
		delete index;
		index = new SpatialIndex(indexLevel,buildLevel);
		domain = SpatialDomain(index);
		domain.add(*rc);
//		for(vector<RangeConvex>::iterator gIter = geometry.begin(); gIter != geometry.end(); ++gIter) {
//			domain.add(*gIter);
//		}
		overlap = domain.intersect(index,&range,varlen_individualHTMIds,&interior,&boundary);

		if(false) {
			includeIntersectRange.purge(); includeInteriorRange.purge(); includeBoundaryRange.purge();
			include = SpatialDomain(index);
			include.add(*rc1);
			bool overlapInclude = include.intersect(
					index,
					&includeIntersectRange,varlen_individualHTMIds,
					&includeInteriorRange, &includeBoundaryRange );
			//
			HtmRange newRange; newRange.purge(); KeyPair p0, p1;
			interior.reset();
			while((indexp = interior.getNext(p0)) > 0) {
				//			cout << "p0: " << p0.lo << " " << p0.hi << " " << flush;
				p1 = includeInteriorRange.intersection(p0);
				//			cout << "p1: " << p1.lo << " " << p1.hi << " " << endl << flush;
				if(p1.set) {
					newRange.addRange(p1.lo,p1.hi);
				}
			}
		}

		HtmRange *addRange = &interior;
//		addRange = &includeInteriorRange;
//		addRange = &newRange;

		cout << "Adding " << addRange->nranges() << " triangles at level " << indexLevel << endl << flush;
		totalTrianglesAdded += addRange->nranges() ;
		addRange->reset();
		k = 0;
		while((indexp = addRange->getNext(kp)) > 0) {
			Key l = rightJustified.leftJustifiedId(kp.lo);
			Key h = rightJustified.leftJustifiedId(kp.hi);
			//#define hexOut(a,b,c) cout << a << " 0x" << hex << setfill('0') << setw(16) << b << ".." << c << dec << endl << flush;
			//			hexOut("lh  ",l,h);
			//#undef hexOut
			hstmRange.addRange(l,h);
			k++;
			if(false) if( hstmRange.range->my_los->getLength() != hstmRange.range->my_his->getLength() ) {
				cout << "addRange error at iteration " << k << endl << flush;
				cout << "nranges "
						<< hstmRange.range->my_los->getLength() << " "
						<< hstmRange.range->my_his->getLength() << endl << flush;
				cout << " lh " << hex << l << " " << h << dec << endl << flush;
				cout << " lh1 " << hex << leftJustified.maskOffLevelBit(l) << " "
						<< leftJustified.maskOffLevelBit(h)  << dec << endl << flush;
				cout << " kp " << hex << kp.lo << " " << kp.hi << dec << endl << flush;
				cout << " lhl " << leftJustified.nameById(l) << " "
						<< leftJustified.nameById(h) << endl << flush;
				cout << " lhr " << rightJustified.nameById(kp.lo) << " "
						<< rightJustified.nameById(kp.hi) << endl << flush;
				int kount = 0;
				hstmRange.reset();
				while( (indexp = hstmRange.getNext(kp)) > 0 ) {
					cout << kount << " kp 0x" << hex << kp.lo << ", 0x" << kp.hi << ",  " << dec
							<< endl << flush;
					kount++;
				}
				cout << "lows" << endl << flush;
				hstmRange.range->print( HtmRangeMultiLevel_NameSpace::HtmRangeMultiLevel::LOWS,cout,true);
				cout << "highs" << endl << flush;
				hstmRange.range->print( HtmRangeMultiLevel_NameSpace::HtmRangeMultiLevel::HIGHS,cout,true);
				cout << " exiting " << endl << flush;
				exit(1);
			}
			//			// Try to remove interior.
			// Broken...
			//			if(indexLevel < 4) {
			//				SpatialVector a,b,c;
			//				for( uint64 numericId = kp.lo; numericId <= kp.hi; numericId ++ ) {
			//					uint64 nodeIndex = index->nodeIndexFromId(numericId);
			//					index->nodeVertex(nodeIndex,a,b,c);
			//					RangeConvex r = RangeConvex(&a,&b,&c);
			//					r.invert();
			//					geometry.push_back(r);
			//				}
		}



		//		if(k>0)break;
	}
	cout << "nranges "
			<< hstmRange.range->my_los->getLength() << " "
			<< hstmRange.range->my_his->getLength() << endl << flush;
	if(false) {
		hstmRange.reset();
		k = 0;
		cout << "Printing hstmRange" << endl << flush;
		while( (indexp = hstmRange.getNext(kp)) > 0 ) {
			cout << k << " kp 0x" << hex << kp.lo << ", 0x" << kp.hi << ",  " << dec
					<< endl << flush;
			k++;
		}
	}

	cout << "Total triangles added: " << totalTrianglesAdded << endl << flush;
	cout << "HstmRange intervals:   " << hstmRange.range->nranges() << endl << flush;

	cout << "constructing hstm visualization" << endl << flush;
	//	viz->addHstmRange(&hstmRange,0.0,1.0,1.0);
	//	viz->addHTMRange(&interior,0.0,1.0,0.0);
	//	viz->addHTMRange(&boundary,1.0,0.0,0.0);
	//	return;

	range.purge(); interior.purge(); boundary.purge();
	indexLevel = 8;
	delete index;
	index = new SpatialIndex(indexLevel,buildLevel);
	domain = SpatialDomain(index);
	domain.add(*rc);
	overlap = domain.intersect(index,&range,varlen_individualHTMIds,&interior,&boundary);

	boundary.reset();
	k = 0;
	while((indexp = boundary.getNext(kp)) > 0) {
		Key l = rightJustified.leftJustifiedId(kp.lo);
		Key h = rightJustified.leftJustifiedId(kp.hi);
		//		cout << 200 << " lh " << hex << l << " " << h << dec << endl << flush;
		//		cout << 201 << " kp " << hex << kp.lo << " " << kp.hi << dec << endl << flush;
		hstmRange.addRange(l,h);
		k++;
		//		if(k>0)break;
	}

	viz->addHstmRange(&hstmRange,0.0,1.0,1.0,-1.0,1.1);

}

void testTransparency(VizHTM *viz) {

	// On Coin3D on Mac OS X
	// DELAYED_BLEND crashes on the following
	// BLEND seems to work. Looks better than ADD.
	// ADD works.
	// SORTED_* crashes

	viz->addFace3(
			0.0,0.0,0.0,
			0.0,1.0,0.0,
			1.0,1.0,0.0,
			1.0,0.0,0.0,
			1.0,0.0,0.0,
			1.0,0.0,0.0,
			0.0,0.0,0.0
			);

	viz->addFace3(
			0.0,0.0,0.15,
			0.0,1.0,0.15,
			1.0,1.0,0.15,
			0.0,1.0,0.0,
			0.0,1.0,0.0,
			0.0,1.0,0.0,
			0.5,0.5,0.5
			);

	viz->addFace3(
			0.0,0.0,0.3,
			0.0,1.0,0.3,
			1.0,1.0,0.3,
			0.0,0.0,1.0,
			0.0,0.0,1.0,
			0.0,0.0,1.0,
			0.6,0.6,0.6
			);


}

void addFiduciaries(VizHTM *viz) {

	string htmName = "N00";
	SpatialIndex index = SpatialIndex(levelOfName(htmName.c_str()));
	viz->addArcFromIndexAndName(
			&index,
			htmName.c_str(),
			1.0,0.0,0.0
	);

	htmName = "N01";
	index = SpatialIndex(levelOfName(htmName.c_str()));
	viz->addArcFromIndexAndName(
			&index,
			htmName.c_str(),
			0.0,1.0,0.0
	);

	htmName = "N02";
	index = SpatialIndex(levelOfName(htmName.c_str()));
	viz->addArcFromIndexAndName(
			&index,
			htmName.c_str(),
			0.0,0.0,1.0
	);

//	htmName = "N03";
//	index = SpatialIndex(levelOfName(htmName.c_str()));
//	viz->addArcFromIndexAndName(
//			&index,
//			htmName.c_str(),
//			1.0,1.0,1.0
//	);

	htmName = "N10";
	index = SpatialIndex(levelOfName(htmName.c_str()));
	viz->addArcFromIndexAndName(
			&index,
			htmName.c_str(),
			1.0,0.0,0.0
	);

	htmName = "N11";
	index = SpatialIndex(levelOfName(htmName.c_str()));
	viz->addArcFromIndexAndName(
			&index,
			htmName.c_str(),
			0.0,1.0,0.0
	);

	htmName = "N12";
	index = SpatialIndex(levelOfName(htmName.c_str()));
	viz->addArcFromIndexAndName(
			&index,
			htmName.c_str(),
			0.0,0.0,1.0
	);


}

vector< string > cellsFiles;
vector< string > cellsDataRange;
vector< string > cellsGridColor;
float cellsGridRed, cellsGridGreen, cellsGridBlue;
bool cellsPlotFace = false;
vector< string > cellsScale;
vector< string > cellsCentroids;

void plotCells(VizHTM *viz) {

	float dataLo=-1, dataHi=-1;
	float scale = 1.001;

	vector< string >::iterator cellsDataRangeIter = cellsDataRange.begin();
	vector< string >::iterator cellsGridColorIter = cellsGridColor.begin();
	vector< string >::iterator cellsScaleIter     = cellsScale.begin();
	vector< string >::iterator cellsCentroidsIter = cellsCentroids.begin();

	float64 lon[4],lat[4],data;
	SpatialVector v[4];

	SpatialVector centroid;
	float centroidH, centroidRadius;
	float centroidR, centroidG, centroidB;
	bool centroidFlag;

	for(    vector< string >::iterator cellsFilesIter = cellsFiles.begin();
			cellsFilesIter != cellsFiles.end();
			++cellsFilesIter ) {
		string fileName = (*cellsFilesIter);
		cout << "plotCells reading " << fileName << endl << flush;
		if( cellsDataRangeIter != cellsDataRange.end() ) {
			if(strcmp((*cellsDataRangeIter).c_str(),"default") != 0) {
				sscanf((*cellsDataRangeIter).c_str(),"%f,%f",&dataLo,&dataHi);
//				cout << " color mapping data range: " << dataLo << " " << dataHi << endl << flush;
				cellsPlotFace = true;
			} else {
//				cout << " using default color mapping" << endl << flush;
				cellsPlotFace = false;
			}
		} else {
			cellsPlotFace = false;
		}
		if( cellsGridColorIter != cellsGridColor.end() ) {
			if(strcmp((*cellsGridColorIter).c_str(),"default") != 0) {
				sscanf((*cellsGridColorIter).c_str(),"%f,%f,%f",
						&cellsGridRed,&cellsGridGreen, &cellsGridBlue);
			} else {
				cellsGridRed = 0.7; cellsGridGreen = 0.7; cellsGridBlue = 0.7;
			}
		} else {
			cellsGridRed = 0.7; cellsGridGreen = 0.7; cellsGridBlue = 0.7;
		}
		if( cellsScaleIter != cellsScale.end() ) {
			if(strcmp((*cellsScaleIter).c_str(),"default") != 0) {
				sscanf((*cellsScaleIter).c_str(),"%f",&scale);
			} else {
				scale += 0.001;
			}
		} else {
			scale += 0.001;
		}
		// H = height; Radius = radius; RGB = color.
		if( cellsCentroidsIter != cellsCentroids.end() ) {
			if(strcmp((*cellsCentroidsIter).c_str(),"default") != 0) {
				sscanf((*cellsCentroidsIter).c_str(),"%f,%f,%f,%f,%f",
						&centroidH, &centroidRadius,
						&centroidR, &centroidG, &centroidB);
				centroidFlag = true;
			} else {
				centroidFlag = false;
			}
		} else {
			centroidFlag = false;
		}

		float r,g,b;
		ifstream cellsIn;
		cellsIn.open(fileName);
		int uniqueHTMs = 0;
		for(string line; getline(cellsIn, line); ) {
//			cout << "line: " << line << endl << flush;
			int64 nPoints;
			float64 data;
			char line1[2048];
			sscanf(line.c_str(),"%llu,%s",&nPoints,&line1);
			sscanf(line1,"%lf,%s",&data,&line1);
//			cout << "n: " << nPoints
//					<< " d: " << data
//					<< " ( "
//					<< flush;
			int ic = 0;
			for(int ip=0; ip<nPoints; ++ip) {
				sscanf(line1,"%lf,%lf,%s",&lon[ip],&lat[ip],&line1);
//				cout << " " << x[ip] << " " << y[ip] << flush;
				v[ip].setLatLonDegrees(lat[ip],lon[ip]);
			}
//			cout << " ) "<< endl << flush;
			if(cellsPlotFace) {
				if(nPoints == 4){
//					cout << "adding point..." << endl << flush;
					float a[4]; for(int i = 0; i < 4; i++) { a[i] = 0.0; }
					colorScheme1_rgb(data,dataLo,dataHi,r,g,b);
					viz->addFace4FromLatLonDegrees(
							lat[0],lon[0],
							lat[1],lon[1],
							lat[2],lon[2],
							lat[3],lon[3],
							r,g,b,
							r,g,b,
							r,g,b,
							r,g,b,
							a[0],a[1],a[2],a[3],
							scale
							);
				} // TODOO add support for other shapes (nPoints != 4).
			}
			for(int ip=0; ip < nPoints; ++ip) {
				v[ip].setLatLonDegrees(lat[ip],lon[ip]);
			}
			{
				float
				r = cellsGridRed,
				g = cellsGridGreen,
				b = cellsGridBlue;
			viz->addArc(v[0],v[1],r,g,b,0.0,scale,20);
			viz->addArc(v[1],v[2],r,g,b,0.0,scale,20);
			viz->addArc(v[2],v[3],r,g,b,0.0,scale,20);
			viz->addArc(v[3],v[0],r,g,b,0.0,scale,20);
			}
			if(centroidFlag) {
				if(nPoints == 4) {
					centroid = v[0];
					for(int ip=1; ip < 4; ++ip ) {
						centroid = centroid + v[ip];
					}
					centroid *= 0.25*centroidH;
//					cout << "centroid: " << centroid << endl << flush;
					if(centroidR>0){
						viz->addSphere(
								centroid,
								centroidR, centroidG, centroidB,
								centroidRadius);
					} else {
						colorScheme1_rgb(data,centroidG,centroidB,r,g,b);
						viz->addSphere(
								centroid,
								r, g, b,  // from data...
								centroidRadius);
					}
				}
			}
		}
		if(cellsDataRangeIter != cellsDataRange.end()) {
			++cellsDataRangeIter;
		}
		if(cellsGridColorIter != cellsGridColor.end()) {
			++cellsGridColorIter;
		}
		if(cellsScaleIter != cellsScale.end()) {
			++cellsScaleIter;
		}
		if(cellsCentroidsIter != cellsCentroids.end()) {
			++cellsCentroidsIter;
		}
		cellsIn.close();
	}
	cout << "plotCells done" << endl << flush;
}



vector< string > csvNames;
vector< string > csvLevels; // TODO maybe get some class
vector< string > csvData;
vector< string > csvGridColor;

void plotCsv(VizHTM *viz) {

	SpatialVector *tmpV = new SpatialVector();

	EmbeddedLevelNameEncoding leftJustified;
	BitShiftNameEncoding      rightJustified;
	int count = 0;

	vector< string >::iterator csvLevelIter = csvLevels.begin();
	int csvLevel = -1;
	SpatialIndex csvIndex;

	vector< string >::iterator csvDataIter = csvData.begin();
	bool csvDataFlag = false;
	float csvDataLo = -1, csvDataHi = -1;

	vector< string >::iterator csvGridColorIter = csvGridColor.begin();
	float csvGridRed=0.7, csvGridGreen=0.7, csvGridBlue=0.7;

	for(vector< string >::iterator csvIter = csvNames.begin();
			csvIter != csvNames.end();
			++csvIter) {
		int64 idLast = -1;
		++count;
		string fileName = (*csvIter);
		cout << "plotCsv reading " << fileName << "..." << endl << flush;

		if(csvLevelIter != csvLevels.end()) {
			cout << "csvLevel string: " << (*csvLevelIter).c_str() << endl;
			if(strcmp((*csvLevelIter).c_str(),"default") != 0) {
				sscanf((*csvLevelIter).c_str(),"%d",&csvLevel);
				csvIndex = SpatialIndex(csvLevel,5);
				cout << "coercing level to " << csvLevel << endl;
			} else {
				cout << "using default level" << endl;
				csvLevel = -1;
			}
		} else {
			csvLevel = -1;
		}

		if(csvDataIter != csvData.end()) {
			if(strcmp((*csvDataIter).c_str(),"default") != 0) {
				csvDataFlag = true;
				cout << "data bounds: " << (*csvDataIter).c_str() << endl;
				sscanf((*csvDataIter).c_str(),"%f,%f",&csvDataLo,&csvDataHi);
				cout << "data bounds: " << csvDataLo << " " << csvDataHi << endl;
			}
		}

		if(csvGridColorIter != csvGridColor.end()) {
			if(strcmp((*csvGridColorIter).c_str(),"default") != 0) {
				sscanf((*csvGridColorIter).c_str(),"%f,%f,%f",
						&csvGridRed, &csvGridGreen, &csvGridBlue);
			} else {
				csvGridRed=0.0; csvGridGreen=0.7; csvGridBlue=0.0;
			}
		} else {
			csvGridRed=0.0; csvGridGreen=0.7; csvGridBlue=0.0;
		}

		ifstream csvIn;
		csvIn.open(fileName);
		int uniqueHTMs = 0;
		for(string line; getline(csvIn, line); ) {
			int64 iHtm;
			int64 iCell;
			float64 data;
			sscanf(line.c_str(),"%llu,%llu,%lf",&iHtm,&iCell,&data);
//			cout << "... " << iHtm << " " << iCell << " " << data << " "; // << endl << flush;
			leftJustified.setIdFromSciDBLeftJustifiedFormat(iHtm);
			rightJustified.setId(leftJustified.rightJustifiedId());
//			cout << "+ " << endl << flush;
//			cout << "+ " << leftJustified.getName() << endl << flush;
//			cout << "- " << rightJustified.getName() << flush;
			int level = leftJustified.getLevel();

//			cout << " l= " << level << endl << flush;
			SpatialIndex index(level,5);

			if(idLast != rightJustified.getId()) {
				idLast = rightJustified.getId();
				++uniqueHTMs;
				HtmRange *range = new HtmRange;
				range->purge();
				if(csvLevel==-1) {
					range->addRange(idLast,idLast);
				} else {
					index.pointByHtmId(*tmpV,idLast);
					int64 tId = csvIndex.idByPoint(*tmpV);
					range->addRange(tId,tId);
				}
				range->reset();
				double r = 1.0,g = 0.8,b = 0.3,a = 0.0,scale = 1.000;
				if(count==2) {
					a = 0.0;
					scale = 1.0015;
				} else if(count==3) {
					a = 0.0;
					scale = 1.003;
				}
				r = csvGridRed;
				g = csvGridGreen;
				b = csvGridBlue;

				if(csvLevel==-1) {
					viz->addHTMRange(
							&index,
							range,
							r,
							g,
							b,
							a,
							scale
					);
				} else {
					viz->addHTMRange(
							&csvIndex,
							range,
							r,
							g,
							b,
							a,
							scale
					);
				}

				if(csvDataFlag) {
					Key lo, hi;
					float r0, g0, b0;
					float r1, g1, b1;
					float r2, g2, b2;

					range->reset(); range->getNext(lo,hi);
					colorScheme1_rgb(data,csvDataLo,csvDataHi,r0,g0,b0);
					colorScheme1_rgb(data,csvDataLo,csvDataHi,r1,g1,b1);
					colorScheme1_rgb(data,csvDataLo,csvDataHi,r2,g2,b2);
					if(csvLevel==-1) {
						viz->addFaceFromIndexAndId(
								&index,
								lo,
								r0, g0, b0,
								r1, g1, b1,
								r2, g2, b2,
								a, scale
						);
					} else {
						viz->addFaceFromIndexAndId(
								&csvIndex,
								lo,
								r0, g0, b0,
								r1, g1, b1,
								r2, g2, b2,
								a, scale
						);
					}
				}


			}
		}
		if(csvLevelIter != csvLevels.end()) {
			++csvLevelIter;
		}
		if(csvDataIter != csvData.end()) {
			++csvDataIter;
		}
		if(csvGridColorIter != csvGridColor.end()) {
			++csvGridColorIter;
		}
		csvIn.close();
		cout << "plotCsv " << fileName << " done added " << uniqueHTMs << " unique HTMs" << endl << flush;
	} // csvIter
	cout << "plotCsv done" << endl << flush;

	delete tmpV;
}

string lookFromArgs;
void setLookFrom(SoSeparator *root, SoSeparator *content, SbViewportRegion *vpRegion) {

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



void loadScene(SoSeparator *root, SoSeparator *content, SbViewportRegion *vpRegion){

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
    SpatialVector p = (*nmq_trmm_vizCenter);
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



int main(int argc, char *argv[]) {

	const char* mainName = "VizHTM-main";

	int c;
	static int verbose_flag=0, help_flag=0, quiet_flag=0;
	static int
		testLevelChunk_flag=0,
		testMultiResolutionHtmRange_flag=0,
		testTriaxis_flag=0,
		testAddRectangle_flag=0,
		testTwoRectangle_flag=0,
		testLatLonSpiral_flag=0,
		testGeorgia_flag=0,
		testPlotDataSetIntersection_flag=0,
		testPlotDataSetIntersectionRangeContains_flag=0,
		testPlotDataSetIntersectionNativeIntersect_flag=0,
		testTenDegreeGrid_flag=0,
		testHTMRangeToLevel15_flag=0,
		testSwathNMQ_flag=0,
		testShapeFiles_flag=0,
		testChunks_flag=0,
		blockingSphere_flag=0,
		shapeFilesBlockingSphereScale_flag=0,
		lookFrom_flag=0;

	static int
		level_     = 0,
		saveLevel_ = 0;

	vector<string> htmNames; htmNames.clear();

	char *cvalue = NULL;
	// TODO Make all of this into classified code.
	static struct option long_options[] =
	{
			// { "name", no/required _argument, flag to set, ShortForm

			// These options set flags.
			// {"verbose", no_argument, &verbose_flag, 1},
			{"help",    no_argument, &help_flag,                  1},
			{"quiet",   no_argument, &quiet_flag,                 2},
			{"testLevelChunk", no_argument, &testLevelChunk_flag, 3},
			{"testMultiResolutionHtmRange", no_argument, &testMultiResolutionHtmRange_flag, 4},
			{"testTriaxis", no_argument, &testTriaxis_flag, 5},
			{"testAddRectangle", no_argument, &testAddRectangle_flag, 6},
			{"testTwoRectangle", no_argument, &testTwoRectangle_flag, 7},
			{"testLatLonSpiral", no_argument, &testLatLonSpiral_flag, 8},
			{"testGeorgia", no_argument, &testGeorgia_flag, 9},
			{"testPlotDataSetIntersection", no_argument, &testPlotDataSetIntersection_flag, 10},
			{"testPlotDataSetIntersectionRangeContains", no_argument, &testPlotDataSetIntersectionRangeContains_flag, 11},
			{"testPlotDataSetIntersectionNativeIntersect", no_argument, &testPlotDataSetIntersectionNativeIntersect_flag, 12},
			{"testTenDegreeGrid", no_argument, &testTenDegreeGrid_flag, 13},
			{"testHTMRangeToLevel15", no_argument, &testHTMRangeToLevel15_flag, 14},

			{"testSwathNMQ", no_argument, &testSwathNMQ_flag, 15},

			{"testChunks", no_argument, &testChunks_flag, 16},
			{"blockingSphere", no_argument, &blockingSphere_flag, 17},

			// These options are not setters.
			// {"add",     no_argument,       0, 'a'},
			// {"append",  no_argument,       0, 'b'},
			// {"delete",  required_argument, 0, 'd'},
			// {"create",  required_argument, 0, 'c'},
			// {"file",    required_argument, 0, 'f'},
			// corresponds to parsing pattern "abc:d:f:"
			// optarg has the data

			{"level",      required_argument, 0, 1000},
			{"saveLevel",  required_argument, 0, 1001},
			{"htmName",   required_argument, 0, 1002},

			{"csv",      required_argument, 0, 1010},
			{"csvLevel", required_argument, 0, 1011},
			{"csvData",  required_argument, 0, 1012},
			{"csvGridColor", required_argument, 0, 1013},

			{"cells", required_argument, 0, 1020},
			{"cellsDataRange", required_argument, 0, 1021},
			{"cellsGridColor", required_argument, 0, 1022},
			{"cellsScale", required_argument, 0, 1023},
			{"cellsCentroid", required_argument, 0, 1024},

			{"testShapeFiles", required_argument, 0, 1100},
			{"shapeFilesBlockingSphereScale", required_argument, 0, 1101},

			{"lineWidth", required_argument, 0, 2000},

			{"lookFrom", required_argument, 0, 3000},

			{"test0",     required_argument,  0, 0},
			{"test",      required_argument,  0, 't'},

			{0, 0, 0, 0}
	};
	string tmpStr;
	size_t nameSize;
	stringstream ss;
	while (1) {
		int option_index = 0;
		c = getopt_long( argc, argv, "t:",
				long_options, &option_index);
		if (c == -1) { break; }

		switch (c) {
		case 0: // Arguments that do not have a short form.
			if(long_options[option_index].flag != 0) {
				break; // A flag was set.
			}
//			cout << long_options[option_index].name << endl << flush;
			if(optarg) {
//				cout << " optarg: " << optarg << endl << flush;
			}
			break;

		case 1000:
			sscanf(optarg,"%d",&level_);
			break;
		case 1001:
			sscanf(optarg,"%d",&saveLevel_);
			break;
		case 1002:
			htmNames.push_back(optarg);
			break;

		case 1010:
			csvNames.push_back(optarg);
			break;
		case 1011:
			nameSize = csvNames.size();
			if(csvLevels.size() != (nameSize-1)) {
				while(csvLevels.size() < (nameSize-1)) {
					csvLevels.push_back("default");
				}
			}
			csvLevels.push_back(optarg);
			break;
		case 1012:
			nameSize = csvNames.size();
			if(csvLevels.size() != (nameSize-1)) {
				while(csvData.size() < (nameSize-1)) {
					csvData.push_back("default");
				}
			}
			csvData.push_back(optarg);
			break;
		case 1013:
			nameSize = csvNames.size();
			if(csvGridColor.size() != (nameSize-1)) {
				while(csvGridColor.size() < (nameSize-1)) {
					csvGridColor.push_back("default");
				}
			}
			csvGridColor.push_back(optarg);
			break;

		case 1020:
			cellsFiles.push_back(optarg);
			break;
		case 1021:
			nameSize = cellsFiles.size();
			if(cellsDataRange.size() < (nameSize-1)) {
				while(cellsDataRange.size() < (nameSize-1)) {
					cellsDataRange.push_back("default");
				}
			}
			cellsDataRange.push_back(optarg);
			break;
		case 1022:
			nameSize = cellsFiles.size();
			if(cellsGridColor.size() < (nameSize-1)) {
				while(cellsGridColor.size() < (nameSize-1)) {
					cellsGridColor.push_back("default");
				}
			}
			cellsGridColor.push_back(optarg);
			break;
		case 1023:
			nameSize = cellsFiles.size();
			if(cellsScale.size() < (nameSize-1)) {
				while(cellsScale.size() < (nameSize-1)) {
					cellsScale.push_back("default");
				}
			}
			cellsScale.push_back(optarg);
			break;
		case 1024:
			nameSize = cellsFiles.size();
			if(cellsCentroids.size() < (nameSize-1)) {
				while(cellsCentroids.size() < (nameSize-1)) {
					cellsCentroids.push_back("default");
				}
			}
			cellsCentroids.push_back(optarg);
			break;

		case 1100:
			cout << "testShapes: " << optarg << endl << flush;
			testShapeFiles_flag = true;
			if(strcmp(optarg,"default")) {
				shapeFile = optarg;
			}
			cout << "testShapes sf: " << shapeFile << endl << flush;
			break;
		case 1101:
			shapeFilesBlockingSphereScale_flag = true;
			sscanf(optarg,"%f",&shapeFilesBlockingSphereScale);
			break;

		case 2000:
			sscanf(optarg,"%f",&lineWidth);
			break;

		case 3000:
			lookFrom_flag = 1;
//			cout << "3000: optarg: " << optarg << endl;
			ss << optarg;
			ss >> lookFromArgs;
			break;

		case 't':
			// something
			cout << "--test: " << optarg << endl << flush;
			break;

		default:
			c = 0;
			abort ();
		}
	}

	if(!quiet_flag){
		cout << mainName << " starting." << endl;
	}
	if(help_flag) {
		cout << "help " << endl << flush;
		int iOpt = 0;
		while(long_options[iOpt].name != 0) {
			cout << " --" << long_options[iOpt].name << flush;
			cout << endl << flush;
			iOpt++;
		}
		return 0;
	}

	if(level_==0) level_ = 4;
	if(saveLevel_==0) saveLevel_ = 5;

	float64 PI = atan2(0,-1); // cout << "PI=" << PI << endl << flush;
	SpatialVector xHat = SpatialVector(1.,0.,0.);
	SpatialVector yHat = SpatialVector(0.,1.,0.);
	SpatialVector zHat = SpatialVector(0.,0.,1.);

	cout << "viz..." << flush;
	VizHTM *viz = new VizHTM(NARRAY_);
	cout << "allocated." << endl << flush;

	examiner_viz = true;

	if(csvNames.size()   > 0) { plotCsv(viz);	}
	if(cellsFiles.size() > 0) { plotCells(viz); }

//	testTenDegreeGrid_flag = true;

	if(blockingSphere_flag) plotBlockingSphere(viz,0.2,0.2,0.2,0.999);

	if(testTenDegreeGrid_flag) testTenDegreeGrid(viz);
	if(false) testTenDegreeGridRGB(viz,0.6,0.6,0.6);
	if(false) testHTMRange(viz,1,"N0","N0");
	if(false) testHTMRange(viz,2,"N033","N033");

	if(false) testHTMRange(viz,1,"N01","N01"); // One green triangle.
	if(false) testHTMRange(viz,2,"N01","N01"); // Subdivide into 4 child cells.

	if(false) testAddEdgesFromIndexAndName(viz,"N01"); // TODO Write tests for latlon & idByLatLon etc.

	if(false) testPlotEdgesFromHTMNameInterval(viz); // Simple subdivision for Kuo. 2016-0317

	// testPlotDataSetIntersection_flag = false; // fix intersect null hrBoundary bug.
	if(testPlotDataSetIntersection_flag) testPlotDataSetIntersection(viz);
	// Claim to fame...
//	testPlotDataSetIntersectionRangeContains_flag = true;s
	if(testPlotDataSetIntersectionRangeContains_flag) testPlotDataSetIntersection0_PlotHtmRangeContains(viz,level_); // 7

//	testPlotDataSetIntersectionNativeIntersect_flag = true;
	if(testPlotDataSetIntersectionNativeIntersect_flag) testPlotDataSetIntersection0_PlotIntersectTwoRectanglesOutput(viz,level_);

	if(testMultiResolutionHtmRange_flag) testMultiResolutionHtmRange(viz);


	if(htmNames.size()>0) {
		for(std::vector<string>::iterator it = htmNames.begin();
				it != htmNames.end(); ++it ) {
			string name = *it;
			testHTMRange(viz,level_,name.c_str(),name.c_str());
		}
	}

	if(false) {
		int level = 2;
		testHTMRange(viz,level,"N01","N01");
		testHTMRange(viz,level,"N11","N11");
		testHTMRange(viz,level,"N21","N21");
		testHTMRange(viz,level,"N31","N31");
	}

	if(false) {
		int level = 5;
		testHTMRange(viz,level,"N01","N01");
		testHTMRange(viz,level+1,"N11","N11");
		testHTMRange(viz,level+2,"N21","N21");
		testHTMRange(viz,level+3,"N31","N31");
	}

	if(testHTMRangeToLevel15_flag){
		char *htmName = new char[32];
		strcpy(htmName,"N0");
		for(int j=1; j<15; j++){
			testHTMRange(viz,j,htmName,htmName);
			strcat(htmName,"3");
		}
	}

	if(testLatLonSpiral_flag) testLatLonSpiral(viz);

	// examiner_viz = true;
	if(false) testAddRectangle0(viz);

	if(testAddRectangle_flag) testAddRectangle(viz,level_);

	if(testTwoRectangle_flag) testTwoRectangle(viz,level_);

	// broken
	if(false) { // Test Georgia multiple times
		int s0 = 49; int s1 = 49; int sDelta = 1;
		for(int s=s0; s<=s1; s+=sDelta){
			cout << "hullStep=" << s << endl << flush;
			float r = (s-s0)/(1.0*(s1-s0));
			float g = (s-s0)/(1.0*(s1-s0));
			float b = (s1-s)/(1.0*(s1-s0));
			r=1; g=1; b=1;
			testGeorgia(viz,8,s,r,g,b);
		}
	}
	if(testGeorgia_flag) { // Test Georgia
			float r=0.5; float g=0.8; float b=0.4;
			testGeorgia(viz,level_,-1,r,g,b);
	}

	if(testTriaxis_flag) testTriaxis(viz);
	if(false) testIJKRGBFace(viz);
	if(false) testAnEdge(viz);
	if(false) testText1(
			viz,
			unitVector(xHat+yHat+zHat));
	if(false) testShowEdgeProjections(
			viz,
			unitVector(xHat+yHat+zHat),
			0.99);

	// examiner_viz = true;
	if(false) testTwoConstraints(viz,5);

	if(testChunks_flag) testChunks(viz);

	if(testLevelChunk_flag) testLevelChunk(viz);

	if(testShapeFiles_flag) testShapeFiles(viz);

	if(false) testTIFF(viz);

	if(false) testDelaunay(viz);

	if(false) testNearestNeighbors(viz);


//	cout << "a8000" << endl << flush;

	// examiner_viz = true;
	// SwathNMQFlag
	if(testSwathNMQ_flag) { // Visualization using real data.
		initialize_nmq_trmm_viz();
		if(true) testSwath(viz);
		if(true) testnmq(viz);
		if(addFiduciaryTriangles) addFiduciaries(viz);
	}

//	cout << "a9000" << endl << flush;

	if(false) {
		swathGeometryHTM_viz = true;
		examiner_viz = true;
		testSwath(viz);
	}

//	cout << "a10000" << endl << flush;

	// current test.
	if(false) {
		if(false) testHstmRange(viz);
		if(true)  testMultiResolution(viz);
	}

	if(false) testTransparency(viz);

	if(false) viz->debug_dump();

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
			loadScene(offscreen->root,selectionRoot,offscreen->vpRegion);
		} else {
			setLookFrom(offscreen->root,selectionRoot,offscreen->vpRegion);
		}
		offscreen->saveImage(1);
	}

	cout << 10200 << endl << flush;

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
				SpatialVector *p = nmq_trmm_vizCenter;
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


