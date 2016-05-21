/*
 * VizHTM_main.C
 *
 *  Created on: Dec 16, 2015
 *      Author: mrilee
 */

#include "VizHTM.h"

#include <shapefil.h>
#include <geompack.h>

#include <unistd.h>
#include <getopt.h>

#include <QtOpenGL/QGL>
#include <QtGui/QImage>

#include <iostream>
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
#include <Inventor/sensors/SoTimerSensor.h>
#include <Inventor/SbBasic.h>

#include "VizHTM_main.h"

using namespace std;

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

	HtmRange *range = new HtmRange();

	range->purge();
	bool varlen_individualHTMIds = false; // true for individuals, false for ranges
	bool overlap = domain.intersect(index,range,varlen_individualHTMIds);

	Key lo = 0; Key hi = 0;
	SpatialVector v1,v2,v3;
//	if(false)
//		cout
//		<< " overlap: " << overlap
//		<< " nRanges: " << range->nranges() << flush;
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
	SpatialVector o = SpatialVector(0.0,0.0,0.0);

	//	SpatialVector a = SpatialVector(1.0,0.0,0.0);
	SpatialVector a = randomVector();
	//	SpatialVector a = SpatialVector(1.0,1.0,0.0); a.normalize();

	SpatialVector as = 1.05*a;

	//	float64 d = 0.9999;
	float64 d = 0.99;
	SpatialConstraint sc = SpatialConstraint(a,d);

	if(false) testConstraintAD(viz,a,d);

	RangeConvex convex = RangeConvex();
	convex.add(sc);

	//	SpatialVector a1 = a+randomVector()*0.15; a1.normalize();
	//	SpatialConstraint sc1 = SpatialConstraint(a1,0.99);

	SpatialVector a1 = -1*a; a1.normalize();
	SpatialConstraint sc1 = SpatialConstraint(a1,0.995);
	sc1.invert();
	// SpatialConstraint(randomVector(),0.99);
	viz->addConstraint(sc1,0.5,1.0,1.0);
	convex.add(sc1);

	testConvexHtmRangeIntersection(viz,convex,htmIdLevel);
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
//	cout << "1" << flush;
	if(true){
		RangeConvex rc = RangeConvex(u0,u1,u2,u3);
		domain1.add(rc);
	}
//	cout << "2" << flush;
	if(true){
		RangeConvex rc = RangeConvex(v0,v1,v2,v3);
		domain2.add(rc);
	}
//		return;
//	cout << "3" << flush;
	bool varlen_individualHTMIds = false; // true for individuals, false for ranges

//	cout << "." << flush;
//	HtmRange range1 = HtmRange();
	HtmRange range1;
	range1.purge();
	bool overlap1 = domain1.intersect(index,&range1,varlen_individualHTMIds);
//	cout << "." << flush;
	range1.defrag();
	range1.reset();
	//	if(range1.nranges()==0)return;
	rangeU->addRange(&range1);
	range1.reset();

//	cout << "." << flush;
//	cout << "4" << flush;
//	cout << " overlap1: " << overlap1 << ";" << flush;
//	cout << endl << flush;

	Key lo = 0, hi = 0;
	uint64 indexp = 0;

//	cout << " range1.ranges(): " << range1.nranges() << endl << flush;
	range1.reset();
	indexp = range1.getNext(lo,hi);
//	cout << " range1indexp " << indexp << endl << flush;
//	cout << " range1.lo,hi " << lo << " " << hi << endl << flush;
//	cout << "        level " << levelOfId(lo) << endl << flush;

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


//	cout << "5" << flush;

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
//	cout << "6" << flush;
//	HtmRange *resultRange = HTMRangeAtLevelFromIntersection(htmIdLevel,&range1,&range2);
	HtmRange *resultRange = range1.HTMRangeAtLevelFromIntersection(&range2,htmIdLevel);
	if(!resultRange) return;
//	cout << "7" << flush;
	rangeIntersection->addRange(resultRange); // Note: resultRange is copied piece by piece here.
//	cout << "8" << flush;
	rangeIntersection->defrag();
//	cout << "9" << flush;
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
	viz->addArc(*x0,target,1,1,1,-1.0,200);
	viz->addArc(target,*x1,1,1,1,-1.0,200);

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
	}

	const SpatialVector zHat = SpatialVector(0.,0.,1.);
	for(int i=0; i<10; i++){
//		cout << "i: " << i << endl << flush;

		double delta = 10.;
		double lat = -80. + 160.0*uniformDouble();
		double lon = 360.0 * uniformDouble();
//		cout << " udll: " << delta << " " << lat << " " << lon << endl << flush;
		u0.setLatLonDegrees(lat,lon);
		u1.setLatLonDegrees(lat+delta,lon);
		u2.setLatLonDegrees(lat+delta,lon+delta);
		u3.setLatLonDegrees(lat,lon+delta);

		double dlat = 12.0*(uniformDouble()-0.5);
		double dlon = 12.0*(uniformDouble()-0.5);
//		cout << "dlat,dlon: " << dlat << " " << dlon << endl << flush;
		double theta = 90.0*uniformDouble();
//		lat = uniformDouble(-80.,80.);
//		lon = uniformDouble(10.,350.);
		lat+=dlat; lon+=dlon;
//		cout << " vdll: " << delta << " " << lat << " " << lon << endl << flush;
		v0.setLatLonDegrees(lat,lon+sin(theta)*delta);
		v1.setLatLonDegrees(lat+cos(theta)*delta,lon+sin(theta)*delta);
		v2.setLatLonDegrees(lat+sin(theta)*delta,lon+cos(theta)*delta);
		v3.setLatLonDegrees(lat+sin(theta)*delta,lon);

		if(i==i){
//		if(i==2) {
//		if(false){
//		if(i==14){
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
				viz->addRectangle(u0,u1,u2,u3,r,g,b);
			}
			if(true){
				float r = 1.0;
				float g = 0.0;
				float b = 0.0;
				viz->addRectangle(v0,v1,v2,v3,r,g,b);
			}
		}
//		cout << "+" << endl << flush;
	}

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

	granuleSets = DataIntersectionDriver(viz);
	granules0 = &(granuleSets[0]);
	granules1 = &(granuleSets[1]);

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

	bool focus = false;
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
//		cout << " latlon count: " << count << endl << flush;
		int buildLevel = 5;
//		cout << 1000 << endl << flush;
		htmInterface *htm = new htmInterface(level,buildLevel);
//		cout << 1100 << endl << flush;
//		HTMRangeValueVector htmRangeVector = htm->convexHull(latlon1,hullSteps);
		HTMRangeValueVector htmRangeVector = htm->convexHull(latlon,hullSteps);

		SpatialVector *xStart = &(htm->polyCorners_[0].c_);
		SpatialVector *x0 = xStart;
//		cout << " polyCorners: 0 " << flush;
		for(int j=1; j<htm->polyCorners_.size();j++){
//			cout << j << " " << flush;
			SpatialVector *x1 = &(htm->polyCorners_[j].c_);
//			cout << *x0 << ", " << *x1 << "; " << flush;
			viz->addEdge((*x0)*1.001,(*x1)*1.001,hullR,hullG,hullB);
			x0=x1;
		}
//		cout << endl << flush;
		viz->addEdge((*x0)*1.001,(*xStart)*1.001,hullR,hullG,hullB);

//		cout << 1200 << endl << flush;
//		cout << " htmRangeVector-size: " << htmRangeVector.size() << endl << flush;
//		cout << " htmRangeVector-j:    ";
		for( int j=0; j<htmRangeVector.size(); j++) {
			htmRange hr = htmRangeVector[j];
//			cout << " (j=" << j << " " << hr.lo << " " << hr.hi << " ) ";
			plotHTMInterval(viz,htm->index(),hr);
		}
//		cout << endl << flush;
	} else {
		cout << "Couldn't open '" << fileName << "'." << endl << flush;
	}
}

void testAddEdgesFromIndexAndName(VizHTM* viz, const char* htmIdName, int saveLevel=5) {
	SpatialIndex *index = new SpatialIndex(htmIdName,saveLevel);
	viz->addEdgesFromIndexAndName(index,htmIdName,0.6,0.6,0.9);
}

void testShapeFiles(VizHTM* viz) {

	bool verbose = false;

	if(true) {

	plotBlockingSphere(viz,0.2,0.2,0.2,0.999);
	testTenDegreeGridRGB(viz,0.6,0.6,0.6);

	//	string shapeFile = "data/ne_110m_coastline/ne_110m_coastline.shp"; // okay
	string shapeFile = "data/ne_50m_coastline/ne_50m_coastline.shp"; // okay
	// string shapeFile = "data/ne_10m_coastline/ne_10m_coastline.shp"; // Doesn't work yet.
	float r = 0.5;
	float g = 0.9;
	float b = 0.9;

	SHPHandle hSHP = SHPOpen(shapeFile.c_str(),"rb");

	if(hSHP>0) {
		cout << "file: " << shapeFile << " found." << endl << flush;

		int nShapeType, nEntities;
		double adfMinBound[4], adfMaxBound[4];
		SHPGetInfo(hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );

		cout << "nEntities: " << nEntities << ", nShapeType: " << nShapeType << endl << flush;

		/* Debugging using string shapeFile = "data/ne_50m_coastline/ne_50m_coastline.shp";
		   Asia and Africa
			int nStart = 1387;
			int nEnd   = 1388;
		*/

		int nStart = 0;
		int nEnd   = nEntities;

		for(int i=nStart; i < nEnd; i++ ) {
			if(verbose){
				if(i % 100 == 0) {
					cout << "on entity " << i << endl << flush;
				}
			}
			SHPObject* psShape = SHPReadObject(hSHP,i);
			if(verbose){
				cout << "<psShape>" << endl
						<< " nSHPType  " << psShape->nSHPType  << " " << endl
						<< " nVertices " << psShape->nVertices << " " << endl
						<< " nParts    " << psShape->nParts    << " " << endl;
				for(int ipart=0; ipart<psShape->nParts; ipart++) {
					cout << "  part [" << ipart << "] = " << psShape->panPartStart[ipart] << endl;
				}
				cout << "</psShape>" << endl;
			}
			if(psShape->nSHPType == SHPT_ARC) {
				double *yA = psShape->padfY;
				double *xA = psShape->padfX;
				int nVerts = psShape->nVertices;
				int jStart = 0;
				if(psShape->nParts != 0) {
					for(int j=0; j < psShape->nParts; j++) {
						jStart = psShape->panPartStart[j];
						yA = psShape->padfY + psShape->panPartStart[j];
						xA = psShape->padfX + psShape->panPartStart[j];
						if(j < psShape->nParts-1) {
							nVerts = psShape->panPartStart[j+1] - psShape->panPartStart[j];
						} else {
							nVerts = psShape->nVertices - psShape->panPartStart[j];
						}
						if(verbose) cout << "<addArcParts j,start=[" << j << ", " << jStart << "] nVerts= " << nVerts << " />" << endl;
						if( true ) {
//							float r_ = j/(psShape->nParts-1.0);
							float r_ = r;
							viz->addArcsFromLatLonDegrees(
									yA, xA, nVerts,
									false,
									r_,g,b,-1.,3
							);
						}
					}
				} else {
					if(verbose) cout << "<addArc 0,start=[ 0, " << jStart << "] nVerts= " << nVerts << " />" << endl;
					if( true ) {
						viz->addArcsFromLatLonDegrees(
								yA, xA, nVerts,
								false,
								0.,g,b,-1.,3
						);
					}
				}
			}
		}
		if(verbose) cout << endl;
		SHPClose(hSHP);
	}

	}

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
		testHTMRangeToLevel15_flag=0;

	static int
		level_     = 0,
		saveLevel_ = 0;

	vector<string> htmNames; htmNames.clear();

	char *cvalue = NULL;
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

			{"test0",     required_argument,  0, 0},
			{"test",      required_argument,  0, 't'},

			{0, 0, 0, 0}
	};
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

	if(testTenDegreeGrid_flag) testTenDegreeGrid(viz);
	if(false) testTenDegreeGridRGB(viz,0.6,0.6,0.6);
	if(false) testHTMRange(viz,1,"N0","N0");
	if(false) testHTMRange(viz,2,"N033","N033");

	if(false) testHTMRange(viz,1,"N01","N01"); // One green triangle.
	if(false) testHTMRange(viz,2,"N01","N01"); // Subdivide into 4 child cells.

	if(false) testAddEdgesFromIndexAndName(viz,"N01"); // TODO Write tests for latlon & idByLatLon etc.

	if(false) testPlotEdgesFromHTMNameInterval(viz); // Simple subdivision for Kuo. 2016-0317

	if(testPlotDataSetIntersection_flag) testPlotDataSetIntersection(viz);
	if(testPlotDataSetIntersectionRangeContains_flag) testPlotDataSetIntersection0_PlotHtmRangeContains(viz,level_); // 7
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
	if(testAddRectangle_flag) testAddRectangle(viz,level_);
	if(testTwoRectangle_flag) testTwoRectangle(viz,level_);

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
	if(false) testTwoConstraints(viz,5);

	if(testLevelChunk_flag) testLevelChunk(viz);

	if(false) testShapeFiles(viz);

	if(true) testDelaunay(viz);

	if(true) testNearestNeighbors(viz);

	if(false) viz->debug_dump();

	/**************** Start Graphics ****************/
	QWidget *window = SoQt::init(argv[0]);
	if (window == NULL) exit(1);

	SoSelection *selectionRoot = new SoSelection;
	selectionRoot->policy = SoSelection::SINGLE;

	SoSeparator *root = new SoSeparator;
	//	root->ref(); // TODO Figure out ->ref();

	root->addChild(viz->makeRoot());

	selectionRoot->addChild(root);

	SoQtExaminerViewer *viewer = new SoQtExaminerViewer(window);
	viewer->setSceneGraph(selectionRoot);
	viewer->setTitle(mainName);
	viewer->show();

	SoQt::show(window);
	SoQt::mainLoop();

	delete viewer;
	//	root->unref();

	cout << mainName << " done." << endl;

	return 0;
}


