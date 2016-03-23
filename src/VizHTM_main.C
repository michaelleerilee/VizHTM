/*
 * VizHTM_main.C
 *
 *  Created on: Dec 16, 2015
 *      Author: mrilee
 */




#include "VizHTM.h"

#include <iostream>
#include <iomanip>
#include <random>
#include <fstream>
#include <vector>

#include "SpatialException.h"
#include "SpatialIndex.h"
#include "SpatialVector.h"
#include "SpatialInterface.h"

#include "BitShiftNameEncoding.h"

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
 * struct KeyPair is a lightweight aggregate a single interval [lo,hi] for HtmRange
 */
struct KeyPair {
	Key lo; Key hi;
};

/**
 * Translate an HtmRange to one at a greater level.  If the desired level is
 * less that the level implicit in the input range (lo & hi), then just return
 * an HtmRange constructed from the input range without modification.
 * @param htmIdLevel
 * @param lo the low end of the range
 * @param hi the high end of the range
 * @return levelAdaptedRange an HtmRange
 */
KeyPair HTMRangeAtLevelFromHTMRange(int htmIdLevel, Key lo, Key hi) {
	// htmIdLevel is used to set maxlevel in the index. aka olevel.
	int levelLo = levelOfId(lo);
	if(levelLo<htmIdLevel) {
		lo = lo << (2*(htmIdLevel-levelLo));
	}
	int levelHi = levelOfId(hi);
	if(levelHi<htmIdLevel) {
		for(int shift=0; shift < (htmIdLevel-levelHi); shift++) {
			hi = hi << 2;
			hi += 3; /// Increment hi by 3 to catch all of the faces at that level.
		}
	}
	KeyPair levelAdaptedRange;
	levelAdaptedRange.lo = lo;
	levelAdaptedRange.hi = hi;
	return levelAdaptedRange;
}

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
		VizHTM *viz,
		const SpatialIndex  *index,
		const SpatialVector *u0,
		const SpatialVector *u1,
		const SpatialVector *u2,
		const SpatialVector *u3,
		const SpatialVector *v0,
		const SpatialVector *v1,
		const SpatialVector *v2,
		const SpatialVector *v3
		) {
	int htmIdLevel = index->getLeafLevel();

	SpatialDomain domain1 = SpatialDomain(index);
//	domain1.setOlevel(htmIdLevel); // Note this sets the olevel on the convexes.

	SpatialDomain domain2 = SpatialDomain(index);
//	domain2.setOlevel(htmIdLevel); // Note this sets the olevel on the convexes.

//	cout << "1" << flush;
	if(true){
//		SpatialVector *v0 = VectorFromLatLonDegrees(30.0,20.0);
//		SpatialVector *v1 = VectorFromLatLonDegrees(30.0,-10.0);
//		SpatialVector *v2 = VectorFromLatLonDegrees(20.0,0.0);
//		SpatialVector *v3 = VectorFromLatLonDegrees(20.0,30.0);
		float r = 0.0;
		float g = 1.0;
		float b = 1.0;
//		viz->addRectangle(*u0,*u1,*u2,*u3,r,g,b);

		RangeConvex rc = RangeConvex(u0,u1,u2,u3);
		//	cout << "nConstraints: " << rc->numConstraints() << endl << flush;
		//	SpatialConstraint *sc = new SpatialConstraint(SpatialVector(0.,0.,1.),0.5);
		//	viz->addConstraint(*sc,1.0,1.0,1.0);
		//	rc->add(*sc);
//		rc.setOlevel(htmIdLevel); // Note this is supposed to be done when added to the domain.
		domain1.add(rc);
	}
//	cout << "2" << flush;
	if(true){
//		SpatialVector *v0 = VectorFromLatLonDegrees(10.0,0.0);
//		SpatialVector *v1 = VectorFromLatLonDegrees(30.0,-10.0);
//		SpatialVector *v2 = VectorFromLatLonDegrees(60.0,10.0);
//		SpatialVector *v3 = VectorFromLatLonDegrees(40.0,20.0);
		float r = 1.0;
		float g = 1.0;
		float b = 0.0;
//		viz->addRectangle(*v0,*v1,*v2,*v3,r,g,b);

		RangeConvex rc = RangeConvex(v0,v1,v2,v3);
		//	cout << "nConstraints: " << rc->numConstraints() << endl << flush;
		//	SpatialConstraint *sc = new SpatialConstraint(SpatialVector(0.,0.,1.),0.5);
		//	viz->addConstraint(*sc,1.0,1.0,1.0);
		//	rc->add(*sc);
//		rc.setOlevel(htmIdLevel); // Note this is supposed to be done when added to the domain. // !!!BUG!!!
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
	if(range1.nranges()==0)return;

//	cout << "." << flush;
//
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

	HtmRange range2 = HtmRange();
	range2.purge();
	bool overlap2 = domain2.intersect(index,&range2,varlen_individualHTMIds);
//	cout << "." << flush;
	range2.defrag();
	range2.reset();
	if(range2.nranges()==0)return;

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
	HtmRange *resultRange = HTMRangeAtLevelFromIntersection(htmIdLevel,&range1,&range2);
	if(!resultRange) return;

//	cout << "5" << flush;
	HtmRange *range = resultRange;

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

	delete resultRange;
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
	cout << "loId,hiId: " << loId << " " << hiId << endl << flush;
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

SpatialVector *pointFromReferenceAndDeltas(SpatialVector target, SpatialVector aheadOfTarget, float64 deltaAlongDegrees, float64 deltaCrossDegrees) {
	float64 k = 2.*M_PI/360.0;
	SpatialVector L = target^aheadOfTarget;
	return new SpatialVector(
			sin(k*deltaCrossDegrees)*L
			+cos(k*deltaCrossDegrees)*
			(sin(k*deltaAlongDegrees)*aheadOfTarget
					+ cos(k*deltaAlongDegrees)*target));
}

vector<SpatialVector> *plotCircularTrack(
		VizHTM *viz,
		float64 latDegrees, float64 lonDegrees,
		float64 groundTrackDegreesFromEast,
		float r, float g, float b, float alpha,
		float64 thetaAlongTrackDegrees0, float64 thetaAlongTrackDegrees1,
		float64 phiCrossTrackDegrees0, float64 phiCrossTrackDegrees1,
		float64 deltaAlongTrackDegrees, float64 deltaCrossTrackDegrees
		) {

	vector<SpatialVector> *granules = new vector<SpatialVector>;

	SpatialVector target; target.setLatLonDegrees(latDegrees,lonDegrees);
	SpatialVector north = SpatialVector(0.,0.,1.);
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

	SpatialVector* x0 =
			pointFromReferenceAndDeltas(target,aheadOfTarget,thetaAlongTrackDegrees0,0.);
//			sin(k*thetaAlongTrackDegrees0)*aheadOfTarget
//			+ cos(k*thetaAlongTrackDegrees0)*target;
	SpatialVector* x1 =
			pointFromReferenceAndDeltas(target,aheadOfTarget,thetaAlongTrackDegrees1,0.);
//			sin(k*thetaAlongTrackDegrees1)*aheadOfTarget
//			+ cos(k*thetaAlongTrackDegrees1)*target;
	viz->addArc(*x0,*x1,1.,1.,1.);

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
//			viz->addArc(*v0,*v1,r,g,b,alpha);
//			viz->addArc(*v1,*v2,r,g,b,alpha);
//			viz->addArc(*v2,*v3,r,g,b,alpha);
//			viz->addArc(*v3,*v0,r,g,b,alpha);
			granules->push_back(*v0);
			granules->push_back(*v1);
			granules->push_back(*v2);
			granules->push_back(*v3);
		}
	}
	return granules;
}


void testPlotDataSetIntersection(VizHTM *viz) {
	SpatialIndex index(5,5);
	SpatialVector u0, u1, u2, u3, v0, v1, v2, v3;

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
				&v0,&v1,&v2,&v3
				);
	}

	const SpatialVector zHat = SpatialVector(0.,0.,1.);
	for(int i=0; i<600; i++){
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
					&v0,&v1,&v2,&v3
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
	}
}

void testPlotDataSetIntersection0(VizHTM *viz) {
	viz->lineWidth = 1.5;
	viz->sphereComplexity = 1.0;
	plotBlockingSphere(viz,0.3,0.3,0.3,0.99);

	vector<SpatialVector> *granules0, *granules1;

	{
	float64 latDegrees = 0.0;
	float64 lonDegrees = 45.0;
	float64 groundTrackDegreesFromEast = 45.0;
	float r = 0.3;
	float g = 0.3;
	float b = 1.0;
	float a = -1;
	float64 thetaAlongTrackDegrees0 = -50;
	float64 thetaAlongTrackDegrees1 =  50;
	float64 phiCrossTrackDegrees0 =   -10;
	float64 phiCrossTrackDegrees1 =    10;
	float64 deltaAlongTrackDegrees = 10;
	float64 deltaCrossTrackDegrees = 10;

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
	float64 latDegrees = 0.0;
	float64 lonDegrees = 45.0;
	float64 groundTrackDegreesFromEast = 5.0;
	float r = 0.0;
	float g = 0.8;
	float b = 0.8;
	float a = -1;
	float64 thetaAlongTrackDegrees0 = -40;
	float64 thetaAlongTrackDegrees1 =  40;
	float64 phiCrossTrackDegrees0 =   -5;
	float64 phiCrossTrackDegrees1 =    5;
	float64 deltaAlongTrackDegrees = 5;
	float64 deltaCrossTrackDegrees = 5;

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

	bool focus = false;
	int iFocus = 0;
	int jFocus = 0;
	SpatialIndex index(5,5);
//	cout << "intersecting" << endl << flush;
	int count=0;
	int i=0, j=0;
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
						&v0,&v1,&v2,&v3
				);
			} else {
				if((iFocus==i)&&(jFocus==j)){
					intersectTwoRectangles(
							viz,
							&index,
							&u0,&u1,&u2,&u3,
							&v0,&v1,&v2,&v3
					);
					viz->addRectangle(u0,u1,u2,u3,1.0,0.,0.);
					viz->addRectangle(v0,v1,v2,v3,0.0,0.,1.);
				}
			}
//			cout <<"+" << flush;

//			count++;
//			if(count > 2000) return;
//			return;
		}
	}

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

int main(int argc, char *argv[]) {
	const char* mainName = "VizHTM-main";
	cout << mainName << " starting." << endl;

	float64 PI = atan2(0,-1); // cout << "PI=" << PI << endl << flush;
	SpatialVector xHat = SpatialVector(1.,0.,0.);
	SpatialVector yHat = SpatialVector(0.,1.,0.);
	SpatialVector zHat = SpatialVector(0.,0.,1.);

	cout << "viz..." << flush;
	VizHTM *viz = new VizHTM(NARRAY_);
	cout << "allocated." << endl << flush;

	if(false) testTenDegreeGrid(viz);
	if(false) testTenDegreeGridRGB(viz,0.6,0.6,0.6);
	if(false) testHTMRange(viz,1,"N0","N0");
	if(false) testHTMRange(viz,2,"N033","N033");

	if(false) testHTMRange(viz,1,"N01","N01"); // One green triangle.
	if(false) testHTMRange(viz,2,"N01","N01"); // Subdivide into 4 child cells.

	if(false) testAddEdgesFromIndexAndName(viz,"N01"); // TODO Write tests for latlon & idByLatLon etc.

	if(false) testPlotEdgesFromHTMNameInterval(viz); // Simple subdivision for Kuo. 2016-0317

	if(true) testPlotDataSetIntersection(viz);

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

	if(false){
		char *htmName = new char[32];
		strcpy(htmName,"N0");
		for(int j=1; j<15; j++){
			testHTMRange(viz,j,htmName,htmName);
			strcat(htmName,"3");
		}
	}

	if(false) testLatLonSpiral(viz);
	if(false) testAddRectangle(viz,4);
	if(false) testTwoRectangle(viz,6);
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
	if(false) { // Test Georgia
			float r=0.5; float g=0.8; float b=0.4;
			testGeorgia(viz,6,-1,r,g,b);
	}

	if(false) testTriaxis(viz);
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


