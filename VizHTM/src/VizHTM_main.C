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

#include "SpatialException.h"
#include "SpatialIndex.h"
#include "SpatialVector.h"
#include "SpatialInterface.h"

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
	cout << " a in sc(a,d)? " << sc.contains(a) << endl << flush;
	cout << " acos(a*a),s " << acos(a*a) << " " << sc.coneAngle() << endl << flush;
	cout << " i in sc(a,d)? " << sc.contains(SpatialVector(1.0,0.0,0.0)) << endl << flush;
	cout << " acos(i*a),s " << acos(SpatialVector(1.0,0.0,0.0)*a) << " " << sc.coneAngle() << endl << flush;
	testConstraint(viz,sc);
}

void testConvexHtmRangeIntersection(VizHTM * viz, RangeConvex convex, int htmIdDepth) {

	int saveDepth = 5;
	SpatialIndex *index = new SpatialIndex(htmIdDepth,saveDepth);
	SpatialDomain domain = SpatialDomain(index);
	domain.add(convex);
	domain.setOlevel(htmIdDepth);

	HtmRange *range = new HtmRange();

	range->purge();
	bool varlen_individualHTMIds = false; // true for individuals, false for ranges
	bool overlap = domain.intersect(index,range,varlen_individualHTMIds);

	Key lo = 0; Key hi = 0;
	SpatialVector v1,v2,v3;
	if(false)
		cout
		<< " overlap: " << overlap
		<< " nRanges: " << range->nranges() << flush;
	range->defrag();
	if(false)
		cout
		<< " nRangeDefrag: " << range->nranges()
		<< " nConvexes: " << domain.numConvexes()
		<< endl << flush;

	//		overlap = false;

	if(overlap) {
		range->reset();
		uint64 indexp = range->getNext(lo,hi);
		//		indexp = range->getNext(lo,hi);
		// cout << " indexp,lo,hi: " << indexp << " " << lo << " " << hi << endl << flush;
		int k = 0; int kMax = 10000;
		int triangles = 0; int tMax = 10000;
		// cout << " index-depth: " << index->getMaxlevel() << " -saveDepth: " << index->getBildLevel() << endl << flush;
		do {
			k++;
			//			for(uint64 nodeIndex=lo; nodeIndex <= min(lo+10,hi); nodeIndex++){
			for(uint64 numericId=lo; numericId <= hi && triangles < tMax; numericId++){
				triangles ++;
				uint64 nodeIndex = index->nodeIndexFromId(numericId);

				if(false){
					cout
					<< "working " << flush
					<< " triangles: " << triangles << flush
					<< " numericId: " << numericId << flush
					<< " nodeIndex: " << nodeIndex << flush;
					cout << " k=" << k << flush;
				}

				index->nodeVertex(nodeIndex,v1,v2,v3);

				if(false){
					cout << "..." << flush;
					cout
					<< " v1: " << v1.x() << " " << v1.y() << " " << v1.z()
					<< " v2: " << v2.x() << " " << v2.y() << " " << v2.z()
					<< " v3: " << v1.x() << " " << v3.y() << " " << v3.z() << flush;
				}

				if(true){
					float size = pow(0.5,htmIdDepth+3);
					float r = 0.4, g = 0.4, b = 0.6;
					SpatialVector x = 3.*v1+v2+v3; x.normalize(); x *= 1.0+1.0e-6;
					SpatialVector x_ = v1+v2+v3; x_.normalize();
					if(false){
						cout
						<< " nI: " << nodeIndex
						<< " x: " << x.x() << " " << x.y() << " " << x.z();
					}
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
				if(false) cout << endl << flush;;
			}
		} while (range->getNext(lo,hi) && k < kMax && triangles < tMax);
	}
	//	cout << "overlap done " << endl << flush;



}

void testTwoConstraints(VizHTM *viz, int htmIdDepth) {
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

	testConvexHtmRangeIntersection(viz,convex,htmIdDepth);
}

void testAddRectangle(VizHTM *viz, int htmIdDepth) {
	int saveDepth = 5;

	SpatialIndex *index = new SpatialIndex(htmIdDepth,saveDepth);
	SpatialDomain domain = SpatialDomain(index);
	domain.setOlevel(htmIdDepth); // Note this sets the olevel on the convexes.

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
	rc->setOlevel(htmIdDepth); // Note this is supposed to be done when added to the domain.
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
	// The indexes lo and hi are not restricted to being at one depth.
	//	cout << "indexp: " << indexp << endl << flush;
	SpatialVector x1,x2,x3;
	do {
//		cout << "original lo,hi= " << lo << " " << hi << " " << flush;

		// Leaves are at a depth of maxlevel+1

		int level = htmIdDepth + 1; // htmIdDepth is used to set maxlevel in the index. aka olevel.

		int depthLo = depthOfId(lo);
		if(depthLo<level) {
			lo = lo << (2*(level-depthLo));
		}
		int depthHi = depthOfId(hi);
		if(depthHi<level) {
			for(int shift=0; shift < (level-depthHi); shift++) {
				hi = hi << 2;
				hi += 3; // TODO Determine if I am overthinking it?
			}
		}
//		cout << "   fixed lo,hi= " << lo << " " << hi << " "
//				<< "depth lo,hi,level= "
//				<< depthLo << " "
//				<< depthHi << " "
//				<< level << " "
//				<< endl << flush;
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
						float size = pow(0.5,htmIdDepth+3);
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

struct KeyPair {
	Key lo; Key hi;
};
KeyPair HTMRangeAtDepthFromHTMRange(int htmIdDepth, Key lo, Key hi) {
	int level = htmIdDepth + 1; // htmIdDepth is used to set maxlevel in the index. aka olevel.
	int depthLo = depthOfId(lo);
	if(depthLo<level) {
		lo = lo << (2*(level-depthLo));
	}
	int depthHi = depthOfId(hi);
	if(depthHi<level) {
		for(int shift=0; shift < (level-depthHi); shift++) {
			hi = hi << 2;
			hi += 3; // TODO Determine if I am overthinking it?
		}
	}
	KeyPair depthAdaptedRange;
	depthAdaptedRange.lo = lo;
	depthAdaptedRange.hi = hi;
	return depthAdaptedRange;
}
HtmRange *HTMRangeAtDepthFromIntersection(int htmIdDepth,HtmRange *range1, HtmRange *range2){
//	cout << "Comparing..." << endl << flush;
	HtmRange *resultRange = new HtmRange();
	resultRange->purge();
	Key lo1,hi1,lo2,hi2;
	range1->reset();
	uint64 indexp1 = range1->getNext(lo1,hi1);
	do {
		KeyPair testRange1 = HTMRangeAtDepthFromHTMRange(htmIdDepth,lo1,hi1);
		range2->reset();
		uint64 indexp2 = range2->getNext(lo2,hi2);
		bool intersects = false;
		do {
			KeyPair testRange2 = HTMRangeAtDepthFromHTMRange(htmIdDepth,lo2,hi2);
			intersects = testRange2.lo <= testRange1.hi
					&& testRange2.hi >= testRange1.lo;
//			cout << "lh1,lh2: "
//					<< lo1 << " " << hi1 << ", "
//					<< lo2 << " " << hi2 << ", "
//					<< intersects << flush;
			if(intersects){
				Key lo_ = max(testRange1.lo,testRange2.lo);
				Key hi_ = min(testRange1.hi,testRange2.hi);
				resultRange->addRange(lo_,hi_);
//				cout << ", added "
//						<< lo_ << " "
//						<< hi_ << flush;
			}
//			cout << "." << endl << flush;
		} while (range2->getNext(lo2,hi2));
	} while (range1->getNext(lo1,hi1));
	resultRange->defrag();
	return resultRange;
}

void testTwoRectangle(VizHTM *viz, int htmIdDepth) {
	int saveDepth = 5;

	SpatialIndex *index = new SpatialIndex(htmIdDepth,saveDepth);

	SpatialDomain domain1 = SpatialDomain(index);
	domain1.setOlevel(htmIdDepth); // Note this sets the olevel on the convexes.

	SpatialDomain domain2 = SpatialDomain(index);
	domain2.setOlevel(htmIdDepth); // Note this sets the olevel on the convexes.

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
	rc->setOlevel(htmIdDepth); // Note this is supposed to be done when added to the domain.
	domain1.add(*rc);

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
	rc->setOlevel(htmIdDepth); // Note this is supposed to be done when added to the domain.
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

	HtmRange *resultRange = HTMRangeAtDepthFromIntersection(htmIdDepth,range1,range2);

	HtmRange *range = resultRange;
	range->reset();
	indexp = range->getNext(lo,hi);
	SpatialVector x1,x2,x3;
//	if(indexp)
	do {
//		cout << "original lo,hi= " << lo << " " << hi << " " << flush;

		// Leaves are at a depth of maxlevel+1

		int level = htmIdDepth + 1; // htmIdDepth is used to set maxlevel in the index. aka olevel.

		int depthLo = depthOfId(lo);
		if(depthLo<level) {
			lo = lo << (2*(level-depthLo));
		}
		int depthHi = depthOfId(hi);
		if(depthHi<level) {
			for(int shift=0; shift < (level-depthHi); shift++) {
				hi = hi << 2;
				hi += 3; // TODO Determine if I am overthinking it?
			}
		}
//		cout << "   fixed lo,hi= " << lo << " " << hi << " "
//				<< "depth lo,hi,level= "
//				<< depthLo << " "
//				<< depthHi << " "
//				<< level << " "
//				<< endl << flush;
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
						float size = pow(0.5,htmIdDepth+3);
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


void testTenDegreeGrid(VizHTM *viz, int htmIdDepth) {
	int saveDepth = 5;
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
	float *last_xyz= xyzFromLatLonDegrees(last_lat,last_lon);

	for(int i=1; i <= nEdges; i++) {
		float lat = lat_(i,nPoints);
		float lon = lon_(i,nPoints);
		float r   = r_  (i,nPoints);
		float g   = 1.0;
		float b   = 1.0;
		float *xyz= xyzFromLatLonDegrees(lat,lon);

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

	if(false) testLatLonSpiral(viz);
	if(false) testAddRectangle(viz,4);
	if(true) testTwoRectangle(viz,6);
	if(true) testTenDegreeGrid(viz,5);

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


