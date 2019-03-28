/*
 * Chunking.C
 *
 *  Created on: Mar 24, 2019
 *      Author: mrilee
 *
 *  Copyright (C) 2019 Rilee Systems Technologies LLC
 */

#include "STARE.h"
#include "HstmRange.h"
#include "tests.h"

#include <mfhdf.h>
#include <hdf.h>
#include <HdfEosDef.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <array>

using namespace std;

void testShapeFiles(VizHTM *viz, float r, float g, float b, float deltaZ);

void AllTrixelsAtLevelAcrossBins(int level, int nbins, HstmRange ranges[]) {
	EmbeddedLevelNameEncoding lj;

	stringstream ss;
	ss << "S0"; for(int i=0; i<level; ++i) ss << "0";
	lj.setName(ss.str().c_str());
	uint64 ljId = lj.getId();

	// ss << "N3"; for(int i=0; i<level; ++i) ss << "3";
	// lj.setName(ss.str().c_str());
	// uint64 lj1 = lj.getId();


	// ljId = lj.increment(ljId,level);
	int nChunks = 8*int(pow(4,level));
	// for( int iChunk=0; iChunk < nChunks; ++iChunk ) {
	for( int iChunk=0; iChunk < nChunks; ++iChunk ) {
		if(false) {
			cout << iChunk << " " << (iChunk % nbins) << " " << hex << ljId << dec << endl << flush;
		}
		if(false){
			uint64 a_ = ljId, b_ = ljId;
			EmbeddedLevelNameEncoding leftJustifiedEncoding;
			Key a = leftJustifiedEncoding.maskOffLevelBit(a_);
			int aLevel = leftJustifiedEncoding.levelById(a_);
			Key b = leftJustifiedEncoding.maskOffLevelAndLevelBit(b_);
			int bLevel = leftJustifiedEncoding.levelById(b_);
			leftJustifiedEncoding.setId(b_);
			b = leftJustifiedEncoding.getIdTerminator_NoDepthBit();
			// cout << "a_,b_ " << hex << a_ << "," << b_ << dec << endl << flush;
			// cout << "a,b " << hex << a << "," << b << dec << endl << flush;
		}
		ranges[(iChunk % nbins)].addRange(ljId,ljId);
		if(false){
			cout << iChunk << " ranges loop nranges " << ranges[iChunk % nbins].range->nranges() << endl << flush;
		}
		ljId = lj.increment(ljId,level);
	}
	if(false){
		for( int ir = 0; ir < nbins; ++ir) {
			cout << ir << " ranges ret nranges " << ranges[ir].range->nranges() << endl << flush;
		}
	}
}

bool Chunk1( VizHTM *viz ) {
	bool ok = false;

	cout << "Hello world" << endl;

	// MLR Doesn't work Fedora+GLX viz->lineWidth = 1;

	testShapeFiles(viz,0,1,1,0.03);

	STARE index;
	HstmRange* range = new HstmRange;
	int resolutionLevel = 4;
	// int resolutionLevel = 5;
	// int resolutionLevel = 6;
	// int resolutionLevel = 7;
	// int resolutionLevel = 8;
	// int resolutionLevel = 14;
	SpatialIndex sIndex = index.getIndex(resolutionLevel);

	SpatialIndex sIndex0 = index.getIndex(0);
	HstmRange* range0 = new HstmRange;

	int ncolor = 16;
	// float color_scale = 0.5;
	// float color_scale = 0.75;
	float color_scale = 1.0;
	float r[ncolor],g[ncolor],b[ncolor];
	for( int ic = 0; ic < ncolor; ++ic ) {
		r[ic] = color_scale*ic/(ncolor-1.0);
		g[ic] = color_scale*ic/(ncolor-1.0);
		b[ic] = color_scale*ic/(ncolor-1.0);
		// b[ic] = color_scale*(ncolor-ic-1)/(ncolor-1.0);
	}
	HstmRange range_edges;
	HstmRange range_faces[ncolor];

	cout << "ncolor: " << ncolor << endl << flush;

	int chunkLevel = 6;
	cout << "chunkLevel, nChunks: " << chunkLevel << ", " << 8*int(pow(4,chunkLevel)) << endl << flush;
//	EmbeddedLevelNameEncoding lj;
//	lj.setName("S00000");
//	uint64 lj0 = lj.getId();
//	lj.setName("N33333");
//	uint64 lj1 = lj.getId();
	// cout << hex;
	// cout << "lj0,lj1: " << lj0 << "," << lj1 << endl << flush;
	// uint64 ljx = lj.decrement(lj1,4);
	// cout << "ljx: " << ljx << endl << flush;
	// ljx = lj.decrement(ljx,4);
	// cout << "ljx: " << ljx << endl << flush;
	// ljx = lj.increment(lj0,4);
	// cout << "ljx: " << ljx << endl << flush;
	// ljx = lj.increment(ljx,4);
	// cout << "ljx: " << ljx << endl << flush;
	// cout << dec;
	// exit(0);
	// 8*4**level
	// int nodeFromChunkNumber[8*int(pow(4,chunkLevel)))];

	if(true){
		HstmRange ranges[ncolor];
		AllTrixelsAtLevelAcrossBins(chunkLevel,ncolor,ranges);
		SpatialIndex sIndex = index.getIndex(chunkLevel);
		float re=0, ge=0.75, be=0;
		if(true) {
			for( int ic = 0; ic < ncolor; ++ic ) {
				// cout << ic << " ranges size " << ranges[ic].range->nranges() << endl << flush;
				viz->addHstmRangeFaces(&(ranges[ic]),r[ic],g[ic],b[ic],0.0,1.0,0.01,&sIndex);
				//viz->addHstmRange(&(ranges[ic]),re,ge,be,0.0,1.0,true,0.02,&sIndex);
			}
		}
		// int ic = 15;
		// re = 1; ge = 0; be = 0;
		// viz->addHstmRange(&(ranges[ic]),re,ge,be,0.0,1.0,true,0.02,&sIndex);
		// re = 0; ge = 1; be = 0;
		// viz->addHstmRange(&(ranges[++ic]),re,ge,be,0.0,1.0,true,0.02,&sIndex);
		// re = 0; ge = 0; be = 1;
		// viz->addHstmRange(&(ranges[++ic]),re,ge,be,0.0,1.0,true,0.02,&sIndex);

//		int ic;
//		ic = 5;	re = 1; ge = 1; be = 0;
//		viz->addHstmRange(&(ranges[ic]),re,ge,be,0.0,1.0,true,0.02,&sIndex);
//		viz->addHstmRangeFaces(&(ranges[ic]),r[ic],g[ic],b[ic],0.0,1.0,0.01,&sIndex);
//		ic = 6; re = 1; ge = 0; be = 1;
//
//		ic = 7; re = 0; ge = 1; be = 1;
//
//
//		ic = 12; re = 0.5; ge = 0.5; be = 1;
//
//		ic = 13; re = 0.5; ge = 1; be = 0.5;
//
//		ic = 14; re = 1; ge = 0.5; be = 0.5;

	}

	if(true){
		int nbins = 1, nlevels = 6;
		float re=1.0, ge=1.0, be=0.5, zlayer = 0.02, rgb_scale = 0.7;
		HstmRange ranges0[nbins];
		for(int level=0; level<nlevels; ++level) {
			for(int ib=0; ib < nbins; ++ib) {
				ranges0[ib].purge();
			}
			AllTrixelsAtLevelAcrossBins(level,nbins,ranges0);
			SpatialIndex sIndex0 = index.getIndex(level);
			for( int ic = 0; ic < nbins; ++ic ) {
				// cout << "ranges0 size " << ranges0[ic].range->nranges() << endl << flush;
				// be = (1.0*ic)/nbins;
				viz->addHstmRange(&(ranges0[ic]),re,ge,be,0.0,1.0,true,zlayer,&sIndex0);
				// re *= 0.75; ge *= 0.75; be *= 0.75; zlayer -= 0.001;
				re *= rgb_scale; ge *= rgb_scale; be *= rgb_scale; zlayer -= 0.001;
			}
		}
	}


//	int nChunks = 8*int(pow(4,chunkLevel));
//	uint64 ljIdFromChunkNumber[nChunks];
//
//	int iChunk = 0;
//	ljIdFromChunkNumber[iChunk] = lj0;
//	range_faces[0].addRange(ljIdFromChunkNumber[iChunk],ljIdFromChunkNumber[iChunk]);
//	for(iChunk=1; iChunk < nChunks; ++iChunk){
//		ljIdFromChunkNumber[iChunk] = lj.increment(ljIdFromChunkNumber[iChunk-1],chunkLevel);
//		range_faces[iChunk % ncolor].addRange(ljIdFromChunkNumber[iChunk],ljIdFromChunkNumber[iChunk]);
//		range->addRange(ljIdFromChunkNumber[iChunk],ljIdFromChunkNumber[iChunk]);
//	}
//	for( int ic = 0; ic < ncolor; ++ic ) {
//		viz->addHstmRangeFaces(&range_faces[ic],r[ic],g[ic],b[ic],0.0,1.0,0.01,&sIndex);
//	}
//	{
//		float re=0, ge=0.75, be=0;
//		viz->addHstmRange(range,re,ge,be,0.0,1.0,true,0.02,&sIndex);
//	}

	if(false) {
	// range->purge();
	int k = 0;
	float64 delta = 0.1;
	delta = 0.5;
	for( float64 lon=0; lon <= 360; lon += delta ) {
		for( float64 lat= -90; lat <= 90; lat += delta ) {

			STARE_ArrayIndexSpatialValue idx =
					index.ValueFromLatLonDegrees( lat, lon, resolutionLevel );

			// cout << k << ": " << hex << a << dec << endl << flush;
			// LatLonDegrees64 latlon = index.LatLonDegreesFromValue(a);

			// << setprecision(17)
			// << setw(20)
			// << scientific
			/*
		cout
		<< " stare = " << hex << a << dec
		<< ", latlon = " << latlon.lat << "," << latlon.lon
		<< endl << flush;
			 */

			EmbeddedLevelNameEncoding leftJustified;
			leftJustified.setIdFromSciDBLeftJustifiedFormat(idx);
			uint64 idx_gfx = leftJustified.getId();
			// range->reset();
			range->purge();
			range->addRange(idx_gfx,idx_gfx);
			++k;

			range_edges.addRange(idx_gfx,idx_gfx);

			float re=0, ge=0.75, be=0;

			viz->addHstmRange(range,re,ge,be,0.0,1.0,true,0.02,&sIndex);
			// viz->addHstmRange(range,0.0,0.75,0.125,0.0,1.0,true,0.02,&sIndex);

			// viz->addHstmRangeFaces(range,0.0,0.75,0.0,0.0,1.0,0.01,&sIndex);
			// int kc = k % ncolor;
			//

			// cout << 100 << endl << flush;
			STARE_ArrayIndexSpatialValue idx0 =
					index.ValueFromLatLonDegrees( lat, lon, 0 );

			// cout << 110 << endl << flush;
			leftJustified.setIdFromSciDBLeftJustifiedFormat(idx0);
			idx_gfx = leftJustified.getId();
			if(false) {
				cout << 120 << endl << flush;
				range0->purge();
				cout << " latlon = " << lat << "," << lon << ", idx0: " << hex << idx0 << ", idx: " << idx << ", idx_gfx: " << idx_gfx << dec << endl << flush;
				range0->addRange(idx_gfx,idx_gfx);
				cout << 130 << endl << flush;
				re=1, ge=1, be=0;
				viz->addHstmRange(range0,re,ge,be,0.0,1.0,true,0.04,&sIndex0);
				cout << 140 << endl << flush;
			} else {
				range0->addRange(idx_gfx,idx_gfx);
			}
		}
	}
	// cout << 500 << endl << flush;
	// range->reset();
	//	viz->addHstmRange(range,0.0,1.0,0.25,0.0,1.0,true,0.02,&sIndex);
	//	viz->addHstmRangeFaces(range,0.0,0.75,0.0,0.0,1.0,0.01,&sIndex);

	if( true ){
		// cout << 200 << endl << flush;
		float re=1, ge=1, be=0;
		viz->addHstmRange(range0,re,ge,be,0.0,1.0,true,0.04,&sIndex0);
		// cout << 210 << endl << flush;
	}
	}

	ok = true;
	return ok;
}
