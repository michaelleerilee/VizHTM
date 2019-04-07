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

#include <fstream>
#include <iostream>

#include <array>

using namespace std;

void plotBlockingSphere(VizHTM* viz, float r, float g, float b, float radius);
void testShapeFiles(VizHTM *viz, float r, float g, float b, float deltaZ);

void AllTrixelsAtLevelAcrossBins(int level, int nbins, HstmRange ranges[], bool randomized) {
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
		int ic = -1;
		if(randomized) { ic = uniformDouble()*nbins; } else { ic = (iChunk % nbins); }
		ranges[ ic ].addRange(ljId,ljId);
		if(false){
			cout << iChunk << " ranges loop nranges " << ranges[ic].range->nranges() << endl << flush;
		}
		ljId = lj.increment(ljId,level);
	}
	if(false){
		for( int ir = 0; ir < nbins; ++ir) {
			cout << ir << " ranges ret nranges " << ranges[ir].range->nranges() << endl << flush;
		}
	}
}

bool Area1( VizHTM *viz, int level0 ) {
	bool ok = false;

	cout << "Area1 level0 = " << level0 << endl;

	/*
	std::ofstream noteFile;
	noteFile.open( base+"/notes.org" );
	noteFile << notes->str();
	noteFile.close();
	*/
	bool areaFile_flag = true;
	ofstream areaFile;
	if(areaFile_flag) {
		areaFile.open("areas-level="+to_string(level0)+".csv");
		areaFile << "area,count,min area,max area,avg area,std dev,total area,total area/4" << endl;
	}

	bool noteStream_flag = false;

	if(noteStream_flag) {
		if( !viz->notestream() ) {
			viz->addNotes(new stringstream);
		}
		*(viz->notestream()) << "Area1 Start" << endl;
		*(viz->notestream()) << "level0 = " << level0 << endl << flush;
	}
	// int level0 = 2;
	int nbins = 1;
	HstmRange *range;
	HstmRange ranges[nbins];
	AllTrixelsAtLevelAcrossBins(level0,nbins,ranges,false);
	range = &(ranges[0]);

	STARE index;
	SpatialIndex sIndex(level0);
	EmbeddedLevelNameEncoding leftJustified;
	KeyPair kp; int indexp;
	range->reset();
	// int count = 0;
	float64 total_area = 0, min_area = 100, max_area = -1;
	float64 average_area = 0, std_deviation2 = 0;
	if(noteStream_flag) { *(viz->notestream()) << "Area1-100" << " range size = " << range->range->nranges() << endl << flush; }
	range->reset();
	int count = 0;
	while((indexp = range->getNext(kp)) > 0) {
		leftJustified.setId(kp.lo);
		int level = leftJustified.getLevel();
		string loName = leftJustified.getName();
		uint64 termId = leftJustified.idFromTerminatorAndLevel_NoDepthBit(kp.hi,level);
		leftJustified.setId(termId);
		string hiName = leftJustified.getName();
		// cout << "099" << " count = " << ++count << endl << flush;
		// cout << 100 << " lo,hi name: " << loName << " " << hiName << endl << flush;
		// cout << 101 << " level:      " << level << endl << flush;
		// cout << 102 << " kp:         " << hex << kp.lo << " " << kp.hi << dec << endl << flush;
		// cout << 103 << " sI mxlevel: " << sIndex.getMaxlevel() << endl << flush;
		if( sIndex.getMaxlevel() != level ) {
			sIndex.setMaxlevel(level);
		}
		uint64 id0 = sIndex.idByName(loName.c_str());
		uint64 id1 = sIndex.idByName(hiName.c_str());
		// cout << 103 << " id0,id1 = " << hex << id0 << " " << id1 << dec << endl << flush;
		for(uint64 id=id0; id<=id1; id++) {
			float64 area = sIndex.areaByHtmId(id);
			total_area += area;
			min_area = (min_area < area) ? min_area : area;
			max_area = (max_area > area) ? max_area : area;
			float64    a0 = average_area;
			average_area  = a0 + ((area-a0)/(count+1));
			float64    q0 = std_deviation2;
			std_deviation2 = q0 + (area-a0)*(area-average_area);
			if(noteStream_flag) { *(viz->notestream()) << "area: " << area << endl << flush; }
			if(areaFile_flag) { areaFile << area << endl; }
			++count;
		}
	}
	if(noteStream_flag) {
		*(viz->notestream())
					<< "count:         " << count << endl << flush
					<< "min,max area:  " << min_area << "," << max_area << endl << flush
					<< "average area:  " << total_area/range->range->nranges() << endl << flush
					<< "std deviation: " << sqrt(std_deviation2)/(count-1) << endl << flush
					<< "total_area:    " << total_area << endl << flush
					<< "total_area/4:  " << 0.25*total_area << endl << flush
					<< "Area1 done" << endl << flush;
	}

	if(areaFile_flag) {
		areaFile
		<< ","
		<< count << ","
		<< min_area << ","
		<< max_area << ","
		<< total_area/range->range->nranges() << ","
		<< sqrt(std_deviation2)/(count-1)  << ","
		<< total_area << ","
		<< total_area*0.25
		<< endl
		;
		areaFile.close();
	}

	ok = true;
	return ok;
}

bool Chunk1( VizHTM *viz ) {
	bool ok = false;

	cout << "Chunk1" << endl;

	if( !viz->notestream() ) {
		viz->addNotes(new stringstream);
	}
	*(viz->notestream()) << "Chunk1 Start" << endl;

	// MLR Doesn't work Fedora+GLX viz->lineWidth = 1;

	testShapeFiles(viz,0,1,1,0.035);
	// float re=1.0, ge=1.0, be=0.5, ae=0.0, zlayer = 0.03, rgb_scale = 0.7;
	float re=1.0, ge=1.0, be=0.5, ae=0.0, zlayer = 0.035, rgb_scale = 0.7;

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

	int  ncolor      = 16;
	int  chunkLevel  = 3; // faces: chunks at level
	bool randomized  = true;
	// bool randomized  = false;
	int  nbins       = 1;
	int  nlevels     = 5; // grid: nlevels=3 => levels 0,1,2


	*(viz->notestream()) << "ncolor = " << ncolor << endl;
	cout << "ncolor: " << ncolor << endl << flush;

	*(viz->notestream()) << "chunkLevel = " << chunkLevel << endl << flush;
	cout << "chunkLevel, nChunks: " << chunkLevel << ", " << 8*int(pow(4,chunkLevel)) << endl << flush;

	float chunkTransparency = 0.2;
	// float chunkTransparency = 0.0;
	*(viz->notestream()) << "chunkTransparency = " << chunkTransparency << endl << flush;

	// float color_scale = 0.5;
	// float color_scale = 0.75;
	// float color_scale = 0.9;
	float color_scale = 1.0;
	*(viz->notestream()) << "chunk color_scale = " << color_scale << endl << flush;

	float r[ncolor],g[ncolor],b[ncolor],a[ncolor];
	for( int ic = 0; ic < ncolor; ++ic ) {
		r[ic] = color_scale*ic/(ncolor-1.0);
		g[ic] = color_scale*ic/(ncolor-1.0);
		b[ic] = color_scale*ic/(ncolor-1.0);
		a[ic] = chunkTransparency;
		// b[ic] = color_scale*(ncolor-ic-1)/(ncolor-1.0);
	}
	HstmRange range_edges;
	HstmRange range_faces[ncolor];

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
		*(viz->notestream()) << endl << "Drawing faces with ncolors bin colors" << endl;
		HstmRange ranges[ncolor];
		AllTrixelsAtLevelAcrossBins(chunkLevel,ncolor,ranges,randomized);
		SpatialIndex sIndex = index.getIndex(chunkLevel);
		if(randomized) {
			*(viz->notestream()) << "chunks randomized: uniformDouble" << endl;
		}
		// float re=0, ge=0.75, be=0;
		if(true) {
			for( int ic = 0; ic < ncolor; ++ic ) {
				// cout << ic << " ranges size " << ranges[ic].range->nranges() << endl << flush;
				int ic0 = ic;
//				if(randomized) {
//					ic0 = uniformDouble()*ncolor;
//				}
//				/*
//				 * Verify randomness
//				 */
//				cout << "ic,ic0: " << ic << "," << ic0 << endl << flush;
//				/**/
				viz->addHstmRangeFaces(&(ranges[ic]),r[ic0],g[ic0],b[ic0],a[ic],1.0,0.01,&sIndex);
				// viz->addHstmRangeFaces(&(ranges[ic]),r[ic],g[ic],b[ic],0.0,1.0,0.01,&sIndex);
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
		*(viz->notestream()) << endl;
		*(viz->notestream()) << "Drawing trixels at various levels" << endl;
		*(viz->notestream()) << "nbins is the number of bins into which to put trixels" << endl;
		*(viz->notestream()) << "nlevels is the number of levels to draw. levels=[0,1,2,...,nlevels-1]" << endl << endl;
		*(viz->notestream()) << "nbins = " << nbins << endl;
		*(viz->notestream()) << "nlevels = " << nlevels << endl;
		HstmRange ranges0[nbins];
		for(int level=0; level<nlevels; ++level) {
			for(int ib=0; ib < nbins; ++ib) {
				ranges0[ib].purge();
			}
			AllTrixelsAtLevelAcrossBins(level,nbins,ranges0,false);
			SpatialIndex sIndex0 = index.getIndex(level);
			for( int ic = 0; ic < nbins; ++ic ) {
				// cout << "ranges0 size " << ranges0[ic].range->nranges() << endl << flush;
				// be = (1.0*ic)/nbins;
				viz->addHstmRange(&(ranges0[ic]),re,ge,be,ae,1.0,true,zlayer,&sIndex0);
				// re *= 0.75; ge *= 0.75; be *= 0.75; zlayer -= 0.001;
				re *= rgb_scale; ge *= rgb_scale; be *= rgb_scale; zlayer -= 0.003;
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

	*(viz->notestream()) << "Chunk1 End" << endl;

	ok = true;
	return ok;
}
