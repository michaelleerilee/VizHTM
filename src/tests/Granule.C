/*
 * Granule.C
 *
 *  Created on: Mar 18, 2019
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
using namespace std;


string file_directory = "/home/mrilee/data/STARE/";

/*
string files[12] =
{
		"MYD09.A2019003.2020.006.2019005020237.hdf",
		"MYD09.A2019003.2025.006.2019005020415.hdf",
		"MYD09.A2019003.2030.006.2019005020607.hdf",
		"MYD09.A2019003.2035.006.2019005020748.hdf",
		"MYD09.A2019003.2040.006.2019005020913.hdf",
		"MYD09.A2019003.2045.006.2019005021033.hdf",
		"MYD09.A2019003.2050.006.2019005021135.hdf",
		"MYD09.A2019003.2135.006.2019005015309.hdf",
		"MYD09.A2019003.2140.006.2019005015427.hdf",
		"MYD09.A2019003.2145.006.2019005015553.hdf",
		"MYD09.A2019003.2150.006.2019005015727.hdf",
		"MYD09.A2019003.2155.006.2019005015909.hdf"
};
*/

string files[22] =
{
		"MYD09.A2019003.1820.006.2019005015616.hdf",

		"MYD09.A2019003.1850.006.2019005020457.hdf",
		"MYD09.A2019003.1855.006.2019005020704.hdf",
		"MYD09.A2019003.1900.006.2019005020825.hdf",
		"MYD09.A2019003.1905.006.2019005020940.hdf",
		"MYD09.A2019003.1910.006.2019005021044.hdf",

		"MYD09.A2019003.2000.006.2019005015631.hdf", // 6
		"MYD09.A2019003.2005.006.2019005015815.hdf",
		"MYD09.A2019003.2010.006.2019005015932.hdf",
		"MYD09.A2019003.2015.006.2019005020109.hdf",
		"MYD09.A2019003.2020.006.2019005020237.hdf",
		"MYD09.A2019003.2025.006.2019005020415.hdf",
		"MYD09.A2019003.2030.006.2019005020607.hdf",
		"MYD09.A2019003.2035.006.2019005020748.hdf",
		"MYD09.A2019003.2040.006.2019005020913.hdf",
		"MYD09.A2019003.2045.006.2019005021033.hdf",
		"MYD09.A2019003.2050.006.2019005021135.hdf", // 16

		"MYD09.A2019003.2135.006.2019005015309.hdf",
		"MYD09.A2019003.2140.006.2019005015427.hdf",
		"MYD09.A2019003.2145.006.2019005015553.hdf",
		"MYD09.A2019003.2150.006.2019005015727.hdf",
		"MYD09.A2019003.2155.006.2019005015909.hdf"
};


bool Granule1( VizHTM *viz ) {

	cout << "Hello world" << endl;

	string swath1_key = "MODIS SWATH TYPE L2";

	// for( int ifile = 0; ifile < 5; ++ifile ) {
	// for( int ifile = 5; ifile < 6; ++ifile ) {
	for( int ifile = 6; ifile < 17; ++ifile ) {
	// for( int ifile = 7; ifile < 8; ++ifile ) {

		// char filename1[128] = "/home/mrilee/data/OMI/OMI-Aura_L2-OMAERO_2004m1001t0003-o01132_v003-2011m0926t170457.he5";

		// An hdf4 file...
		// char filename1[128] = "/home/mrilee/data/STARE/MYD09.A2019003.2040.006.2019005020913.hdf";
		string filename1 = file_directory + files[ifile];

		// char filename1[128] = "/home/mrilee/data/TestData/AMSR_E_L2A_BrightnessTemperatures_V09_200206190029_D.hdf";
		// string swath1_key = "Low_Res_Swath";

		/* Open using swath API */

		cout << "Opening swathfile " << filename1 << endl;
		int32 swathfile1;
		if ((swathfile1 = SWopen(filename1.c_str(), DFACC_RDONLY)) == -1) {
			fprintf(stderr, "error: cannot open swath '%s'\n",filename1.c_str());
			return -1;
		}

		/* Open a swath */
		int32 swath1;
		cout << "Attaching to " << swath1_key.c_str();
		if ((swath1 = SWattach(swathfile1, swath1_key.c_str())) == -1) {
			fprintf(stderr, "error: cannot attach to '%s'\n",swath1_key.c_str());
			return -1;
		} else {
			cout << " - ATTACHED" << endl << flush;
		}

		/*
	short 1km_Surface_Reflectance_Band_2(1km_Data_Lines=2030, 1km_Data_Samples=1354);
	  :long_name = "1km Surface Reflectance Band 2";
	  :units = "reflectance";
	  :Nadir_Data_Resolution = "1km";
	  :valid_range = -100S, 16000S; // short
	  :_FillValue = -28672S; // short
	  :scale_factor = 10000.0; // double
	  :scale_factor_err = 0.0; // double
	  :add_offset = 0.0; // double
	  :add_offset_err = 0.0; // double
	  :calibrated_nt = 5; // int
		 */
		// track vs. x-track?
		// was lines=1997, samples=243.
		int lines = 2030, samples = 1354;

		int32 datafield1rank;
		int32 datafield1dimsize[32];
		int32 datafield1type;
		char datafield1dimname[512];
		int16 *datafield1data;

		string variable_key = "1km Surface Reflectance Band 2";

		/* Retrieve information about '23.8H_Approx._Res.3_TB_(not-resampled)' datafield */
		if ((SWfieldinfo(swath1, variable_key.c_str(), &datafield1rank, datafield1dimsize, &datafield1type, datafield1dimname)) == -1) {
			fprintf(stderr, "error: cannot get the field info for '%s'\n",variable_key.c_str());
			return -1;
		} else {
			cout
			<< " variable_key: " << variable_key << endl
			<< "  datafield1rank:    " << datafield1rank << endl;
			for(int i=0; i < datafield1rank; ++i) {
				cout
				<< "  datafield1dimsize: " << i << ": " << datafield1dimsize[i] << endl;
			}
			cout
			<< "  datafield1type:    " << datafield1type << endl
			<< "  datafield1dimname: " << datafield1dimname << endl;

			/*
			 * Check hntdefs.h to identify what the type codes mean. datafield1type...
			 *  DFNT_FLOAT64  6
			 *  DFNT_INT16   22
			 *
			 */
		}

		lines   = datafield1dimsize[0];
		samples = datafield1dimsize[1];

		/* Allocate buffer for '23.8H_Approx._Res.3_TB_(not-resampled)' */
		if ((datafield1data = (int16 *)malloc(sizeof(int16) * lines * samples)) == NULL) {
			fprintf(stderr, "error: cannot allocate memory for '%s'\n",variable_key.c_str());
			return -1;
		}
		/* Read data from '23.8H_Approx._Res.3_TB_(not-resampled)' */
		if ((SWreadfield(swath1, variable_key.c_str(), NULL, NULL, NULL, datafield1data)) == -1) {
			fprintf(stderr, "error: cannot read field '%s'\n",variable_key.c_str());
			return -1;
		}
		/* Dump data from '23.8H_Approx._Res.3_TB_(not-resampled)'
		for (int i = 0; i < 10; ++i) {
			for (int j = 0; j < 2; ++j) {
				printf("%d ", datafield1data[j + 243 * i]);
			}
			printf("\n");
		}
		*/

		/*
	float Longitude(1km_Data_Lines=2030, 1km_Data_Samples=1354);
  :long_name = "Longitude";
  :units = "degrees_east";
  :Nadir_Data_Resolution = "1km";
  :valid_range = -180.0f, 180.0f; // float
  :_FillValue = 0.0f; // float
  :scale_factor = 1.0; // double
  :scale_factor_err = 0.0; // double
  :add_offset = 0.0; // double
  :add_offset_err = 0.0; // double
  :calibrated_nt = 5; // int
  :_CoordinateAxisType = "Lon";
		 */

		int32 geofield1rank;
		int32 geofield1dimsize[32];
		int32 geofield1type;
		char geofield1dimname[512];
		float32 *geofield1data;

		/* Retrieve information about 'Longitude' geolocation field */
		if ((SWfieldinfo(swath1, "Longitude", &geofield1rank, geofield1dimsize, &geofield1type, geofield1dimname)) == -1) {
			fprintf(stderr, "error: cannot get the field info for 'Longitude'\n");
			return -1;
		}
		/* Allocate buffer for 'Longitude' */
		// if ((geofield1data = malloc(sizeof(float32) * 1997 * 243)) == NULL) {
		if ((geofield1data = (float32 *)malloc(sizeof(float32) * lines * samples )) == NULL) {
			fprintf(stderr, "error: cannot allocate memory for 'Longitude'\n");
			return -1;
		}
		/* Read data from 'Longitude' */
		if ((SWreadfield(swath1, "Longitude", NULL, NULL, NULL, geofield1data)) == -1) {
			fprintf(stderr, "error: cannot read field 'Longitude'\n");
			return -1;
		}
		/* Dump data from 'Longitude'
		for (int i = 0; i < 10; ++i) {
			for (int j = 0; j < 2; ++j) {
				printf("%f ", geofield1data[j + samples * i]);
			}
			printf("\n");
		}
		*/

		/*
		 *
		 */

		int32 geofield2rank;
		int32 geofield2dimsize[32];
		int32 geofield2type;
		char geofield2dimname[512];
		float32 *geofield2data;

		/* Retrieve information about 'Latitude' geolocation field */
		if ((SWfieldinfo(swath1, "Latitude", &geofield2rank, geofield2dimsize, &geofield2type, geofield2dimname)) == -1) {
			fprintf(stderr, "error: cannot get the field info for 'Latitude'\n");
			return -1;
		}
		/* Allocate buffer for 'Latitude' */
		if ((geofield2data = (float32 *)malloc(sizeof(float32) * lines * samples)) == NULL) {
			fprintf(stderr, "error: cannot allocate memory for 'Latitude'\n");
			return -1;
		}
		/* Read data from 'Latitude' */
		if ((SWreadfield(swath1, "Latitude", NULL, NULL, NULL, geofield2data)) == -1) {
			fprintf(stderr, "error: cannot read field 'Latitude'\n");
			return -1;
		}
		/* Dump data from 'Latitude'
		for (int i = 0; i < 10; ++i) {
			for (int j = 0; j < 2; ++j) {
				printf("%f ", geofield2data[j + samples * i]);
			}
			printf("\n");
		}
		*/

		/*
		 *
		 */

		// cout << 100 << endl << flush;
		STARE index;
		HstmRange* range = new HstmRange;
		// int resolutionLevel = 5;
		int resolutionLevel = 6;
		// int resolutionLevel = 7;
		// int resolutionLevel = 8;
		// int resolutionLevel = 14;
		SpatialIndex sIndex = index.getIndex(resolutionLevel);

		/*
		cout << 200 << endl << flush;
		cout
		<< setprecision(17)
		<< setw(20)
		<< scientific;
		*/

		range->purge();
		int k = 0;
		// for( int iline = 0; iline < lines; ++iline ) {
		//	for( int isample = 0; isample < samples; ++isample ) {
		for( int iline = 0; iline < lines; iline += 10 ) {
			for( int isample = 0; isample < samples; isample += 10 ) {

				// cout << 300 << " " << k << " iline,isample = " << iline << "," << isample << endl << flush;
				/*
			cout << k << " " << isample << " " << iline << " kij lat,lon = "
			<< geofield2data[isample + samples * iline] << "," << geofield1data[isample + samples * iline]
			<< endl << flush;
				 */

				STARE_ArrayIndexSpatialValue a =
						index.ValueFromLatLonDegrees(
								geofield2data[isample + samples * iline],
								geofield1data[isample + samples * iline],
								resolutionLevel);
				// cout << k << ": " << hex << a << dec << endl << flush;
				LatLonDegrees64 latlon = index.LatLonDegreesFromValue(a);

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
				leftJustified.setIdFromSciDBLeftJustifiedFormat(a);
				uint64 b = leftJustified.getId();
				// range->reset();
				range->addRange(b,b);
				++k;
			}
		}
		// cout << 500 << endl << flush;
		// range->reset();
		viz->addHstmRange(range,0.0,1.0,0.25,0.0,1.0,true,0.02,&sIndex);

		viz->addHstmRangeFaces(range,0.0,0.75,0.0,0.0,1.0,0.01,&sIndex);
		// cout << 600 << endl << flush;
		/*
	void VizHTM::addHstmRange(
			HstmRange *range,
			float r, float g, float b, float a, float scale, bool arcFlag
		 */

		/*
		 *
		 */

		/* Release the buffer for '23.8H_Approx._Res.3_TB_(not-resampled)' */
		free(datafield1data);
		/* Release the buffer for 'Longitude' */
		free(geofield1data);
		/* Release the buffer for 'Latitude' */
		free(geofield2data);

		/*
		 *
		 */

		/* Close the swath named 'Low_Res_Swath' */
		cout << "Detaching..." << endl;
		if ((SWdetach(swath1)) == -1) {
			fprintf(stderr, "error: cannot detach from '%s'\n",swath1_key.c_str());
			return -1;
		}
		cout << "Closing swathfile " << filename1 << endl;
		/* Close */
		if ((SWclose(swathfile1)) == -1) {
			fprintf(stderr, "error: cannot close swath\n");
			return -1;
		}
	}
	cout << "Done...";
	return 0;
}



