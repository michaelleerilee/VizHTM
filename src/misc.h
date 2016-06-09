/*
 * misc.h
 *
 *  Created on: Mar 11, 2016
 *      Author: mrilee
 */

#ifndef MISC_H_
#define MISC_H_

#include <string>
#include <regex>
#include <vector>

std::string executableNameFromPathAsCString(const char* arg0);
std::string formattedZeroPaddedInteger(int iOut, const int fieldWidth=3, const char padChar='0');
std::string formattedOutFileName(
		const std::string baseName="output",
		const std::string itemName="",
		const std::string extension=".tmp");
std::string formattedDateTime(std::string format="%Y-%m%d-%H%M%S-%Z");

#endif /* MISC_H_ */
