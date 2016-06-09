/*
Copyright (C) 2015-20 S.Tolga Yildiran.
All rights reserved.

This file is part of Tolga Yildiran Video Search library and is made available under
the terms of the BSD license (see the COPYING file).
*/
/************************************************************************/
/* Tolga Yildiran														*/
/* 24/05/2016 - 2020															*/
/* stolgayldrn@gmail.com												*/
/************************************************************************/

#ifndef ISE_HELPERS_H
#define ISE_HELPERS_H

#include "helpers2.h"
#include "postProc.h"
#include <vector>
#include "uchar_descriptors.h"
#include <omp.h>
#include "ES_image.h"
#include "TVoctreeVLFeat.h"
#include <array>
#include "vld.h"
#include <ctime>
#include <chrono>

#define  MIN_GOOD_MATCHES 10
#define GOOD_MATCHES_DISTANCE__AKAZE 150.0
#define GOOD_MATCHES_DISTANCE__EZSIFT 80.0
//void writeToCSV(std::vector<std::vector<std::string>> testSet, std::string CSV_Path, int fileNum);
//void writeToCSV(std::vector<std::vector<float>> dataVV, char* fileName, int fileNum);
//void writeToCSV(std::vector<std::vector<int>> dataVV, char* fileName, int fileNum);
void releaseMemoryForRawIm(unsigned* vwI, json_t* myJSON, std::string& vwS, const char* ES_id);
void initializeParameters(Path &myPath, ELKParams &myES);
void configureImageConfig(std::vector<std::string> dscList, int m, std::string imagePath, ImageInfo &myIm, bool imgPath);
void configureImageConfig(std::string imagePath, std::string imageFileName, ImageInfo& myIm, bool imgPath);
void importRawImage(Path myPath, ELKParams myES, TVoctreeVLFeat* VT, std::string dscPath, ImageInfo myIm,
	UcharDescriptors * my_desc);
void calculateNormalizedELKScore(std::vector<float>& scoresELK, int totalNumELK, std::vector<float>& normELKScore);
void queryRawImage(Path myPath, ELKParams myES, TVoctreeVLFeat* VT, std::string dscPath, ImageInfo myIm,
                   UcharDescriptors *my_desc, std::vector<std::string> & returnedFileNames, std::vector<float> &scorePP, std::vector<float> &scoresELK);
int postProcess(UcharDescriptors *query, std::vector<std::string> fileNameVec, std::vector<float> & scoresV, int numELK);
void splitUpNewspaperImages(const cv::Mat Input, std::vector<cv::Mat> &OutputVector, std::vector<std::string> &OutputNames, int maxSize);
void scoreWeighting(std::vector<float> scoresPP, std::vector<float> scoresELK, std::vector<float> & scoresW);
std::string currentDateTime();
void printErrorToLog(Path myPath, std::ofstream &ofs, std::basic_string<char> imgPathDUMP, std::string errCmd);

#endif