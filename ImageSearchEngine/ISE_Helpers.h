#ifndef ISE_HELPERS_H
#define ISE_HELPERS_H

#include "helpers2.h"
#include "postProc.h"
#include <vector>
#include "descriptors.h"
#include <omp.h>
#include "ES_image.h"
#include "TVoctreeVLFeat.h"
#include <array>
#include "vld.h"

void WriteCSV(vector<vector<string>> testSet, char* fileName, int fileNum);
void WriteCSV(vector<vector<float>> dataVV, char* fileName, int fileNum);
void ReleaseAll__ImRawIm(unsigned* vwI, json_t* myJSON, string& vwS, const char* ES_id);
void ReleaseAll__ImRawIm(unsigned* vwI, unsigned* vwI_low, json_t* myJSON, string& vwS, string& vwS_low, const char* ES_id);
void paramsConfig(Path &myPath, ES_params &myES);
void ImageConfig(vector<string> dscList, int m, string imagePath, Image_Info &myIm);
void ImportRawImage(Path myPath, ES_params myES, TVoctreeVLFeat* VT, string dscPath, Image_Info myIm,
	uchar_descriptors * my_desc);
void ImportRawImage(Path myPath, ES_params myES, TVoctreeVLFeat* VT, string dscPath, Image_Info myIm,
	uchar_descriptors * my_desc, TVoctreeVLFeat* VT_low, string dsc2Path, uchar_descriptors * my_desc2);

void QueryRawImage(Path myPath, ES_params myES, TVoctreeVLFeat* VT, string dscPath, Image_Info myIm,
	uchar_descriptors *my_desc, vector<string> &testSet, vector<float> &scorePP, vector<float> &scoresELK, string &returnFileName);
int postProcess(uchar_descriptors *query, vector<string> dscV, vector<string> fnV, vector<float> &scoresV, int numEsReturns);

#endif