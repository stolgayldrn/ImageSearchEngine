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


void WriteCSV(std::vector<std::vector<std::string>> testSet, std::string CSV_Path, int fileNum);
void WriteCSV(std::vector<std::vector<float>> dataVV, char* fileName, int fileNum);
void WriteCSV(std::vector<std::vector<int>> dataVV, char* fileName, int fileNum);
void ReleaseAll__ImRawIm(unsigned* vwI, json_t* myJSON, std::string& vwS, const char* ES_id);
void ReleaseAll__ImRawIm(unsigned* vwI, unsigned* vwI_low, json_t* myJSON, std::string& vwS, std::string& vwS_low, const char* ES_id);
void paramsConfig(Path &myPath, ELK_params &myES);
void ImageConfig(std::vector<std::string> dscList, int m, std::string imagePath, Image_Info &myIm, bool imgPath);
void ImportRawImage(Path myPath, ELK_params myES, TVoctreeVLFeat* VT, std::string dscPath, Image_Info myIm,
	uchar_descriptors * my_desc);
void ImportRawImage(Path myPath, ELK_params myES, TVoctreeVLFeat* VT, std::string dscPath, Image_Info myIm,
	uchar_descriptors * my_desc, TVoctreeVLFeat* VT_low, std::string dsc2Path, uchar_descriptors * my_desc2);
void QueryRawImage(Path myPath, ELK_params myES, TVoctreeVLFeat* VT, std::string dscPath, Image_Info myIm,
	uchar_descriptors *my_desc, std::vector<std::string> &testSet, std::vector<float> &scorePP, std::vector<float> &scoresELK);
int postProcess(uchar_descriptors *query, std::vector<std::string> dscV, std::vector<std::string> fnV, std::vector<float> &scoresV, int numELK);
void ImageSpliter(const cv::Mat Input, std::vector<cv::Mat> &OutputVector, std::vector<std::string> &OutputNames, int maxSize);
void scoreWeighting(std::vector<float> scoresPP, std::vector<float> scoresELK, std::vector<float> & scoresW);


#endif