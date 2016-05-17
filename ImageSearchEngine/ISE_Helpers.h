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
#include <ctime>
#include <chrono>

#define  MIN_GOOD_MATCHES 10
#define GOOD_MATCHES_DISTANCE__AKAZE 150.0
#define GOOD_MATCHES_DISTANCE__EZSIFT 80.0
//void WriteCSV(std::vector<std::vector<std::string>> testSet, std::string CSV_Path, int fileNum);
//void WriteCSV(std::vector<std::vector<float>> dataVV, char* fileName, int fileNum);
//void WriteCSV(std::vector<std::vector<int>> dataVV, char* fileName, int fileNum);
void ReleaseAll__ImRawIm(unsigned* vwI, json_t* myJSON, std::string& vwS, const char* ES_id);
void paramsConfig(Path &myPath, ELK_params &myES);
void ImageConfig(std::vector<std::string> dscList, int m, std::string imagePath, Image_Info &myIm, bool imgPath);
void ImageConfig(std::string imagePath, std::string imageFileName, Image_Info& myIm, bool imgPath);
void ImportRawImage(Path myPath, ELK_params myES, TVoctreeVLFeat* VT, std::string dscPath, Image_Info myIm,
	uchar_descriptors * my_desc);
void scoreNormELK(std::vector<float>& scoresELK, int totalNumELK, std::vector<float>& normELKScore);
void QueryRawImage(Path myPath, ELK_params myES, TVoctreeVLFeat* VT, std::string dscPath, Image_Info myIm,
                   uchar_descriptors *my_desc, std::vector<std::string> & returnedFileNames, std::vector<float> &scorePP, std::vector<float> &scoresELK);
int postProcess(uchar_descriptors *query, std::vector<std::string> dscV, std::vector<std::string> fnV, std::vector<float> &scoresV, int numELK);
void ImageSpliter(const cv::Mat Input, std::vector<cv::Mat> &OutputVector, std::vector<std::string> &OutputNames, int maxSize);
void scoreWeighting(std::vector<float> scoresPP, std::vector<float> scoresELK, std::vector<float> & scoresW);
std::string currentDateTime();
void printErrorToLog(Path myPath, std::ofstream &ofs, std::vector<std::string> imgList, int i, std::basic_string<char> imgPathDUMP, uchar_descriptors myDesc, std::string errCmd);
void findMatches(uchar_descriptors &descriptor_1, uchar_descriptors &descriptor_2, std::vector<DMatch >& good_matches);
void findIntersectedFeatures(std::string imgPath, cv::Mat img1, uchar_descriptors& descriptor_1, std::vector<DMatch >& inMatches);

#endif