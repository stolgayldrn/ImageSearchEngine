//#include "helpers2.h"
//#include "postProc.h"
//#include <vector>
//#include "descriptors.h"
//#include <omp.h>
//#include "ES_image.h"
//#include "TVoctreeVLFeat.h"
//#include <array>
//#include "vld.h"
#include "ISE_Helpers.h"

using namespace std;
using namespace cv;
//
//void WriteCSV(vector<vector<string>> testSet, char* fileName);
//void ReleaseAll__ImRawIm(unsigned* vwI, json_t* myJSON, string& vwS, const char* ES_id);
//void ReleaseAll__ImRawIm(unsigned* vwI, unsigned* vwI_low, json_t* myJSON, string& vwS, string& vwS_low, const char* ES_id);
//void paramsConfig(Path &myPath, ES_params &myES);
//void ImageConfig(vector<string> dscList, int m, string imagePath, Image_Info &myIm);
//void ImportRawImage(Path myPath, ES_params myES, TVoctreeVLFeat* VT, string dscPath, Image_Info myIm, 
//	uchar_descriptors * my_desc);
//void ImportRawImage(Path myPath, ES_params myES, TVoctreeVLFeat* VT, string dscPath, Image_Info myIm,
//	uchar_descriptors * my_desc, TVoctreeVLFeat* VT_low, string dsc2Path, uchar_descriptors * my_desc2);
//
//void QueryRawImage(Path myPath, ES_params myES, TVoctreeVLFeat* VT, string dscPath, Image_Info myIm,
//	uchar_descriptors * my_desc, vector<string> &testSet, string &returnFileName);
//int post_process_step(uchar_descriptors *query, vector<string> dscV, vector<string> fnV, vector<float> scoresV, int numEsReturns);

int main()
{
	Path myPath;
	ES_params myES;
	TVoctreeVLFeat* VT = new TVoctreeVLFeat;
	TVoctreeVLFeat* VT_low = new TVoctreeVLFeat;
	vector<string> imgList, dscList, foldList;
	vector<vector<string>> testSet;
	long start = clock();
	double duration;
	unsigned int err_count = 0;

	paramsConfig(myPath, myES);
	VocTreeInit(VT, VT_low, myPath);
	GET_FolderList(myPath.DataSet.c_str(), foldList);

	/*
	for (unsigned int f = 0; f < foldList.size(); f++)
	{
		vector <string> subfoldList;
		String dscFoldPath, imgFoldPath, dscFold2Path;
		GET_FolderList((myPath.DataSet + "/" + foldList[f] + "/" + myPath.dscFoldName).c_str(), subfoldList);
		//GET_FolderList((myPath.DataSet + "/" + foldList[f] + "/" + myPath.imgFoldName).c_str(), subfoldList);
		for (unsigned int sf = 0; sf < subfoldList.size(); sf++)
		{
			dscFoldPath = myPath.DataSet + "/" + foldList[f] + "/" + myPath.dscFoldName + "/" + subfoldList[sf];;
			dscFold2Path = myPath.DataSet + "/" + foldList[f] + "/" + myPath.dscFoldName2 + "/" + subfoldList[sf];;
			imgFoldPath = myPath.DataSet + "/" + foldList[f] + "/" + myPath.imgFoldName + "/" + subfoldList[sf];;
			printf("\nMain:::extracted akaze for folder:%s", subfoldList[sf].c_str());
			//GET_DirectoryImages(imgFoldPath.c_str(), imgList);
			//PathControl(dscFoldPath);
			//PathControl(dscFold2Path);
			GET_DirectoryDSCs(dscFoldPath.c_str(), dscList);

			omp_set_dynamic(0);     // Explicitly disable dynamic teams
			omp_set_num_threads(8); // Use 4 threads for all consecutive parallel regions
			//for (int m = 0; m < imgList.size(); m++)
			//#pragma omp parallel for
			for (int m = 0; m < dscList.size(); m++)
			{
				string imgPath = imgFoldPath + "/" + dscList[m];
				string dscPath = dscFoldPath + "/" + dscList[m];
				//string dsc2Path = dscFold2Path + "/" + imgList[m] + ".dsc";
				Image_Info myIm;
				ImageConfig(dscList, m, dscPath.c_str(), myIm);

				if (IS_DscFile(dscPath.c_str()))
				{
					try
					{
						uchar_descriptors my_desc(imgPath.c_str(), dscPath.c_str(), AKAZE_FEATS);
						my_desc.read_dsc();
						//my_desc.extract_AKAZE_feats();
						
						//uchar_descriptors my_desc2(imgPath.c_str(), dsc2Path.c_str(), AKAZE_FEATS);
						//my_desc2.extract_AKAZE_low_feats();

						//ImportRawImage(myPath, myES, VT, dscPath, myIm, &my_desc);

					}
					catch (exception e)
					{
						printf("\nMain:::extract akaze and write:%s", e.what());
					}
				}

				if (m % 100 == 0)
				{
					double k = dscList.size() / 100;
					double num = (m/k);
					printf("\rProcess Rate = %4f\n", num);
				}
			}
			imgList.clear();
			imgList.shrink_to_fit();
		}
		subfoldList.clear();
		subfoldList.shrink_to_fit();
	}
	foldList.clear();
	foldList.shrink_to_fit();
	//*/
	
	///*
	for (unsigned int f = 0; f < foldList.size(); f++)
	{
		String imgFoldPath;
		imgFoldPath = myPath.DataSet + "/" + foldList[f]  ;
		GET_DirectoryImages(imgFoldPath.c_str(), imgList);
		vector<string> testSetFold;
		testSetFold.push_back(foldList[f] + ".jpg");
		for (int i = 0; i < 11*25; i++)
			testSetFold.push_back("-1");

			omp_set_dynamic(0);     // Explicitly disable dynamic teams
			omp_set_num_threads(8); // Use 4 threads for all consecutive parallel regions
			//#pragma omp parallel for
			for (int m = 0; m < imgList.size(); m++)
			{
				string imgPath = imgFoldPath + "/" + imgList[m];
				Image_Info myIm;
				ImageConfig(imgList, m, imgPath.c_str(), myIm);
				string returnImage;

				if (IS_ImageFile(imgPath.c_str()))
				{
					try
					{
						vector<string> testIm;
						uchar_descriptors my_desc(imgPath.c_str(), "NoDSC", AKAZE_FEATS);
						my_desc.extract_AKAZE_feats();
						//uchar_descriptors my_desc2(imgPath.c_str(), "NoDSC", AKAZE_FEATS);
						//my_desc2.extract_AKAZE_low_feats();

						QueryRawImage(myPath, myES, VT, "NoDSC", myIm, my_desc, testIm, returnImage);
						if (testIm.size() == 11)
						{
							for (int i = 0; i < 11; i++)
								testSetFold[ (m*11) + i + 1 ] = (testIm[i]);
							
							testIm.clear();
							testIm.shrink_to_fit();
						}
						else
						{
							for (int i = 0; i < 11; i++)
								testSetFold.push_back("-1");
						}
					}
					catch (exception e)
					{
						printf("\nMain:::extract akaze and write:%s", e.what());
					}
				}
			}
			testSet.push_back(testSetFold);
			testSetFold.clear();
			testSetFold.shrink_to_fit();
			imgList.clear();
			imgList.shrink_to_fit();
			if (f % 5 == 0)
			{
				double num = (f*100.00 / 1000);
				printf("\rProcess Rate = %.2f\n", num);
			}
	}
	foldList.clear();
	foldList.shrink_to_fit();
	WriteCSV(testSet,"flicker10K_test_akaze__33_LEAK.csv");
	//*/
	duration = clock() - start;
	delete VT;
	delete VT_low;
	printf(":::Duration: %.2f, Error count: %d", duration, err_count);
	waitKey(0);
	return 0;
}
//
//void ReleaseAll__ImRawIm(unsigned* vwI, json_t* myJSON, string& vwS, const char* ES_id)
//{
//	delete[] vwI;
//	json_decref(myJSON);
//	vwS = "";
//	delete ES_id;
//}
//
//void ReleaseAll__ImRawIm(unsigned* vwI, unsigned* vwI_low, json_t* myJSON, string& vwS, string& vwS_low, const char* ES_id)
//{
//	delete[] vwI;
//	delete[] vwI_low;
//	json_decref(myJSON);
//	vwS = "";
//	vwS_low = "";
//	delete ES_id;
//}
//
//void WriteCSV(vector<vector<string>> testSet, char* fileName)
//{
//	ofstream myfile;
//	string CSV_Path = fileName;
//	myfile.open(CSV_Path.c_str());
//	for (int i = 0; i < testSet.size(); i++){
//		string lineStr = "";
//		for (int ii = 0; ii < 276 ; ii++){
//			lineStr += testSet[i][ii] + ";";
//		}
//		lineStr += "\n";
//		myfile << lineStr;
//	}
//}
//
//void paramsConfig(Path &myPath, ES_params &myES)
//{
//	myPath.dscFoldName = "dsc_akaze2";
//	myPath.dscFoldName2 = "dsc_akaze_low";
//	myPath.DataSet = "C:/ImageSearch/TEST_IM_33";
//	myPath.imgFoldName = "images";
//	myPath.subFolderingLevel = 2;
//	myPath.VocTree = "C:/ImageSearch/VT_Trees/VT_flicker500K_AKAZE_middle_tree_S2_P";
//	myPath.VocTreeLow = "C:/ImageSearch/VT_Trees/VT_flicker500K_AKAZE_small_tree_S2_P";
//
//	myES.index = "flicker1m_test2";
//	myES.type = "akaze";
//	myES.url = "http://172.16.10.202:9200";
//	myES.userPWD = "";
//}
//
//void ImageConfig(vector<string> vList, int m, string imagePath, Image_Info& myIm)
//{
//	myIm.dataSet = "flicker1M";
//	myIm.dataSubSet = "";
//	myIm.descriptorType = "akaze";
//	myIm.encoding = "jpg";
//	//myIm.fileName = vList[m];
//	myIm.fileName = vList[m].substr(0, vList[m].length()-4);
//	myIm.height = 0.0;
//	myIm.width = 0.0;
//	myIm.Import = "true";
//	myIm.Query = "false";
//	myIm.path = imagePath;
//	myIm.source_type = "flicker1M";
//}
//
//void ImportRawImage(Path myPath, ES_params myES, TVoctreeVLFeat* VT, string dscPath, 
//	Image_Info myIm, uchar_descriptors * my_desc)
//{
//	unsigned int * vwI = new unsigned int[my_desc->get_num_descriptors()];
//	json_t* myJSON = json_object();
//	string vwS = "";
//	const char * ES_id = new    char;
//
//	if (my_desc->get_num_descriptors() > 0)
//	{
//		try
//		{
//			VT->quantize_multi(vwI, my_desc->get_data(), my_desc->get_num_descriptors(), 61);
//			for (unsigned int s = 0; s < my_desc->get_num_descriptors(); s++)
//				vwS += " " + int2string(int(vwI[s]));
//
//			if (vwS != "")
//			{
//				getJSON_new_image(myIm, myPath, myJSON, vwS);
//				ES_commit(myES, myJSON, ES_id, myIm.fileName.c_str());
//			}
//		}
//		catch (exception e)
//		{
//			printf("\nElasticSearch:::commit error:%s", e.what());
//		}
//	}
//	try
//	{
//		ReleaseAll__ImRawIm(vwI, myJSON, vwS, ES_id);
//	}
//	catch (exception e)
//	{
//		printf("\nElasticSearch:::release error:%s", e.what());
//	}
//}
//
//void ImportRawImage(Path myPath, ES_params myES, TVoctreeVLFeat* VT, string dscPath, Image_Info myIm, 
//	uchar_descriptors * my_desc, TVoctreeVLFeat* VT_low, string dsc2Path, uchar_descriptors * my_desc2)
//{
//	unsigned int * vwI = new unsigned int[my_desc->get_num_descriptors()];
//	unsigned int * vwI_low = new unsigned int[my_desc2->get_num_descriptors()];
//	json_t* myJSON = json_object();
//	string vwS = "";
//	string vwS_low = "";
//	const char * ES_id = new char;
//
//	if (my_desc->get_num_descriptors() > 0 && my_desc2->get_num_descriptors() > 0)
//	{
//		try
//		{
//			VT->quantize_multi(vwI, my_desc->get_data(), my_desc->get_num_descriptors(), 61);
//			VT_low->quantize_multi(vwI_low, my_desc2->get_data(), my_desc2->get_num_descriptors(), 61);
//			for (unsigned int s = 0; s < my_desc->get_num_descriptors(); s++)
//				vwS += " " + int2string(int(vwI[s]));
//
//			for (unsigned int s = 0; s < my_desc2->get_num_descriptors(); s++)
//				vwS_low += " " + int2string(int(vwI_low[s]));
//
//			if (vwS != "")
//			{
//				getJSON_new_image(myIm, myPath, myJSON, vwS, vwS_low);
//				ES_commit(myES, myJSON, ES_id, myIm.fileName.c_str());
//			}
//		}
//		catch (exception e)
//		{
//			printf("\nElasticSearch:::commit error:%s", e.what());
//		}
//	}
//	try
//	{		
//		ReleaseAll__ImRawIm(vwI, vwI_low, myJSON, vwS, vwS_low, ES_id);
//	}
//	catch (exception e)
//	{
//		printf("\nElasticSearch:::release error:%s", e.what());
//	}
//}
//
//void QueryRawImage(Path myPath, ES_params myES, TVoctreeVLFeat* VT, string dscPath, Image_Info myIm,
//	uchar_descriptors * my_desc, vector<string> &testSet, string &returnFileName)
//{
//	unsigned int * vwI = new unsigned int[my_desc->get_num_descriptors()];
//	json_t* myJSON = json_object();
//	string vwS = "";
//	const char * ES_id = new char;
//	vector<string> dscPathsV;
//	vector<float> scoresV;
//
//	if (my_desc->get_num_descriptors() > 0 )
//	{
//		try
//		{
//			VT->quantize_multi(vwI, my_desc->get_data(), my_desc->get_num_descriptors(), 61);
//			for (unsigned int s = 0; s < my_desc->get_num_descriptors(); s++)
//				vwS += " " + int2string(int(vwI[s]));
//
//			if (vwS != "")
//			{
//				printf("nok \n");
//				getJSON_query_image(myJSON, vwS, "words_string");
//				//TODO: add scores at ES_post_query
//				ES_post_query(myES, myJSON, myIm, testSet, dscPathsV, scoresV);
//				post_process_step(my_desc, dscPathsV, testSet, scoresV, 10);
//
//				printf("ok \n");
//			}
//		}
//		catch (exception e)
//		{
//			printf("\nElasticSearch:::commit error:%s", e.what());
//		}
//	}
//	try
//	{
//		ReleaseAll__ImRawIm(vwI, myJSON, vwS, ES_id);
//	}
//	catch (exception e)
//	{
//		printf("\nElasticSearch:::release error:%s", e.what());
//	}
//}
//
//int post_process_step(uchar_descriptors *query,vector<string> dscV, vector<string> fnV, vector<float> scoresV, int numEsReturns)
//{
//	vector<float> scores;
//	std::vector<float> scoresSorted;
//	std::vector<int> scoreRank;
//	Mat descriptorQuery; 
//	query->get_descriptors(descriptorQuery);
//	std::vector<Point2f> coordsQuery = query->getCoords();
//	//std::vector<float> oriQuery, scaleQuery;
//	//TODO: runOptions Create
//	//TODO: runOptions Create
//	//int numCorr = runOptions.numMatches;
//	//scores.reserve(numCorr);
//	//scores.resize(numCorr);
//
//#if defined _OPENMP
////#pragma omp parallel for ordered schedule(dynamic)
//#endif
//	for (int i = 0; i<numEsReturns; i++)
//	{
//		//int threadId = omp_get_thread_num();
//		// read the signature of the descriptor
//		if (IS_DscFile(dscV[i].c_str()))
//		{
//			uchar_descriptors matchDsc(dscV[i].c_str(), AKAZE_FEATS);
//			matchDsc.read_dsc();
//			Mat descriptorMatch;
//			matchDsc.get_descriptors(descriptorMatch);
//			vector<Point2f> coordsMatch = matchDsc.getCoords();
//
//			// Match with FLANN
//			std::vector<DMatch > matches, good_matches;
//			cv_FLANN_Matcher(descriptorQuery, descriptorMatch, matches, good_matches);
//			// TODO : compute H using estimate_homography - shown below
//			// DONE : compute similarity
//			double score;
//			cv_GeoRR_Scoring_Location(coordsQuery, coordsMatch, good_matches, score, T_SCORE_FGC_WEIGHTED);
//#if defined _OPENMP
//#pragma omp ordered
//#endif
//			{
//				scores.push_back(score);
//			}
//
//			matches.clear();
//			coordsMatch.clear();
//			good_matches.clear();
//		}
//		else
//		{
//#pragma omp ordered
//			{
//				scores.push_back(-1);
//			}
//			printf("Unable to read file %s. Skipping...\n", dscV[i].c_str());
//		}
//		// DONE : fill in the score
//		// DONE : check that score sorting is after this point
//	}
//
//	//indexed_sort(scores, scoresSorted, scoreRank, 1);
//
//	//std::vector<TMatch> sortedMatchList;
//	//sortedMatchList.resize(numCorr);
//	//for (int i = 0; i<numCorr; i++)
//	//{
//	//	int newInd = scoreRank[i];
//	//	sortedMatchList[i].score = scoresSorted[i];
//	//	sortedMatchList[i].fileName = runOptions.matchList[newInd].fileName;
//	//	sortedMatchList[i].label = runOptions.matchList[newInd].label;
//	//}
//
//
//	descriptorQuery.release();
//	coordsQuery.clear();
//	scoresSorted.clear();
//	scoreRank.clear();
//	scores.clear();
//
//	return 0;
//}
