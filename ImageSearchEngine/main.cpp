#include "ISE_Helpers.h"

int main(int argc, char *argv[])
{
	Path myPath;
	ELK_params myES;
	TVoctreeVLFeat* VT = new TVoctreeVLFeat;
	TVoctreeVLFeat* VT_low = new TVoctreeVLFeat;
	std::vector<std::vector<std::string>> allEstimates, allFN;
	std::vector<std::vector<float>> allScorePP, allScoreELK, allScoreW;
	long start = clock();
	double duration;
	unsigned int err_count = 0;

	paramsConfig(myPath, myES);
	VocTreeInit(VT, myPath);

	//Face Reco Import New 
	/*
	GET_FolderList(myPath.DataSet.c_str(), foldList);
	for (unsigned int f = 0; f < foldList.size(); f++)
	{
	vector <std::string> imgList;
	String imgFoldPath;
	imgFoldPath = myPath.DataSet + "/" + foldList[f];
	GET_DirectoryImages(imgFoldPath.c_str(), imgList);

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(8); // Use 4 threads for all consecutive parallel regions
	#pragma omp parallel for
	for (int i = 0; i < imgList.size(); i++)
	{
	string imgPath = imgFoldPath + "/" + imgList[i];
	Image_Info myIm;
	ImageConfig(imgList, i, imgPath.c_str(), myIm, true);

	if (IS_ImageFile(imgPath.c_str()))
	{
	try
	{
	uchar_descriptors my_desc(imgPath.c_str(), "", AKAZE_FEATS);
	my_desc.ExtractAKAZE();
	printf("nok\n");
	ImportRawImage(myPath, myES, VT, imgPath, myIm, &my_desc);
	printf("ok\n");
	}
	catch (exception e)
	{
	printf("\nMain:::extract akaze and write:%s", e.what());
	}
	}
	}
	}//*/
	//Face Reco Query
	/*
	GET_FolderList(myPath.DataSet.c_str(), foldList);
	//for (unsigned int f = 0; f < foldList.size(); f++)
	for (unsigned int f = 0; f < 1; f++)
	{
	std::vector <std::string> imgList;
	cv::String imgFoldPath;
	imgFoldPath = myPath.DataSet + "/" + foldList[f];
	GET_DirectoryImages(imgFoldPath.c_str(), imgList);

	std::vector<std::string> estimatesFold, fnFold;
	std::vector<float> scoresPP__F, scoresELK__F, scoresW__F;
	estimatesFold.push_back(foldList[f] + ".jpg");
	scoresELK__F.push_back(atof(foldList[f].c_str()));
	scoresPP__F.push_back(atof(foldList[f].c_str()));
	scoresW__F.push_back(atof(foldList[f].c_str()));
	for (int i = 0; i < 10 * 25; i++)
	{
	estimatesFold.push_back("-1");
	scoresELK__F.push_back(-1);
	scoresPP__F.push_back(-1);
	scoresW__F.push_back(-1);
	}
	for (int i = 0; i < 25; i++)
	fnFold.push_back("NoName");

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(8); // Use 4 threads for all consecutive parallel regions
	#pragma omp parallel for
	for (int m = 0; m < imgList.size(); m++)
	{
	std::string imgPath = imgFoldPath + "/" + imgList[m];
	Image_Info myIm;
	ImageConfig(imgList, m, imgPath.c_str(), myIm, true);
	fnFold[m] = imgList[m];

	if (IS_ImageFile(imgPath.c_str()))
	{
	try
	{
	std::vector<std::string> testIm;
	std::vector<float> scoresPP, scoresELK, scoresW;
	uchar_descriptors my_desc(imgPath.c_str(), "", AKAZE_FEATS);
	my_desc.ExtractAKAZE();
	QueryRawImage(myPath, myES, VT, "NoDSC", myIm, &my_desc, testIm,
	scoresPP, scoresELK);
	//score weighting    w1*ELKscore + w2*PPscore
	for (unsigned int l = 0; l < scoresPP.size(); l++)
	{
	float wScore = (0.3 * scoresELK[l]) + (0.7 * scoresPP[l]);
	scoresW.push_back(wScore);
	}
	printf("%d\n", m);
	if (testIm.size() == 10)
	{
	for (int i = 0; i < scoresW.size(); i++)
	{
	estimatesFold[(m * 10) + i + 1] = testIm[i];
	scoresELK__F[(m * 10) + i + 1] = scoresELK[i];
	scoresPP__F[(m * 10) + i + 1] = scoresPP[i];
	scoresW__F[(m * 10) + i + 1] = scoresW[i];
	}
	testIm.clear(); testIm.shrink_to_fit();
	scoresPP.clear(); scoresPP.shrink_to_fit();
	scoresELK.clear(); scoresELK.shrink_to_fit();
	scoresW.clear(); scoresW.shrink_to_fit();
	}
	else
	{
	for (int i = 0; i < 11; i++)
	{
	estimatesFold.push_back("-1");
	scoresELK__F.push_back(-1);
	scoresPP__F.push_back(-1);
	scoresW__F.push_back(-1);
	}
	}
	printf("ok\n");
	}
	catch (std::exception e)
	{
	printf("\nMain:::extract akaze and write:%s", e.what());
	}
	}
	}
	allEstimates.push_back(estimatesFold);
	estimatesFold.clear(); estimatesFold.shrink_to_fit();
	allFN.push_back(fnFold);
	fnFold.clear(); fnFold.shrink_to_fit();
	//Write from folder to dataset
	allScoreELK.push_back(scoresELK__F);
	allScorePP.push_back(scoresPP__F);
	allScoreW.push_back(scoresW__F);
	scoresELK__F.clear(); scoresELK__F.shrink_to_fit();
	scoresPP__F.clear(); scoresPP__F.shrink_to_fit();
	scoresW__F.clear(); scoresW__F.shrink_to_fit();
	printf("\rFolder No:  = %d\n", f);
	}
	foldList.clear();
	foldList.shrink_to_fit();
	WriteCSV(allEstimates, "FaceSub01.csv", 251);
	WriteCSV(allScoreELK, "FaceSub01_scoreELK.csv", 251);
	WriteCSV(allScorePP, "FaceSub01_scorePP.csv", 251);
	WriteCSV(allScoreW, "FaceSub01_scoreW.csv", 251);
	WriteCSV(allFN, "FaceSub01_fn.csv", 25);
	//*/
	//Image Search Import New 
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
	my_desc.ReadDSC();
	//my_desc.ExtractAKAZE();

	//uchar_descriptors my_desc2(imgPath.c_str(), dsc2Path.c_str(), AKAZE_FEATS);
	//my_desc2.ExtractAKAZE_low();

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
	//Image search query
	/*
	//for (unsigned int f = 0; f < foldList.size(); f++)
	for (unsigned int f = 0; f < 50; f++)
	{
	String imgFoldPath;
	imgFoldPath = myPath.DataSet + "/" + foldList[f]  ;
	GET_DirectoryImages(imgFoldPath.c_str(), imgList);
	vector<string> testSetFold;
	vector<float> scoresPP__F, scoresELK__F, scoresW__F;
	testSetFold.push_back(foldList[f] + ".jpg");
	scoresELK__F.push_back(atof(foldList[f].c_str()));
	scoresPP__F.push_back(atof(foldList[f].c_str()));
	scoresW__F.push_back(atof(foldList[f].c_str()));
	for (int i = 0; i < 10 * 25; i++)
	{
	testSetFold.push_back("-1");
	scoresELK__F.push_back(-1);
	scoresPP__F.push_back(-1);
	scoresW__F.push_back(-1);
	}

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(8); // Use 4 threads for all consecutive parallel regions
	#pragma omp parallel for
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
	vector<float> scoresPP, scoresELK, scoresW;

	uchar_descriptors my_desc(imgPath.c_str(), "NoDSC", AKAZE_FEATS);
	my_desc.ExtractAKAZE();

	QueryRawImage(myPath, myES, VT, "NoDSC", myIm, &my_desc, testIm,
	scoresPP, scoresELK);
	//score weighting    w1*ELKscore + w2*PPscore
	for (unsigned int l = 0; l < scoresPP.size(); l++)
	{
	float wScore = (0.3 * scoresELK[l]) + (0.7 * scoresPP[l]);
	scoresW.push_back(wScore);
	}

	printf("%d\n", m);

	if (testIm.size() == 10)
	{
	for (int i = 0; i < scoresW.size(); i++)
	{
	testSetFold[(m * 10) + i + 1]	=	testIm[i];
	scoresELK__F[(m * 10) + i + 1]	=	scoresELK[i];
	scoresPP__F[(m * 10) + i + 1] = scoresPP[i];
	scoresW__F[(m * 10) + i + 1] = scoresW[i];
	}

	testIm.clear();testIm.shrink_to_fit();
	scoresPP.clear();scoresPP.shrink_to_fit();
	scoresELK.clear(); scoresELK.shrink_to_fit();
	scoresW.clear(); scoresW.shrink_to_fit();
	}
	else
	{
	for (int i = 0; i < 11; i++)
	{
	testSetFold.push_back("-1");
	scoresELK__F.push_back(-1);
	scoresPP__F.push_back(-1);
	scoresW__F.push_back(-1);
	}
	}
	}
	catch (exception e)
	{
	printf("\nMain:::extract akaze and write:%s", e.what());
	}
	}
	}

	testSet.push_back(testSetFold);
	testSetFold.clear();testSetFold.shrink_to_fit();
	//Write from folder to dataset
	scoreELK_set.push_back(scoresELK__F);
	scorePP_set.push_back(scoresPP__F);
	scoreW_set.push_back(scoresW__F);
	scoresELK__F.clear(); scoresELK__F.shrink_to_fit();
	scoresPP__F.clear(); scoresPP__F.shrink_to_fit();
	scoresW__F.clear(); scoresW__F.shrink_to_fit();
	printf("\rFolder No:  = %d\n", f);

	imgList.clear();imgList.shrink_to_fit();
	if (f % 5 == 0)
	{
	double num = (f*100.00 / 1000);
	printf("\rProcess Rate = %.2f\n", num);
	}
	}
	foldList.clear();
	foldList.shrink_to_fit();
	WriteCSV(testSet, "TESTIM33.csv", 251);
	WriteCSV(scoreELK_set, "TESTIM33_scoreELK.csv", 251);
	WriteCSV(scorePP_set, "TESTIM33_scorePP.csv", 251);
	WriteCSV(scoreW_set, "TESTIM33_scoreW.csv", 251);
	//*/
	//Istanbul Import New
	/*
	std::vector<std::string> imgList;
	GET_DirectoryImages(myPath.DataSet.c_str(), imgList);
	for (int i = 0; i < imgList.size(); i++)
	{
		std::string imgPath = myPath.DataSet + "/" + imgList[i];
		Image_Info myIm;
		ImageConfig(imgList, i, imgPath.c_str(), myIm, true);
		if (IS_ImageFile(imgPath.c_str()))
		{
			try
			{
				uchar_descriptors my_desc(imgPath.c_str(), "", AKAZE_FEATS);
				my_desc.setResizeImage(true);
				my_desc.ExtractAKAZE();
				myIm.height = my_desc.GetImageHeight();
				myIm.width = my_desc.GetImageWidth();
				printf("nok\n");
				ImportRawImage(myPath, myES, VT, imgPath, myIm, &my_desc);
				printf("ok\n");
			}
			catch (std::exception e)
			{
				printf("\nMain:::extract akaze and write:%s", e.what());
			}
		}
	}
	imgList.clear();
	imgList.shrink_to_fit();
	//*/
	//Istanbul Query for large sizes
	/*
	std::vector<std::string> imgList;
	GET_DirectoryImages(myPath.DataSet.c_str(), imgList);
	for (int i = 0; i < imgList.size(); i++)
	{
		printf("file: %d\n", i);
		std::string imgPath = myPath.DataSet + "/" + imgList[i];
		Image_Info myIm;
		ImageConfig(imgList, i, imgPath.c_str(), myIm, true);
		if (IS_ImageFile(imgPath.c_str()))
		{
			try
			{
				cv::Mat Image = cv::imread(imgPath, cv::IMREAD_GRAYSCALE);
				std::vector<cv::Mat> imagesVector;
				std::vector<std::string> imagesName;
				if (Image.rows > 800 || Image.cols > 800)
					ImageSpliter(Image, imagesVector, imagesName, 800);
				
				uchar_descriptors my_desc(imgPath.c_str(), Image, "", AKAZE_FEATS);
				my_desc.setResizeImage(true);
				my_desc.ExtractAKAZE();
				myIm.height = my_desc.GetImageHeight();
				myIm.width = my_desc.GetImageWidth();

				for (int f = 0; f < imagesVector.size(); f++)
				{
					printf("nok\n");
					uchar_descriptors my_desc_f(imagesName[f].c_str(), imagesVector[f], "", AKAZE_FEATS);
					//my_desc.setResizeImage(true);
					my_desc_f.ExtractAKAZE();
					myIm.height = my_desc_f.GetImageHeight();
					myIm.width = my_desc_f.GetImageWidth();
					printf("ok\n");
				}

			}
			catch (std::exception e)
			{
				printf("\nMain:::extract akaze and write:%s", e.what());
			}
		}
	}
	imgList.clear();
	imgList.shrink_to_fit();
	//*/
	// AA Import Images
	/*
	if (argc >= 1)
		myPath.DataSet = argv[1];
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	int numThread = omp_get_max_threads();
	omp_set_num_threads(numThread); // Use 4 threads for all consecutive parallel regions
	std::vector<std::string> monthList;
	GET_FolderList(myPath.DataSet.c_str(), monthList);
	for (size_t m = 0; m < monthList.size(); m++)
	{

		auto mPath = myPath.DataSet + "/" + monthList[m];
		std::vector<std::string> dayList;
		GET_FolderList(mPath.c_str(), dayList);
		for (size_t d = 0; d < dayList.size(); d++)
		{
			auto dPath = mPath + "/" + dayList[d];
			std::vector<std::string> imgList;
			GET_DirectoryImages(dPath.c_str(), imgList);

#pragma omp parallel for
			for (int i = 0; i < imgList.size(); i++)
			{
				auto imgPath = dPath + "/" + imgList[i];
				Image_Info myIm;
				ImageConfig(imgList, i, imgPath.c_str(), myIm, true);

				if (IS_ImageFile(imgPath.c_str()))
				{
					try
					{
						printf("nok\n");
						uchar_descriptors myDesc(imgPath.c_str(), "", AKAZE_FEATS);
						myDesc.setResizeImage(true);
						myDesc.ExtractAKAZE();
						myIm.height = myDesc.GetImageHeight();
						myIm.width = myDesc.GetImageWidth();
						ImportRawImage(myPath, myES, VT, imgPath, myIm, &myDesc);
						printf("ok\n");
					}
					catch (std::exception e)
					{
						printf("\nMain:::extract akaze and write:%s", e.what());
					}
				}
			}
		}
	}
	//*/
	// AA Query Images
	///*
	if (argc >= 1 && argv[1]!= nullptr)
		myPath.DataSet = argv[1];
	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	int numThread = omp_get_max_threads();
	omp_set_num_threads(numThread); // Use 4 threads for all consecutive parallel regions
	std::vector<std::string> imgList;
	std::vector<std::vector<std::string>> outputsV;
	GET_DirectoryImages(myPath.DataSet.c_str(), imgList);

//#pragma omp parallel for
	for (int i = 0; i < imgList.size(); i++)
	//for (int i = 0; i <10; i++)
	{
		auto imgPath = myPath.DataSet + "/" + imgList[i];
		Image_Info myIm;
		ImageConfig(imgList, i, imgPath.c_str(), myIm, true);

		if (IS_ImageFile(imgPath.c_str()))
		{
			try
			{
				cv::Mat Image = cv::imread(imgPath, cv::IMREAD_GRAYSCALE);
				std::vector<cv::Mat> imagesVector;
				std::vector<std::string> imageNameList;
				// variables for scoring
				std::vector<std::string> testIm;
				std::vector<float> scoresPP, scoresELK, scoresW;
				if (Image.rows > 800 || Image.cols > 800)
					ImageSpliter(Image, imagesVector, imageNameList, 800);
				else
				{
					imagesVector.push_back(Image);
					imageNameList.push_back("raw");
				}
				//#pragma omp parallel for
				for (int f = 0; f < imagesVector.size(); f++)
				{
					printf("nok\n");
					std::vector<std::string> outputLine;
					// 1. Feature Extraction
					uchar_descriptors my_desc_f(imageNameList[f].c_str(), imagesVector[f], "", AKAZE_FEATS);
					my_desc_f.ExtractAKAZE();
					myIm.height = my_desc_f.GetImageHeight();
					myIm.width = my_desc_f.GetImageWidth();
					// 2. ELK Query and Response
					QueryRawImage(myPath, myES, VT, "NoDSC", myIm, &my_desc_f, testIm, scoresPP, scoresELK);
					// 3. Score Weighting   ( w1*ELKscore + w2*PPscore )
					if (testIm.size() > 0)
						scoreWeighting(scoresPP, scoresELK, scoresW);
					printf("%d::%d\n", i, f);
					// 4. Score
					outputLine.push_back(imgList[i].c_str());
					outputLine.push_back(imageNameList[f].c_str());
					outputLine.push_back(int2string(int(my_desc_f.GetNumOfDescriptors())));
					
					for (int p = 0; p < QUERY_RETURN_SIZE; p++)
					{
						if (p < scoresW.size()){
							outputLine.push_back(float2string(scoresPP[p]));
							outputLine.push_back(float2string(scoresELK[p]));
							outputLine.push_back(float2string(scoresW[p]));
							outputLine.push_back(testIm[p].c_str());
						}
						else{ for (int n = 0; n < 4; n++)outputLine.push_back("NoScore"); }
					}
					
					outputsV.push_back(outputLine);
					outputLine.clear(); outputLine.shrink_to_fit();
					// 5. Releases
					testIm.clear(); testIm.shrink_to_fit();
					scoresPP.clear(); scoresPP.shrink_to_fit();
					scoresELK.clear(); scoresELK.shrink_to_fit();
					scoresW.clear(); scoresW.shrink_to_fit();
					printf("ok\n");
				}
			}
			catch (std::exception e)
			{
				printf("\nMain:::extract akaze and write:%s", e.what());
			}
		}
		else 
			printf("\nNot an Image File: %s", imgPath.c_str());
	}
	std::string wscPath = "Test2015.csv";
	if (argc > 1) wscPath = argv[2];
	WriteCSV(outputsV, wscPath, 3 + (4*QUERY_RETURN_SIZE));

//*/
	duration = clock() - start;
	delete VT;
	//delete VT_low;
	printf(":::Duration: %.2f, Error count: %d", duration, err_count);
	return 0;
}

//ISE HELPERS
////#include "helpers2.h"
//#include "postProc.h"
//#include <vector>
//#include "descriptors.h"
//#include <omp.h>
//#include "ES_image.h"
//#include "TVoctreeVLFeat.h"
//#include <array>
//#include "vld.h"
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
//void paramsConfig(Path &myPath, ELK_params &myES)
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
//void ImportRawImage(Path myPath, ELK_params myES, TVoctreeVLFeat* VT, string dscPath, 
//	Image_Info myIm, uchar_descriptors * my_desc)
//{
//	unsigned int * vwI = new unsigned int[my_desc->GetNumOfDescriptors()];
//	json_t* myJSON = json_object();
//	string vwS = "";
//	const char * ES_id = new    char;
//
//	if (my_desc->GetNumOfDescriptors() > 0)
//	{
//		try
//		{
//			VT->quantize_multi(vwI, my_desc->GetUCHAR_descriptors(), my_desc->GetNumOfDescriptors(), 61);
//			for (unsigned int s = 0; s < my_desc->GetNumOfDescriptors(); s++)
//				vwS += " " + int2string(int(vwI[s]));
//
//			if (vwS != "")
//			{
//				GetJSON__NewImage(myIm, myPath, myJSON, vwS);
//				ELK__Commit(myES, myJSON, ES_id, myIm.fileName.c_str());
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
//void ImportRawImage(Path myPath, ELK_params myES, TVoctreeVLFeat* VT, string dscPath, Image_Info myIm, 
//	uchar_descriptors * my_desc, TVoctreeVLFeat* VT_low, string dsc2Path, uchar_descriptors * my_desc2)
//{
//	unsigned int * vwI = new unsigned int[my_desc->GetNumOfDescriptors()];
//	unsigned int * vwI_low = new unsigned int[my_desc2->GetNumOfDescriptors()];
//	json_t* myJSON = json_object();
//	string vwS = "";
//	string vwS_low = "";
//	const char * ES_id = new char;
//
//	if (my_desc->GetNumOfDescriptors() > 0 && my_desc2->GetNumOfDescriptors() > 0)
//	{
//		try
//		{
//			VT->quantize_multi(vwI, my_desc->GetUCHAR_descriptors(), my_desc->GetNumOfDescriptors(), 61);
//			VT_low->quantize_multi(vwI_low, my_desc2->GetUCHAR_descriptors(), my_desc2->GetNumOfDescriptors(), 61);
//			for (unsigned int s = 0; s < my_desc->GetNumOfDescriptors(); s++)
//				vwS += " " + int2string(int(vwI[s]));
//
//			for (unsigned int s = 0; s < my_desc2->GetNumOfDescriptors(); s++)
//				vwS_low += " " + int2string(int(vwI_low[s]));
//
//			if (vwS != "")
//			{
//				GetJSON__NewImage(myIm, myPath, myJSON, vwS, vwS_low);
//				ELK__Commit(myES, myJSON, ES_id, myIm.fileName.c_str());
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
//void QueryRawImage(Path myPath, ELK_params myES, TVoctreeVLFeat* VT, string dscPath, Image_Info myIm,
//	uchar_descriptors * my_desc, vector<string> &testSet, string &returnFileName)
//{
//	unsigned int * vwI = new unsigned int[my_desc->GetNumOfDescriptors()];
//	json_t* myJSON = json_object();
//	string vwS = "";
//	const char * ES_id = new char;
//	vector<string> dscPathsV;
//	vector<float> scoresV;
//
//	if (my_desc->GetNumOfDescriptors() > 0 )
//	{
//		try
//		{
//			VT->quantize_multi(vwI, my_desc->GetUCHAR_descriptors(), my_desc->GetNumOfDescriptors(), 61);
//			for (unsigned int s = 0; s < my_desc->GetNumOfDescriptors(); s++)
//				vwS += " " + int2string(int(vwI[s]));
//
//			if (vwS != "")
//			{
//				printf("nok \n");
//				GetJSON__QueryImage(myJSON, vwS, "words_string");
//				//TODO: add scores at ELK_PostQuery
//				ELK_PostQuery(myES, myJSON, myIm, testSet, dscPathsV, scoresV);
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
//	query->Get_CVDescriptors(descriptorQuery);
//	std::vector<Point2f> coordsQuery = query->GetCoords();
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
//			matchDsc.ReadDSC();
//			Mat descriptorMatch;
//			matchDsc.Get_CVDescriptors(descriptorMatch);
//			vector<Point2f> coordsMatch = matchDsc.GetCoords();
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
