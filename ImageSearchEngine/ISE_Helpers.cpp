#include "ISE_Helpers.h"


void ReleaseAll__ImRawIm(unsigned* vwI, json_t* myJSON, std::string& vwS, const char* ES_id)
{
	delete[] vwI;
	json_decref(myJSON);
	vwS = "";
	delete ES_id;
}

void ReleaseAll__ImRawIm(unsigned* vwI, unsigned* vwI_low, json_t* myJSON, std::string& vwS, std::string& vwS_low, const char* ES_id)
{
	delete[] vwI;
	delete[] vwI_low;
	json_decref(myJSON);
	vwS = "";
	vwS_low = "";
	delete ES_id;
}

void WriteCSV(std::vector<std::vector<std::string>> dataVV, char* fileName, int fileNum)
{
	std::ofstream myfile;
	std::string CSV_Path = fileName;
	myfile.open(CSV_Path.c_str());
	for (int i = 0; i < dataVV.size(); i++){
		std::string lineStr = "";
		for (int ii = 0; ii < fileNum; ii++){
			lineStr += dataVV[i][ii] + ";";
		}
		lineStr += "\n";
		myfile << lineStr;
	}
}

void WriteCSV(std::vector<std::vector<float>> dataVV, char* fileName, int fileNum)
{
	std::ofstream myfile;
	std::string CSV_Path = fileName;
	myfile.open(CSV_Path.c_str());
	for (int i = 0; i < dataVV.size(); i++){
		std::string lineStr = "";
		for (int ii = 0; ii < fileNum; ii++){
			lineStr += std::to_string(dataVV[i][ii]) + ";";
		}
		lineStr += "\n";
		myfile << lineStr;
	}
}

void paramsConfig(Path &myPath, ELK_params &myES)
{
	myPath.dscFoldName = "dsc_akaze2";
	myPath.dscFoldName2 = "dsc_akaze_low";
	myPath.DataSet = "Z:/2016";
	myPath.imgFoldName = "images";
	myPath.subFolderingLevel = 2;
	myPath.VocTree = "D:/v2/voctree";
	//myPath.VocTreeLow = "D:\v2\voctree";

	myES.index = "akaze_test";
	myES.type = "2016";
	//myES.url = "http://172.16.10.202:9200"; //ImageServer2 from local network
	//myES.url = "http://85.105.103.135:9202"; //ImageServer2 from different network
	myES.url = "http://10.254.101.171:3000"; //AA 
	myES.userPWD = "";
}

void ImageConfig(std::vector<std::string> vList, int m, std::string imagePath, Image_Info& myIm, bool imgPath)
{
	myIm.dataSet = "2016";
	myIm.dataSubSet = "";
	myIm.descriptorType = "aa test";
	myIm.encoding = "jpg";
	myIm.fileName = (imgPath) ? vList[m] : vList[m].substr(0, vList[m].length() - 4);
	//myIm.fileName = vList[m];
	//myIm.fileName = vList[m].substr(0, vList[m].length() - 4);
	myIm.height = 0.0;
	myIm.width = 0.0;
	myIm.Import = "true";
	myIm.Query = "false";
	myIm.path = imagePath;
	myIm.source_type = "AA test";
}

void ImportRawImage(Path myPath, ELK_params myES, TVoctreeVLFeat* VT, std::string dscPath,
	Image_Info myIm, uchar_descriptors * my_desc)
{
	unsigned int * vwI = new unsigned int[my_desc->GetNumOfDescriptors()];
	json_t* myJSON = json_object();
	std::string vwS = "";
	const char * ES_id = new    char;

	if (my_desc->GetNumOfDescriptors() > 0)
	{
		try
		{
			VT->quantize_multi(vwI, my_desc->GetUCHAR_descriptors(), my_desc->GetNumOfDescriptors(), 61);
			for (unsigned int s = 0; s < my_desc->GetNumOfDescriptors(); s++)
				vwS += " " + int2string(int(vwI[s]));

			if (vwS != "")
			{
				GetJSON__NewImage(myIm, myPath, myJSON, vwS);
				ELK__Commit(myES, myJSON, ES_id, myIm.fileName.c_str());
			}
		}
		catch (std::exception e)
		{
			printf("\nElasticSearch:::commit error:%s", e.what());
		}
	}
	try
	{
		ReleaseAll__ImRawIm(vwI, myJSON, vwS, ES_id);
	}
	catch (std::exception e)
	{
		printf("\nElasticSearch:::release error:%s", e.what());
	}
}

void ImportRawImage(Path myPath, ELK_params myES, TVoctreeVLFeat* VT, std::string dscPath, Image_Info myIm,
	uchar_descriptors * my_desc, TVoctreeVLFeat* VT_low, std::string dsc2Path, uchar_descriptors * my_desc2)
{
	unsigned int * vwI = new unsigned int[my_desc->GetNumOfDescriptors()];
	unsigned int * vwI_low = new unsigned int[my_desc2->GetNumOfDescriptors()];
	json_t* myJSON = json_object();
	std::string vwS = "";
	std::string vwS_low = "";
	const char * ES_id = new char;

	if (my_desc->GetNumOfDescriptors() > 0 && my_desc2->GetNumOfDescriptors() > 0)
	{
		try
		{
			VT->quantize_multi(vwI, my_desc->GetUCHAR_descriptors(), my_desc->GetNumOfDescriptors(), 61);
			VT_low->quantize_multi(vwI_low, my_desc2->GetUCHAR_descriptors(), my_desc2->GetNumOfDescriptors(), 61);
			for (unsigned int s = 0; s < my_desc->GetNumOfDescriptors(); s++)
				vwS += " " + int2string(int(vwI[s]));

			for (unsigned int s = 0; s < my_desc2->GetNumOfDescriptors(); s++)
				vwS_low += " " + int2string(int(vwI_low[s]));

			if (vwS != "")
			{
				GetJSON__NewImage(myIm, myPath, myJSON, vwS, vwS_low);
				ELK__Commit(myES, myJSON, ES_id, myIm.fileName.c_str());
			}
		}
		catch (std::exception e)
		{
			printf("\nElasticSearch:::commit error:%s", e.what());
		}
	}
	try
	{
		ReleaseAll__ImRawIm(vwI, vwI_low, myJSON, vwS, vwS_low, ES_id);
	}
	catch (std::exception e)
	{
		printf("\nElasticSearch:::release error:%s", e.what());
	}
}

void QueryRawImage(Path myPath, ELK_params myES, TVoctreeVLFeat* VT, std::string dscPath, Image_Info myIm,
	uchar_descriptors *my_desc, std::vector<std::string> &testSet, std::vector<float> &scoresPP, std::vector<float> &scoresELK)
{
	unsigned int * vwI = new unsigned int[my_desc->GetNumOfDescriptors()];
	json_t* myJSON = json_object();
	std::string vwS = "";
	const char * ES_id = new char;
	std::vector<std::string> dscPathsV;
	int totalNumELK;

	if (my_desc->GetNumOfDescriptors() > 0)
	{
		try
		{
			VT->quantize_multi(vwI, my_desc->GetUCHAR_descriptors(), my_desc->GetNumOfDescriptors(), 61);
			for (unsigned int s = 0; s < my_desc->GetNumOfDescriptors(); s++)
				vwS += " " + int2string(int(vwI[s]));

			if (vwS != "")
			{
				printf("nok \n");
				GetJSON__QueryImage(myJSON, vwS, "words_string");
				//TODO: add scores at ELK_PostQuery
				ELK_PostQuery(myES, myJSON, myIm, testSet, dscPathsV, scoresELK, totalNumELK);
				
				//ELK score filtering//////////////////////////
				double sum = 0;
				int forLimit = 0;
				std::vector<float> normELKScore;
				forLimit = totalNumELK < 10 ? totalNumELK : 10;
				for (int i = 0; i < forLimit; i++)
					sum += scoresELK[i];
				for (int i = 0; i < forLimit; i++)
				{
					float scoreNorm = scoresELK[i] / sum;
					if (scoreNorm > 0.08)
						normELKScore.push_back(scoreNorm);
				}
				////////////////////////////////////////////////
				postProcess(my_desc, dscPathsV, testSet, scoresPP, normELKScore.size());

				printf("ok \n");
			}
		}
		catch (std::exception e)
		{
			printf("\nElasticSearch:::commit error:%s", e.what());
		}
	}
	try
	{
		dscPathsV.clear();
		dscPathsV.shrink_to_fit();
		ReleaseAll__ImRawIm(vwI, myJSON, vwS, ES_id);
	}
	catch (std::exception e)
	{
		printf("\nElasticSearch:::release error:%s", e.what());
	}
}

int postProcess(uchar_descriptors *query, std::vector<std::string> dscV, std::vector<std::string> fnV,
                std::vector<float> &scores, int numELK)
{
	std::vector<float> scoresSorted;
	std::vector<int> scoreRank;
	Mat descriptorQuery;
	query->Get_CVDescriptors(descriptorQuery);
	std::vector<Point2f> coordsQuery = query->GetCoords();
	//std::vector<float> oriQuery, scaleQuery;
	//TODO: runOptions Create
	//TODO: runOptions Create
	//int numCorr = runOptions.numMatches;
	//scores.reserve(numCorr);
	//scores.resize(numCorr);

#if defined _OPENMP
	#pragma omp parallel for ordered schedule(dynamic)
#endif
	for (int i = 0; i<numELK; i++)
	{
		//int threadId = omp_get_thread_num();
		// read the signature of the descriptor
		if (IS_DscFile(dscV[i].c_str()) || IS_ImageFile(dscV[i].c_str()))
		{
			uchar_descriptors matchDsc(dscV[i].c_str(), dscV[i].c_str(), AKAZE_FEATS);
			Mat descriptorMatch;
			if (IS_DscFile(dscV[i].c_str()))
			{
				matchDsc.ReadDSC();
				descriptorMatch = Mat(matchDsc.GetNumOfDescriptors(), matchDsc.GetFeatureSize(), CV_8UC1);;
				matchDsc.GetReadModeDescriptors(descriptorMatch);
			}
			else if (IS_ImageFile(dscV[i].c_str()))
			{
				matchDsc.ExtractAKAZE();
				matchDsc.Get_CVDescriptors(descriptorMatch);
			}

			std::vector<Point2f> coordsMatch = matchDsc.GetCoords();

			// Match with FLANN
			std::vector<DMatch > matches, good_matches;
			if (descriptorQuery.type() == 0)  descriptorQuery.convertTo(descriptorQuery, CV_32F);
			if (descriptorMatch.type() == 0)  descriptorMatch.convertTo(descriptorMatch, CV_32F);
			
			if (descriptorQuery.type() == descriptorMatch.type() && descriptorQuery.cols == descriptorMatch.cols)
				cv_FLANN_Matcher(descriptorQuery, descriptorMatch, matches, good_matches);
			// TODO : compute H using estimate_homography - shown below
			// DONE : compute similarity
			double score;
			cv_GeoRR_Scoring_Location(coordsQuery, coordsMatch, good_matches, score, T_SCORE_FGC_WEIGHTED);
#if defined _OPENMP
#pragma omp ordered
#endif
			{
				scores.push_back(score);
			}

			matches.clear();
			coordsMatch.clear();
			good_matches.clear();
		}
		else
		{
#pragma omp ordered
			{
				scores.push_back(-1);
			}
			printf("Unable to read file %s. Skipping...\n", dscV[i].c_str());
		}
		// DONE : fill in the score
		// DONE : check that score sorting is after this point
	}

	//indexed_sort(scores, scoresSorted, scoreRank, 1);

	//std::vector<TMatch> sortedMatchList;
	//sortedMatchList.resize(numCorr);
	//for (int i = 0; i<numCorr; i++)
	//{
	//	int newInd = scoreRank[i];
	//	sortedMatchList[i].score = scoresSorted[i];
	//	sortedMatchList[i].fileName = runOptions.matchList[newInd].fileName;
	//	sortedMatchList[i].label = runOptions.matchList[newInd].label;
	//}
	descriptorQuery.release();
	coordsQuery.clear();
	scoresSorted.clear();
	scoreRank.clear();

	return 0;
}

void ImageSpliter(const cv::Mat Input, std::vector<cv::Mat>& OutputVector, std::vector<std::string>& OutputNames, int maxSize)
{
	int height = Input.rows;
	int widht = Input.cols;

	int hNum = (height / (maxSize*0.7)) + 1;
	int wNum = (widht / (maxSize*0.7)) + 1;

	int hStart = 0;
	int hStep = maxSize;

	for (int h_i = 0; h_i < hNum ; h_i++)
	{
		if ((hStart + hStep) > height)
			break;
		int wStart = 0;
		int wStep = maxSize;
		for (int w_i = 0; w_i < wNum ; w_i++)
		{
			if ((wStart + wStep) > widht)
				break;
			cv::Rect myROI(wStart, hStart, wStep, hStep);
			cv::Mat croppedImage = Input(myROI);

			OutputVector.push_back(croppedImage);
			OutputNames.push_back(int2string(h_i) +"_"+int2string(w_i));

			wStart += (maxSize*0.7);
			if (w_i == wNum - 2)
				wStep = widht - wStart;
		}
		hStart += (maxSize*0.7);
		if (h_i == hNum - 2)
			hStep = height - hStart;
	}


}

