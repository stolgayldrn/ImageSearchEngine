#include "ISE_Helpers.h"


void ReleaseAll__ImRawIm(unsigned* vwI, json_t* myJSON, std::string& vwS, const char* ES_id)
{
	delete[] vwI;
	json_decref(myJSON);
	vwS = "";
	delete ES_id;
}

//void WriteCSV(std::vector<std::vector<std::string>> dataVV, std::string CSV_Path, int fileNum)
//{
//	std::ofstream myfile;
//	//std::string CSV_Path = fileName;
//	myfile.open(CSV_Path.c_str());
//	for (int i = 0; i < dataVV.size(); i++){
//		std::string lineStr = "";
//		for (int ii = 0; ii < fileNum; ii++){
//			lineStr += dataVV[i][ii] + ";";
//		}
//		lineStr += "\n";
//		myfile << lineStr;
//	}
//}
//
//void WriteCSV(std::vector<std::vector<float>> dataVV, char* fileName, int fileNum)
//{
//	std::ofstream myfile;
//	std::string CSV_Path = fileName;
//	myfile.open(CSV_Path.c_str());
//	for (int i = 0; i < dataVV.size(); i++){
//		std::string lineStr = "";
//		for (int ii = 0; ii < fileNum; ii++){
//			lineStr += std::to_string(dataVV[i][ii]) + ";";
//		}
//		lineStr += "\n";
//		myfile << lineStr;
//	}
//}
//
//void WriteCSV(std::vector<std::vector<int>> dataVV, char* fileName, int fileNum)
//{
//	std::ofstream myfile;
//	std::string CSV_Path = fileName;
//	myfile.open(CSV_Path.c_str());
//	for (int i = 0; i < dataVV.size(); i++){
//		std::string lineStr = "";
//		for (int ii = 0; ii < fileNum; ii++){
//			lineStr += std::to_string(dataVV[i][ii]) + ";";
//		}
//		lineStr += "\n";
//		myfile << lineStr;
//	}
//}

void paramsConfig(Path &myPath, ELK_params &myES)
{
	myPath.dscFoldName = "dsc_akaze2";
	myPath.DataSet = "C:/ImageSearch/ImageDataSet/03CustomText";
	myPath.imgFoldName = "images";
	myPath.subFolderingLevel = 2;
	myPath.corruptedFilesFolder = "D:/v2/SharedV2/corruptedFiles";// for AA server
	//myPath.corruptedFilesFolder = ""; //for local PC
	myPath.logPath = "D:/v2/SharedV2";
	
	//myPath.VocTree = "D:/v2/voctree/VT_flicker500K_AKAZE_middle_tree_S2_P"; // for AA server
	myPath.VocTree = "C:/ImageSearch/VT_Trees/VT_flicker500K_AKAZE_middle_tree_S2_P"; // PC 1
	//myPath.VocTree = "C:/ImageSearch/VT_Trees/AA_VT_Middle.dat"; // PC 2

	//myES.index = "flicker1m_test2";
	//myES.type = "istanbul2";
	//myES.url = "http://172.16.10.202:9200"; //ImageServer2 from local network
	//myES.url = "http://85.105.103.135:9202"; //ImageServer2 from different network

	myES.index = "akaze"; // AA
	myES.type = "v1"; // AA
	myES.url = "http://10.254.101.171:3000"; // AA 

	myES.userPWD = "";
}

void ImageConfig(std::vector<std::string> vList, int m, std::string imagePath, Image_Info& myIm, bool imgPath)
{
	myIm.dataSet = "2015";
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
	myIm.numDescs = 0;
	myIm.path = imagePath;
	myIm.source_type = "AA test";
}

void ImageConfig(std::string imagePath, std::string imageFileName, Image_Info& myIm, bool imgPath)
{
	myIm.dataSet = "2015";
	myIm.dataSubSet = "";
	myIm.descriptorType = "aa test";
	myIm.encoding = "jpg";
	myIm.fileName = imageFileName;
	
	myIm.height = 0.0;
	myIm.width = 0.0;
	myIm.Import = "true";
	myIm.Query = "false";
	myIm.numDescs = 0;
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
			VT->quantize_multi(vwI, my_desc->GetUCHAR_descriptors(), my_desc->GetNumOfDescriptors(), my_desc->GetFeatureSize());
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

void scoreNormELK(std::vector<float>& scoresELK, int totalNumELK, std::vector<float>& normELKScore)
{
	double sum = 0;
	int forLimit;
	forLimit = totalNumELK < QUERY_RETURN_SIZE ? totalNumELK : QUERY_RETURN_SIZE;
	for (int i = 0; i < forLimit; i++)
		sum += scoresELK[i];
	for (int i = 0; i < forLimit; i++)
	{
		float scoreNorm = scoresELK[i] / sum;
		if (scoreNorm > 0.02 && scoresELK[i] > 0.1)
			normELKScore.push_back(scoreNorm);
		else
			break;
	}
}

void QueryRawImage(Path myPath, ELK_params myES, TVoctreeVLFeat* VT, std::string dscPath, Image_Info myIm,
                   uchar_descriptors *my_desc, std::vector<std::string> & returnedFileNames, std::vector<float> &scoresPP, std::vector<float> &scoresELK)
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
			VT->quantize_multi(vwI, my_desc->GetUCHAR_descriptors(), my_desc->GetNumOfDescriptors(), my_desc->GetFeatureSize());
			for (unsigned int s = 0; s < my_desc->GetNumOfDescriptors(); s++)
				vwS += " " + int2string(int(vwI[s]));
			/*Temporary words writing*/
			/*std::vector<std::vector<int>> tempWords;
			for (int v = 0; v < my_desc->GetNumOfDescriptors(); v++)
			{
				std::vector<int>tempW;
				tempW.push_back(vwI[v]);
				tempWords.push_back(tempW);
			}
			std::string tempPath = myIm.path + ".csv";
			WriteCSV(tempWords,"a.csv", 1);*/
			/*************************/
			if (vwS != "")
			{
				GetJSON__QueryImage(myJSON, vwS, "words_string");
				
				//TODO: add scores at ELK_PostQuery
				ELK_PostQuery(myES, myJSON, myIm, returnedFileNames, dscPathsV, scoresELK, totalNumELK);
				//ELK score filtering//////////////////////////
				if (totalNumELK > 0)
				{
					std::vector<float> normELKScore;
					scoreNormELK(scoresELK, totalNumELK, normELKScore);
					// Post Process
					postProcess(my_desc, dscPathsV, returnedFileNames, scoresPP, normELKScore.size());
				}
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
	query->CopyOpencvDescriptors(descriptorQuery);
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
				matchDsc.setResizeImage(true);
				matchDsc.ExtractAKAZE();
				matchDsc.CopyOpencvDescriptors(descriptorMatch);
			}

			std::vector<Point2f> coordsMatch = matchDsc.GetCoords();

			// Match with FLANN
			std::vector<DMatch > matches, good_matches;
			if (descriptorQuery.type() == 0)  descriptorQuery.convertTo(descriptorQuery, CV_32F);
			if (descriptorMatch.type() == 0)  descriptorMatch.convertTo(descriptorMatch, CV_32F);
			float gm_dist = (query->GetFeatureType() == AKAZE_FEATS) ? GOOD_MATCHES_DISTANCE__AKAZE : GOOD_MATCHES_DISTANCE__EZSIFT;
			if (descriptorQuery.type() == descriptorMatch.type() && descriptorQuery.cols == descriptorMatch.cols)
				cv_FLANN_Matcher(descriptorQuery, descriptorMatch, matches, good_matches, gm_dist);
			// TODO : compute H using estimate_homography - shown below
			// DONE : compute similarity
			double score = 0 ;
			if (good_matches.size()> MIN_GOOD_MATCHES)
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
	if (Input.rows > maxSize || Input.cols > maxSize)
	{
		int h = Input.rows;
		int w = Input.cols;
		double hs = h*1.0 / maxSize;
		double ws = w*1.0 / maxSize;
		cv::Mat reszImg;
		if (h>w)
			cv::resize(Input, reszImg, cv::Size(w / hs, h / hs), 0, 0, CV_INTER_AREA);
		else
			cv::resize(Input, reszImg, cv::Size(w / ws, h / ws), 0, 0, CV_INTER_AREA);

		OutputVector.push_back(reszImg);
		OutputNames.push_back("resized");
	}


	/*
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
	}*/


}

void scoreWeighting(std::vector<float> scoresPP, std::vector<float> scoresELK, std::vector<float> & scoresW)
{
	for (unsigned int l = 0; l < scoresPP.size(); l++)
	{
		float wScore = (0.3 * scoresELK[l]) + (0.7 * scoresPP[l]);
		scoresW.push_back(wScore);
	}
}

std::string currentDateTime()
{
	std::chrono::time_point<std::chrono::system_clock>  now;
	now = std::chrono::system_clock::now();
	std::time_t end_time = std::chrono::system_clock::to_time_t(now);
	
	char cr[26];
	ctime_s(cr, sizeof cr, &end_time);
	cr[24] = '\0';
	std::string str(cr);
	return str;
}

void printErrorToLog(Path myPath, std::ofstream &ofs, std::vector<std::string> imgList, int i, std::basic_string<char> imgPathDUMP, uchar_descriptors myDesc, std::string errCmd)
{
	try
	{
		ofs << currentDateTime() << " : " << imgPathDUMP << " : " << errCmd <<"\n";
	}
	catch (std::exception e)
	{
		printf("\nMain:::extract akaze and write::log writer error::%s", e.what());
	}
}

void findMatches(uchar_descriptors &descriptor_1, uchar_descriptors &descriptor_2, std::vector<DMatch >& good_matches)
{
	Mat feats_1, feats_2;
	feats_1 = descriptor_1.GetOpencvDescriptors();
	feats_2 = descriptor_2.GetOpencvDescriptors();

	// Match with FLANN
	std::vector<DMatch > matches;
	if (feats_1.type() == 0)  feats_1.convertTo(feats_1, CV_32F);
	if (feats_2.type() == 0)  feats_2.convertTo(feats_2, CV_32F);

	float gm_dist = 400;

	if (feats_1.type() == feats_2.type() && feats_1.cols == feats_2.cols)
		cv_FLANN_Matcher(feats_1, feats_2, matches, good_matches, gm_dist);
}

void findIntersectedFeatures(std::string imgPath, cv::Mat img1, uchar_descriptors& descriptor_1, std::vector<DMatch >& inMatches)
{

	auto img2 = imread(imgPath.c_str());
	auto img3 = imread(imgPath.c_str());

	uchar_descriptors::resizeImage(&img2, 800);
	uchar_descriptors::resizeImage(&img3, 600);

	uchar_descriptors descriptor_2(imgPath.c_str(), img2, "", AKAZE_FEATS);
	uchar_descriptors descriptor_3(imgPath.c_str(), img3, "", AKAZE_FEATS);

	descriptor_1.ExtractAKAZE();
	descriptor_2.ExtractAKAZE();
	descriptor_3.ExtractAKAZE();

	std::vector<DMatch > gm_12, gm_13, gm_23;
	findMatches(descriptor_1, descriptor_2, gm_12);
	findMatches(descriptor_1, descriptor_3, gm_13);
	findMatches(descriptor_2, descriptor_3, gm_23);

	std::vector<int> intersectedMatches;
	int last = 0;
	for (int i = 0; i < gm_12.size(); i++)
	{
		if (gm_12[i].queryIdx < gm_13[last].queryIdx)
			continue;
		for (int k = last; k < gm_13.size(); k++)
		{
			if (gm_12[i].queryIdx == gm_13[k].queryIdx)
			{
				intersectedMatches.push_back(gm_12[i].queryIdx);
				inMatches.push_back(gm_12[i]);
				last = k;
				continue;
			}
		}
	}
}