#include "ISE_Helpers.h"


void ReleaseAll__ImRawIm(unsigned* vwI, json_t* myJSON, string& vwS, const char* ES_id)
{
	delete[] vwI;
	json_decref(myJSON);
	vwS = "";
	delete ES_id;
}

void ReleaseAll__ImRawIm(unsigned* vwI, unsigned* vwI_low, json_t* myJSON, string& vwS, string& vwS_low, const char* ES_id)
{
	delete[] vwI;
	delete[] vwI_low;
	json_decref(myJSON);
	vwS = "";
	vwS_low = "";
	delete ES_id;
}

void WriteCSV(vector<vector<string>> testSet, char* fileName)
{
	ofstream myfile;
	string CSV_Path = fileName;
	myfile.open(CSV_Path.c_str());
	for (int i = 0; i < testSet.size(); i++){
		string lineStr = "";
		for (int ii = 0; ii < 276; ii++){
			lineStr += testSet[i][ii] + ";";
		}
		lineStr += "\n";
		myfile << lineStr;
	}
}

void paramsConfig(Path &myPath, ES_params &myES)
{
	myPath.dscFoldName = "dsc_akaze2";
	myPath.dscFoldName2 = "dsc_akaze_low";
	myPath.DataSet = "C:/ImageSearch/TEST_IM_33";
	myPath.imgFoldName = "images";
	myPath.subFolderingLevel = 2;
	myPath.VocTree = "C:/ImageSearch/VT_Trees/VT_flicker500K_AKAZE_middle_tree_S2_P";
	myPath.VocTreeLow = "C:/ImageSearch/VT_Trees/VT_flicker500K_AKAZE_small_tree_S2_P";

	myES.index = "flicker1m_test2";
	myES.type = "akaze";
	myES.url = "http://172.16.10.202:9200";
	myES.userPWD = "";
}

void ImageConfig(vector<string> vList, int m, string imagePath, Image_Info& myIm)
{
	myIm.dataSet = "flicker1M";
	myIm.dataSubSet = "";
	myIm.descriptorType = "akaze";
	myIm.encoding = "jpg";
	//myIm.fileName = vList[m];
	myIm.fileName = vList[m].substr(0, vList[m].length() - 4);
	myIm.height = 0.0;
	myIm.width = 0.0;
	myIm.Import = "true";
	myIm.Query = "false";
	myIm.path = imagePath;
	myIm.source_type = "flicker1M";
}

void ImportRawImage(Path myPath, ES_params myES, TVoctreeVLFeat* VT, string dscPath,
	Image_Info myIm, uchar_descriptors * my_desc)
{
	unsigned int * vwI = new unsigned int[my_desc->get_num_descriptors()];
	json_t* myJSON = json_object();
	string vwS = "";
	const char * ES_id = new    char;

	if (my_desc->get_num_descriptors() > 0)
	{
		try
		{
			VT->quantize_multi(vwI, my_desc->get_data(), my_desc->get_num_descriptors(), 61);
			for (unsigned int s = 0; s < my_desc->get_num_descriptors(); s++)
				vwS += " " + int2string(int(vwI[s]));

			if (vwS != "")
			{
				getJSON_new_image(myIm, myPath, myJSON, vwS);
				ES_commit(myES, myJSON, ES_id, myIm.fileName.c_str());
			}
		}
		catch (exception e)
		{
			printf("\nElasticSearch:::commit error:%s", e.what());
		}
	}
	try
	{
		ReleaseAll__ImRawIm(vwI, myJSON, vwS, ES_id);
	}
	catch (exception e)
	{
		printf("\nElasticSearch:::release error:%s", e.what());
	}
}

void ImportRawImage(Path myPath, ES_params myES, TVoctreeVLFeat* VT, string dscPath, Image_Info myIm,
	uchar_descriptors * my_desc, TVoctreeVLFeat* VT_low, string dsc2Path, uchar_descriptors * my_desc2)
{
	unsigned int * vwI = new unsigned int[my_desc->get_num_descriptors()];
	unsigned int * vwI_low = new unsigned int[my_desc2->get_num_descriptors()];
	json_t* myJSON = json_object();
	string vwS = "";
	string vwS_low = "";
	const char * ES_id = new char;

	if (my_desc->get_num_descriptors() > 0 && my_desc2->get_num_descriptors() > 0)
	{
		try
		{
			VT->quantize_multi(vwI, my_desc->get_data(), my_desc->get_num_descriptors(), 61);
			VT_low->quantize_multi(vwI_low, my_desc2->get_data(), my_desc2->get_num_descriptors(), 61);
			for (unsigned int s = 0; s < my_desc->get_num_descriptors(); s++)
				vwS += " " + int2string(int(vwI[s]));

			for (unsigned int s = 0; s < my_desc2->get_num_descriptors(); s++)
				vwS_low += " " + int2string(int(vwI_low[s]));

			if (vwS != "")
			{
				getJSON_new_image(myIm, myPath, myJSON, vwS, vwS_low);
				ES_commit(myES, myJSON, ES_id, myIm.fileName.c_str());
			}
		}
		catch (exception e)
		{
			printf("\nElasticSearch:::commit error:%s", e.what());
		}
	}
	try
	{
		ReleaseAll__ImRawIm(vwI, vwI_low, myJSON, vwS, vwS_low, ES_id);
	}
	catch (exception e)
	{
		printf("\nElasticSearch:::release error:%s", e.what());
	}
}

void QueryRawImage(Path myPath, ES_params myES, TVoctreeVLFeat* VT, string dscPath, Image_Info myIm,
	uchar_descriptors my_desc, vector<string> &testSet, string &returnFileName)
{
	unsigned int * vwI = new unsigned int[my_desc.get_num_descriptors()];
	json_t* myJSON = json_object();
	string vwS = "";
	const char * ES_id = new char;
	vector<string> dscPathsV;
	vector<float> scoresV;

	if (my_desc.get_num_descriptors() > 0)
	{
		try
		{
			VT->quantize_multi(vwI, my_desc.get_data(), my_desc.get_num_descriptors(), 61);
			for (unsigned int s = 0; s < my_desc.get_num_descriptors(); s++)
				vwS += " " + int2string(int(vwI[s]));

			if (vwS != "")
			{
				printf("nok \n");
				getJSON_query_image(myJSON, vwS, "words_string");
				//TODO: add scores at ES_post_query
				ES_post_query(myES, myJSON, myIm, testSet, dscPathsV, scoresV);
				post_process_step(my_desc, dscPathsV, testSet, scoresV, 10);

				printf("ok \n");
			}
		}
		catch (exception e)
		{
			printf("\nElasticSearch:::commit error:%s", e.what());
		}
	}
	try
	{
		ReleaseAll__ImRawIm(vwI, myJSON, vwS, ES_id);
	}
	catch (exception e)
	{
		printf("\nElasticSearch:::release error:%s", e.what());
	}
}

int post_process_step(uchar_descriptors query, vector<string> dscV, vector<string> fnV, vector<float> scoresV, int numEsReturns)
{
	vector<float> scores;
	std::vector<float> scoresSorted;
	std::vector<int> scoreRank;
	Mat descriptorQuery;
	query.get_descriptors(descriptorQuery);
	std::vector<Point2f> coordsQuery = query.getCoords();
	//std::vector<float> oriQuery, scaleQuery;
	//TODO: runOptions Create
	//TODO: runOptions Create
	//int numCorr = runOptions.numMatches;
	//scores.reserve(numCorr);
	//scores.resize(numCorr);

#if defined _OPENMP
	//#pragma omp parallel for ordered schedule(dynamic)
#endif
	for (int i = 0; i<numEsReturns; i++)
	{
		//int threadId = omp_get_thread_num();
		// read the signature of the descriptor
		if (IS_DscFile(dscV[i].c_str()))
		{
			uchar_descriptors matchDsc(dscV[i].c_str(), AKAZE_FEATS);
			matchDsc.read_dsc();
			Mat descriptorMatch;
			matchDsc.get_descriptors(descriptorMatch);
			vector<Point2f> coordsMatch = matchDsc.getCoords();

			// Match with FLANN
			std::vector<DMatch > matches, good_matches;
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
	scores.clear();

	return 0;
}
