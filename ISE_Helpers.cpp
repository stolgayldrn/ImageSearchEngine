#include "ISE_Helpers.h"


/** @brief Release all pointers and these are initialized for raw image


@param vwI pointers for visual words 
@param myJSON Pointer
@param vwS String with concatenated visual words 
@param ES_id The pointer for Elastic Search Object ID
*/
void releaseMemoryForRawIm(unsigned* vwI, json_t* myJSON, std::string& vwS, const char* ES_id)
{
	delete[] vwI;
	json_decref(myJSON);
	vwS = "";
	delete ES_id;
}

/** @brief Initialize parameters for address in Path type variable and varibles for Elastic Search

@param myPath Struct for holding all address
@param myES Struct for holding variables of Elastic Search
*/
void initializeParameters(Path &myPath, ELKParams &myES)
{
	myPath.dscFoldName = "dsc_akaze2";
	myPath.DataSet = "C:/ImageSearch/ImageDataSet/03CustomText";
	myPath.imgFoldName = "images";
	myPath.subFolderingLevel = 2;
	myPath.logPath = "D:/v2/SharedV2";
	myPath.corruptedFilesFolder = "D:/v2/SharedV2/corruptedFiles";// for AA server
	//myPath.corruptedFilesFolder = ""; //for local PC
	
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
/** @brief Configure ImageInfo type variable

@param vList The Vector for holding image name list
@param m The index number for current image in image name list
@param imagePath Address of image file
@param myIm ImageInfo variable
@param imgPath Is it image file or dsc file? True for image, False for dsc file
*/
void configureImageConfig(std::vector<std::string> vList, int m, std::string imagePath, ImageInfo& myIm, bool imgPath)
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
/** @brief Configure ImageInfo type variable

@param imagePath Address of image file
@param imageFileName Image file name
@param myIm ImageInfo variable
@param imgPath Is it image file or dsc file? True for image, False for dsc file
*/
void configureImageConfig(std::string imagePath, std::string imageFileName, ImageInfo& myIm, bool imgPath)
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
/** @brief Import descriptor which is extracted from image or read from dsc file to Elastic Search index.

1. Quantize descriptors 
2. Create JSON 
3. Commit Elasticsearch
4. Release memory

@param myPath Path type variable
@param myES ELKParams type variable
@param VT The vocabulary tree for quantization
@param dscPath The address of descriptor ( Use "" if the descriptor is not read from a dsc file)
@param myIm ImageInfo variable
@param my_desc The descriptor variable
*/
void importRawImage(Path myPath, ELKParams myES, TVoctreeVLFeat* VT, std::string dscPath,
	ImageInfo myIm, UcharDescriptors * my_desc)
{
	unsigned int * vwI = new unsigned int[my_desc->getNumOfDescriptors()];
	json_t* myJSON = json_object();
	std::string vwS = "";
	const char * ES_id = new    char;

	if (my_desc->getNumOfDescriptors() > 0)
	{
		try
		{
			VT->quantize_multi(vwI, my_desc->getUcharDescriptors(), my_desc->getNumOfDescriptors(), my_desc->getFeatureSize());
			for (unsigned int s = 0; s < my_desc->getNumOfDescriptors(); s++)
				vwS += " " + int2string(int(vwI[s]));

			if (vwS != "")
			{
				getJsonForNewImage(myIm, myPath, vwS, myJSON);
				commitJsonToELK(myES, myJSON, myIm.fileName.c_str(), ES_id);
			}
		}
		catch (std::exception e)
		{
			printf("\nElasticSearch:::commit error:%s", e.what());
		}
	}
	try
	{
		releaseMemoryForRawIm(vwI, myJSON, vwS, ES_id);
	}
	catch (std::exception e)
	{
		printf("\nElasticSearch:::release error:%s", e.what());
	}
}
/** @brief Create the normalized Elastic Search Score from scoresELK

@param scoresELK INPUT::Score's list from Elastic Search. 
@param totalNumELK INPUT::Total Number of Returned Scores from Elasticsearch.
@param normELKScore OUTPUT::Score'S list of normalized Elastic Search Scores.
*/
void calculateNormalizedELKScore(std::vector<float>& scoresELK, int totalNumELK, std::vector<float>& normELKScore)
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
/** @brief Create a query from the descriptor and send to Elasticsearch.

1. Quantize the descriptor with Vocabulary Tree
2. Create JSON 
3. POST Query to Elastic Search
4. Calculate Post Process score
5. Push the result to result vector list
6. Release memory

@param myPath Path type variable for holding all address about descriptor
@param myES Elastic Search parameter to access the Elastic Search Index
@param VT The vocabulary tree for quantization
@param dscPath The address of descriptor ( Use "" if the descriptor is not read from a dsc file)
@param myIm ImageInfo variable
@param my_desc The descriptor variable
@param returnedFileNames OUTPUT::File name list of returned images from Elastic Search in order to descanding Elastic Search score
@param scoresPP			 OUTPUT::Score list of post process
@param scoresELK		 OUTPUT::Score list of Elastic Search
*/
void queryRawImage(Path myPath, ELKParams myES, TVoctreeVLFeat* VT, std::string dscPath, ImageInfo myIm,
                   UcharDescriptors *my_desc, std::vector<std::string> & returnedFileNames, std::vector<float> &scoresPP, std::vector<float> &scoresELK)
{
	unsigned int * vwI = new unsigned int[my_desc->getNumOfDescriptors()];
	json_t* myJSON = json_object();
	std::string vwS = "";
	const char * ES_id = new char;
	std::vector<std::string> dscPathsV;
	int totalNumELK;

	if (my_desc->getNumOfDescriptors() > 0)
	{
		try
		{
			VT->quantize_multi(vwI, my_desc->getUcharDescriptors(), my_desc->getNumOfDescriptors(), my_desc->getFeatureSize());
			for (unsigned int s = 0; s < my_desc->getNumOfDescriptors(); s++)
				vwS += " " + int2string(int(vwI[s]));
			/*Temporary words writing*/
			/*std::vector<std::vector<int>> tempWords;
			for (int v = 0; v < my_desc->getNumOfDescriptors(); v++)
			{
				std::vector<int>tempW;
				tempW.push_back(vwI[v]);
				tempWords.push_back(tempW);
			}
			std::string tempPath = myIm.path + ".csv";
			writeToCSV(tempWords,"a.csv", 1);*/
			/*************************/
			if (vwS != "")
			{
				getJsonForQueryImage(vwS, "words_string", myJSON);
				
				//TODO: add scores at postQueryToELK
				postQueryToELK(myES, myJSON, myIm, returnedFileNames, dscPathsV, scoresELK, totalNumELK);
				//ELK score filtering//////////////////////////
				if (totalNumELK > 0)
				{
					std::vector<float> normELKScore;
					calculateNormalizedELKScore(scoresELK, totalNumELK, normELKScore);
					// Post Process
					postProcess(my_desc, returnedFileNames, scoresPP, normELKScore.size());
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
		releaseMemoryForRawIm(vwI, myJSON, vwS, ES_id);
	}
	catch (std::exception e)
	{
		printf("\nElasticSearch:::release error:%s", e.what());
	}
}
/** @brief Calculate Post Process scoress based on paper "FAST GEOMETRIC RE-RANKING FOR IMAGE-BASED RETRIEVAL"
by Sam S. Tsai1, David Chen1, Gabriel Takacs1, Vijay Chandrasekhar1, Ramakrishna Vedantham2, Radek Grzeszczuk2, and Bernd Girod1

@param query The Descriptor for query
@param fileNameVec File name list whose post process scores are calculated by query
@param scores OUTPUT:: Post process scores list
@param numELK INPUT:: total number of returned images from Elastic Search
*/
int postProcess(UcharDescriptors *query, std::vector<std::string> fileNameVec, std::vector<float> & scores, int numELK)
{
	std::vector<float> scoresSorted;
	std::vector<int> scoreRank;
	Mat descriptorQuery;
	query->getCopyOfOpencvDescriptors(descriptorQuery);
	std::vector<Point2f> coordsQuery = query->getCoords();
	//std::vector<float> oriQuery, scaleQuery;

#if defined _OPENMP
	#pragma omp parallel for ordered schedule(dynamic)
#endif
	for (int i = 0; i<numELK; i++)
	{
		if (isDscFile(fileNameVec[i].c_str()) || isImageFile(fileNameVec[i].c_str()))
		{
			UcharDescriptors matchDsc(fileNameVec[i].c_str(), fileNameVec[i].c_str(), AKAZE_FEATS);
			Mat descriptorMatch;
			if (isDscFile(fileNameVec[i].c_str()))
			{
				matchDsc.readDSC();
				descriptorMatch = Mat(matchDsc.getNumOfDescriptors(), matchDsc.getFeatureSize(), CV_8UC1);;
				matchDsc.getReadModeDescriptors(descriptorMatch);
			}
			else if (isImageFile(fileNameVec[i].c_str()))
			{
				matchDsc.setResizeImage(true);
				matchDsc.extractAKAZE();
				matchDsc.getCopyOfOpencvDescriptors(descriptorMatch);
			}
			std::vector<Point2f> coordsMatch = matchDsc.getCoords();

			// Match with FLANN
			std::vector<DMatch > matches, good_matches;
			if (descriptorQuery.type() == 0)  descriptorQuery.convertTo(descriptorQuery, CV_32F);
			if (descriptorMatch.type() == 0)  descriptorMatch.convertTo(descriptorMatch, CV_32F);
			float gm_dist = (query->getFeatureType() == AKAZE_FEATS) ? GOOD_MATCHES_DISTANCE__AKAZE : GOOD_MATCHES_DISTANCE__EZSIFT;
			if (descriptorQuery.type() == descriptorMatch.type() && descriptorQuery.cols == descriptorMatch.cols)
				cv_FLANN_Matcher(descriptorQuery, descriptorMatch, matches, good_matches, gm_dist);
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
			printf("Unable to read file %s. Skipping...\n", fileNameVec[i].c_str());
		}
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
/** @brief Create Images by splitting input images by maxSize parameter.

1. If image size is bigger than maxSize, then the first image is downsampled image of hole 
image and named as "resized".

2. The rest of images are cropped image in order to left to right, up to down.


@note The 2. type of images are comment out for a speacial case, 
please clear the comment for proper use and remove this note.

@param Input Raw newspaper image
@param OutputVector List of created images
@param OutputNames List of created image names
@param maxSize The input parameter for splitting images
*/
void splitUpNewspaperImages(const cv::Mat Input, std::vector<cv::Mat>& OutputVector, std::vector<std::string>& OutputNames, int maxSize)
{
	if (Input.rows > maxSize || Input.cols > maxSize)
	{
		int h = Input.rows;
		int w = Input.cols;
		double hs = h*1.0 / maxSize;
		double ws = w*1.0 / maxSize;
		cv::Mat reszImg;
		if (h>w)
			cv::resize(Input, reszImg, cv::Size(w / hs, h / hs), 0, 0, CV_INTER_LINEAR);
		else
			cv::resize(Input, reszImg, cv::Size(w / ws, h / ws), 0, 0, CV_INTER_LINEAR);

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
/** @brief Caluculate weighted scores

	W = ( w1 * scoresELK ) + ( w2 * scoresW )

	w1 = 0.3
	
	w2 = 0.7

@param scoresPP Post Process score
@param scoresELK Elastic Search score
@param scoresW OUTPUT:: Weighted Score
*/
void scoreWeighting(std::vector<float> scoresPP, std::vector<float> scoresELK, std::vector<float> & scoresW)
{
	for (unsigned int l = 0; l < scoresPP.size(); l++)
	{
		float wScore = (0.3 * scoresELK[l]) + (0.7 * scoresPP[l]);
		scoresW.push_back(wScore);
	}
}
/** @brief Returns current date and time 


*/
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
/** @brief Write status of error with current date time to the log.txt

@param myPath Struct for holding all address
@param ofs Output File Stream variable
@param imgPathDUMP Address of the file with error
@param errCmd The error command
*/
void printErrorToLog(Path myPath, std::ofstream &ofs, std::basic_string<char> imgPathDUMP, std::string errCmd)
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
