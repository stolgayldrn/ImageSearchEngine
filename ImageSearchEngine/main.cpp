#include "ISE_Helpers.h"

void CreateThreeLowTree(uchar_descriptors descs_1, unsigned* vwI_1, unsigned* vwI_1_t5, unsigned* vwI_1_t4, unsigned* vwI_1_t3);
void calcTop100VT(Mat hist, int top100[100], int top100M[100]);
void computeHistOfVT(cv::Mat descs1_mat, Mat hist);
void sortTopVT(cv::Mat hist_1, int top100_d1[100], int top100M_d1[100]);
cv::Mat DescriptorAnalyzer(std::string imgPath, std::string imgOrgPath, FeatureType feature_type, TVoctreeVLFeat* VT);

void swapBlocks(std::vector<std::string> & outputLine, std::vector<float> & scoresPP, int m);
void blockSorting(std::vector<float> & scoresPP, std::vector<std::string> & outputLine);
void findMatches(uchar_descriptors &descriptor_1, uchar_descriptors &descriptor_2, std::vector<DMatch >& good_matches);
void downsamplingAnalyzer(std::string imgPah, double size1, double size2, std::vector<DMatch >& good_matches);
void findIntersectedFeatures(std::string imgPath, cv::Mat img1, uchar_descriptors& descriptor_1, std::vector<DMatch >& inMatches);

int main(int argc, char *argv[])
{
	Path myPath;
	ELK_params myES;
	TVoctreeVLFeat* VT = new TVoctreeVLFeat;
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

	std::string myPathDUMP;
	if (argc >= 2)
	myPathDUMP = argv[2];

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	int numThread = omp_get_max_threads();
	omp_set_num_threads(numThread); // Use max threads for all consecutive parallel regions
	std::vector<std::string> monthList;
	GET_FolderList(myPath.DataSet.c_str(), monthList);
	for (size_t m = 0; m < monthList.size(); m++)
	{
	try{
	auto mPath = myPath.DataSet + "/" + monthList[m];
	auto mPathDUMP = myPathDUMP + "/" + monthList[m];
	std::vector<std::string> dayList;
	GET_FolderList(mPath.c_str(), dayList);
	for (size_t d = 0; d < dayList.size(); d++)
	{
	std::ofstream ofs(myPath.logPath + "/log.txt", std::ios::app);
	auto dPath = mPath + "/" + dayList[d];
	auto dPathDUMP = mPathDUMP + "/" + dayList[d];
	std::vector<std::string> imgList;
	GET_DirectoryImages(dPath.c_str(), imgList);

	#pragma omp parallel for
	for (int i = 0; i < imgList.size(); i++)
	{
	auto imgPath = dPath + "/" + imgList[i];
	auto imgPathDUMP = dPathDUMP + "/" + imgList[i];
	Image_Info myIm;
	//ImageConfig(imgList, i, imgPath.c_str(), myIm, true);
	ImageConfig(imgList, i, imgPathDUMP.c_str(), myIm, true);
	if (!FileExist(imgPathDUMP.c_str()))
	continue;
	if (IS_ImageFile(imgPathDUMP.c_str()))
	{
	try
	{
	printf("nok\n");
	uchar_descriptors myDesc(imgPathDUMP.c_str(), "", AKAZE_FEATS);
	myDesc.setResizeImage(true);
	int excCode = myDesc.ExtractAKAZE();
	if (excCode == 1){
	myIm.height = myDesc.GetImageHeight();
	myIm.width = myDesc.GetImageWidth();
	myIm.numDescs = int(myDesc.GetNumOfDescriptors());
	ImportRawImage(myPath, myES, VT, imgPathDUMP, myIm, &myDesc);
	printf("ok\n");
	}
	else if (excCode == 3)
	{
	// 1. Save corrupted files
	printErrorToLog(myPath, ofs, imgList, i, imgPathDUMP, myDesc,
	"Failed at first attempt: corrupted file or small image");
	cv::Mat imgToWrite;
	myDesc.GetImage__Copy(imgToWrite);
	cv::imwrite(myPath.corruptedFilesFolder + "/" + imgList[i], imgToWrite);
	// 2. Try to extract feature again
	for (int r = 0; r < 6; r++)
	{
	uchar_descriptors myDesc2(imgPathDUMP.c_str(), "", AKAZE_FEATS);
	myDesc2.setResizeImage(true);
	int excCode2 = myDesc2.ExtractAKAZE();
	if (excCode2 == 1){
	myIm.height = myDesc2.GetImageHeight();
	myIm.width = myDesc2.GetImageWidth();
	myIm.numDescs = int(myDesc2.GetNumOfDescriptors());
	ImportRawImage(myPath, myES, VT, imgPathDUMP, myIm, &myDesc2);
	printf("ok at %d. try\n", r + 2);
	break;
	}
	else if (excCode2 == 3)
	{
	printf("failed at %d. try:   %s\n", r + 2, imgList[i].c_str());
	printErrorToLog(myPath, ofs, imgList, i, imgPathDUMP, myDesc,
	"Failed at " + int2string(r + 2) + " attemp: corrupted file or small image");
	}
	}
	}
	}
	catch (std::exception e)
	{
	printf("\nMain:::extract akaze and write:%s", e.what());
	}
	}
	}
	ofs.close();
	}
	}
	catch (std::exception e)
	{
	printf("\nMain:::%s", e.what());
	}
	}

	//*/
	// AA Query Images
	/*
	if (argc >= 1 && argv[1] != nullptr)
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
				if (Image.rows > MAX_IMAGE_SIZE || Image.cols > MAX_IMAGE_SIZE)
					ImageSpliter(Image, imagesVector, imageNameList, MAX_IMAGE_SIZE);
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
					my_desc_f.setResizeImage(true);
					int excCode = my_desc_f.ExtractAKAZE();
					myIm.height = my_desc_f.GetImageHeight();
					myIm.width = my_desc_f.GetImageWidth();
					myIm.numDescs = int(my_desc_f.GetNumOfDescriptors());
					if (excCode == 1)
					{
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

						blockSorting(scoresPP, outputLine);
						
						outputsV.push_back(outputLine);
						outputLine.clear(); outputLine.shrink_to_fit();
						// 5. Releases
						testIm.clear(); testIm.shrink_to_fit();
						scoresPP.clear(); scoresPP.shrink_to_fit();
						scoresELK.clear(); scoresELK.shrink_to_fit();
						scoresW.clear(); scoresW.shrink_to_fit();
					}
					else if (excCode == 3)
					{
						outputLine.push_back(imgList[i].c_str());
						outputLine.push_back(imageNameList[f].c_str());
						outputLine.push_back("excCode3");
						for (int e = 0; e < QUERY_RETURN_SIZE * 4; e++)
							outputLine.push_back(int2string(-3));
						outputsV.push_back(outputLine);
						outputLine.clear(); outputLine.shrink_to_fit();
					}
					printf("ok\n");
				}
			}
			catch (std::exception e)
			{
				printf("\nMain:::extract akaze and write: %s", e.what());
			}
		}
		else
			printf("\nNot an Image File: %s", imgPath.c_str());
	}
	try
	{
		//std::string wscPath = "D:/v2/SharedV2/custom03.csv";
		//if (argc >= 2)
		//	wscPath = argv[2];
		WriteCSV(outputsV, "Test03.csv", 3 + (4 * QUERY_RETURN_SIZE));
	}
	catch (std::exception e)
	{
		printf("\nCSV writing error: %s", e.what());
	}
	//*/
	// AA Query One Image
	/*

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	int numThread = omp_get_max_threads();
	omp_set_num_threads(numThread); // Use 4 threads for all consecutive parallel regions

	std::string imgPath = "C:/ImageSearch/ImageDataSet/DebugImg/";
	std::string imgFN = "8459726.jpg";
	imgPath += imgFN;
	Image_Info myIm;
	ImageConfig(imgPath, imgFN, myIm, true);
	std::vector<std::vector<std::string>> writeCsvFN;
	std::vector<std::vector<float>> writeCsvScoreElk;

	if (IS_ImageFile(imgPath.c_str()))
	{
	try
	{
	cv::Mat Image = cv::imread(imgPath, cv::IMREAD_GRAYSCALE);
	std::vector<std::string> returnedFN;
	std::vector<float> scoresPP, scoresELK, scoresW;
	printf("nok\n");

	// 1. Feature Extraction
	for (size_t i = 0; i < 10; i++)
	{
	uchar_descriptors my_desc_f(imgPath.c_str(), "", AKAZE_FEATS);
	my_desc_f.setResizeImage(true);
	int excCode = my_desc_f.ExtractAKAZE();
	myIm.height = my_desc_f.GetImageHeight();
	myIm.width = my_desc_f.GetImageWidth();
	if (excCode == 1)
	{
	// 2. ELK Query and Response
	QueryRawImage(myPath, myES, VT, "NoDSC", myIm, &my_desc_f, returnedFN, scoresPP, scoresELK);

	}
	else if (excCode == 3)
	{
	printf("excution error\n");
	}
	printf("ok\n");
	writeCsvFN.push_back(returnedFN);
	writeCsvScoreElk.push_back(scoresELK);
	}
	WriteCSV(writeCsvFN, "C:/Users/m.alsadi/Desktop/test2.csv", QUERY_RETURN_SIZE);
	WriteCSV(writeCsvScoreElk, "C:/Users/m.alsadi/Desktop/test2ElkSc.csv", QUERY_RETURN_SIZE);

	}
	catch (std::exception e)
	{
	printf("\nMain:::extract akaze and write: %s", e.what());
	}
	}
	else
	printf("\nNot an Image File: %s", imgPath.c_str());


	//*/
	// Feature Marking
	/*
	std::string imgPath = "C:\\Users\\m.alsadi\\Desktop\\downsample\\BL.png";
	std::string imgOrgPath = "C:\\Users\\m.alsadi\\Desktop\\downsample\\BL5.png";
	std::vector<int> compression_params;
	compression_params.push_back(CV_IMWRITE_PNG_COMPRESSION);
	compression_params.push_back(9);
	FeatureType feature_type;
	feature_type = AKAZE_FEATS;
	cv::Mat outImg =  DescriptorAnalyzer(imgPath, imgOrgPath, feature_type, VT);

	try
	{
		imwrite("red4_akaze.png", outImg, compression_params);
	}
	catch (std::runtime_error& ex)
	{
		fprintf(stderr, "Exception converting image to PNG format: %s\n", ex.what());
		return 1;
	}

	//*/
	// Downsampling Analyzer
	///*
	std::vector<DMatch > gm_12, gm_13, gm_23;
	std::string imgPath = "C:\\Users\\m.alsadi\\Desktop\\DownSamplingAnalyzer\\1.jpg";
	std::vector<DMatch > inMatches;
	auto img1 = imread(imgPath.c_str());
	uchar_descriptors::resizeImage(&img1, 1200);
	uchar_descriptors descriptor_1(imgPath.c_str(), img1, "", AKAZE_FEATS);
	findIntersectedFeatures(imgPath.c_str(),img1, descriptor_1, inMatches);
	/*downsamplingAnalyzer(imgPath, 1200, 800, gm_12);
	downsamplingAnalyzer(imgPath, 1200, 600, gm_13);
	downsamplingAnalyzer(imgPath, 800,  600, gm_23);*/
    //*/
	duration = clock() - start;
	delete VT;
	printf(":::Duration: %.2f, Error count: %d", duration, err_count);
	return 0;
}


void CreateThreeLowTree(uchar_descriptors descs_1, unsigned* vwI_1, unsigned* vwI_1_t5, unsigned* vwI_1_t4, unsigned* vwI_1_t3)
{
	for (int i = 0; i < descs_1.GetNumOfDescriptors(); i++)
	{
		unsigned int temp = vwI_1[i];
		unsigned int length = ceil(log10(temp + 1));
		switch (length)
		{
		case 6:
			vwI_1_t5[i] = temp / 10;
			vwI_1_t4[i] = temp / 100;
			vwI_1_t3[i] = temp / 1000;
			break;
		case 5:
			vwI_1_t5[i] = temp;
			vwI_1_t4[i] = temp / 10;
			vwI_1_t3[i] = temp / 100;
			break;
		case 4:
			vwI_1_t5[i] = temp;
			vwI_1_t4[i] = temp;
			vwI_1_t3[i] = temp / 10;
			break;
		case 3:
			vwI_1_t5[i] = temp;
			vwI_1_t4[i] = temp;
			vwI_1_t3[i] = temp;
			break;
		}
	}
}

void calcTop100VT(Mat hist, int top100[100], int top100M[100])
{
	for (int i = 0; i < 1000000; i++)
	{
		int m = -1;
		for (int t = 99; t >= 0; t--)
		{
			if (int(hist.at<float>(i)) >= top100[t] && int(hist.at<float>(i))>0)
			{
				m = t;
				break;
			}
		}
		if (m > 0)
		{
			for (int t = 0; t < m; t++)
			{
				top100[t] = top100[t + 1];
				top100M[t] = top100M[t + 1];
			}
			top100[m] = int(hist.at<float>(i));
			top100M[m] = i;
		}
	}
}

void computeHistOfVT(cv::Mat descs1_mat, Mat hist)
{
	int nimages = 1; // Only 1 image, that is the Mat scene.
	int channels[] = { 0 }; // Index for hue channel
	int dims = 1;// Only 1 channel, the hue channel
	int histSize = 1000000; // 9 bins, 1 each for Red, RY, Yellow, YG etc.
	float hranges[] = { 0, 1000000 }; // hue varies from 0 to 179, see cvtColor
	const float *ranges[] = { hranges };

	// Compute the histogram.
	cv::calcHist(&descs1_mat,
		nimages,
		channels,
		Mat(), // No mask
		hist, dims, &histSize, ranges, true);
}

void sortTopVT(cv::Mat hist_1, int top100_d1[100], int top100M_d1[100])
{
	for (int i = 0; i < 1000000; i++)
	{
		int m = -1;
		for (int t = 99; t >= 0; t--)
		{
			if (int(hist_1.at<float>(i)) >= top100_d1[t] && int(hist_1.at<float>(i))>0)
			{
				m = t;
				break;
			}
		}
		if (m > 0)
		{
			for (int t = 0; t < m; t++)
			{
				top100_d1[t] = top100_d1[t + 1];
				top100M_d1[t] = top100M_d1[t + 1];
			}
			top100_d1[m] = int(hist_1.at<float>(i));
			top100M_d1[m] = i;
		}
	}
}

cv::Mat DescriptorAnalyzer(std::string imgPath, std::string imgOrgPath, FeatureType feature_type, TVoctreeVLFeat* VT)
{
	uchar_descriptors descs_1(imgPath.c_str(), "", feature_type);
	uchar_descriptors descs_2(imgOrgPath.c_str(), "", feature_type);

	switch (feature_type)
	{
	case AKAZE_FEATS:
		descs_1.setResizeImage(true);
		descs_1.ExtractAKAZE();
		descs_2.setResizeImage(true);
		descs_2.ExtractAKAZE();
		break;
	case EZ_SIFT:
		//descs_1.setResizeImage(true);
		descs_1.ExtractEZSIFT();
		descs_1.ConvertEzsiftToOpencv();
		descs_2.setResizeImage(true);
		descs_2.ExtractEZSIFT();
		descs_2.ConvertEzsiftToOpencv();
		break;
	default:
		//descs_1.setResizeImage(true);
		descs_1.ExtractAKAZE();
		descs_2.setResizeImage(true);
		descs_2.ExtractAKAZE();
		break;
	}

	cv::Mat descs1_dsc, descs2_dsc;
	descs1_dsc = descs_1.GetOpencvDescriptors();
	descs2_dsc = descs_2.GetOpencvDescriptors();
	// Match with FLANN
	std::vector<DMatch > matches, good_matches;
	if (descs1_dsc.type() == 0)  descs1_dsc.convertTo(descs1_dsc, CV_32F);
	if (descs2_dsc.type() == 0)  descs2_dsc.convertTo(descs2_dsc, CV_32F);

	float gm_dist = (feature_type == AKAZE_FEATS) ? 0 : 10;

	if (descs1_dsc.type() == descs2_dsc.type() && descs1_dsc.cols == descs2_dsc.cols)
		cv_FLANN_Matcher(descs1_dsc, descs2_dsc, matches, good_matches, gm_dist);
	/*******************************/
	unsigned int * vwI_1 = new unsigned int[descs_1.GetNumOfDescriptors()];
	/*unsigned int * vwI_1_t5 = new unsigned int[descs_1.GetNumOfDescriptors()];
	unsigned int * vwI_1_t4 = new unsigned int[descs_1.GetNumOfDescriptors()];
	unsigned int * vwI_1_t3 = new unsigned int[descs_1.GetNumOfDescriptors()];*/
	VT->quantize_multi(vwI_1, descs_1.GetUCHAR_descriptors(), descs_1.GetNumOfDescriptors(), descs_1.GetFeatureSize());
	//CreateThreeLowTree(descs_1, vwI_1, vwI_1_t5, vwI_1_t4, vwI_1_t3);


	unsigned int * vwI_2 = new unsigned int[descs_2.GetNumOfDescriptors()];
	/*unsigned int * vwI_2_t5 = new unsigned int[descs_2.GetNumOfDescriptors()];
	unsigned int * vwI_2_t4 = new unsigned int[descs_2.GetNumOfDescriptors()];
	unsigned int * vwI_2_t3 = new unsigned int[descs_2.GetNumOfDescriptors()];*/
	VT->quantize_multi(vwI_2, descs_2.GetUCHAR_descriptors(), descs_2.GetNumOfDescriptors(), descs_2.GetFeatureSize());
	//CreateThreeLowTree(descs_2, vwI_2, vwI_2_t5, vwI_2_t4, vwI_2_t3);

	cv::Mat descs1_mat = cv::Mat(descs_1.GetNumOfDescriptors(), 1, CV_32S, vwI_1);
	cv::Mat descs2_mat = cv::Mat(descs_2.GetNumOfDescriptors(), 1, CV_32S, vwI_2);

	if (descs1_mat.type() == 4)
		descs1_mat.convertTo(descs1_mat, CV_32F);
	if (descs2_mat.type() == 4)
		descs2_mat.convertTo(descs2_mat, CV_32F);

	//computeHistOfVT(descs1_mat, hist_1);
	cv::Mat hist_1, hist_2;
	int nimages = 1; // Only 1 image, that is the Mat scene.
	int channels[] = { 0 }; // Index for hue channel
	int dims = 1;// Only 1 channel, the hue channel
	int histSize = 1000000; // 9 bins, 1 each for Red, RY, Yellow, YG etc.
	float hranges[] = { 0, 1000000 }; // hue varies from 0 to 179, see cvtColor
	const float *ranges[] = { hranges };

	// Compute the histogram.
	cv::calcHist(&descs1_mat, nimages, channels, Mat(), hist_1, dims, &histSize, ranges, true);
	cv::calcHist(&descs2_mat, nimages, channels, Mat(), hist_2, dims, &histSize, ranges, true);

	//computeHistOfVT(descs1_mat, hist_2);

	// Show the calculated histogram in command window

	/*for (int h = 0; h < histSize; h++)
	{
	float binVal = hist.at<float>(h);
	std::cout << " " << binVal;
	}*/

	//// Plot the histogram
	//int hist_w = 2000; int hist_h = 400;
	//int bin_w = cvRound((double)hist_w / histSize);

	//Mat histImage(hist_h, hist_w, CV_8UC1, Scalar(0, 0, 0));
	//cv::normalize(hist, hist, 0, histImage.rows, NORM_MINMAX, -1, Mat());

	//for (int i = 1; i < histSize; i++)
	//{
	//	line(histImage, Point(bin_w*(i - 1), hist_h - cvRound(hist.at<float>(i - 1))),
	//		Point(bin_w*(i), hist_h - cvRound(hist.at<float>(i))),
	//		Scalar(255, 0, 0), 2, 8, 0);
	//}

	//cv::namedWindow("Result", 1);    
	//cv::imshow("Result", histImage);
	//cv::waitKey();
	// Now hist will contain the counts in each bin.
	// Lets just print out the values. Note that you can output Mat using std::cout
	//std::cout << "Histogram: " << std::endl << hist << std::endl;

	int top100_d1[100] = { 0 };
	int top100M_d1[100] = { 0 };

	int top100_d2[100] = { 0 };
	int top100M_d2[100] = { 0 };
	//calcTop100VT(*hist_1, top100, top100M);

	sortTopVT(hist_1, top100_d1, top100M_d1);
	sortTopVT(hist_2, top100_d2, top100M_d2);

	FlannBasedMatcher fl_matcher;
	std::vector< DMatch > fl_matches, fl_good_matches;
	if (descs1_mat.type() == descs2_mat.type() && descs1_mat.cols == descs2_mat.cols)
		cv_FLANN_Matcher(descs1_dsc, descs2_dsc, fl_matches, fl_good_matches, 400);

	/*BFMatcher bf_matcher(NORM_L2);
	std::vector<std::vector<DMatch> > bf_matches;
	bf_matcher.knnMatch(descs1_mat, descs2_mat, bf_matches, 2);*/


	cv::Mat outImg4DrawMatches;
	cv::drawMatches(descs_1.GetImageMat(), descs_1.GetOpencvKeypoints(), descs_2.GetImageMat(),
		descs_2.GetOpencvKeypoints(), good_matches, outImg4DrawMatches);
	cv::imshow("Draw matches", outImg4DrawMatches);

	cv::Mat outImg4DrawMatches_flann;
	cv::drawMatches(descs_1.GetImageMat(), descs_1.GetOpencvKeypoints(), descs_2.GetImageMat(),
		descs_2.GetOpencvKeypoints(), fl_good_matches, outImg4DrawMatches_flann);
	int img1_w = descs_1.GetImageWidth();
	for (int k = 0; k < 15; k++)
	{
		//1
		cv::putText(outImg4DrawMatches_flann, float2string(top100M_d1[99 - k]), cv::Point(10, 10 + (k * 20)), 1, 1, cv::Scalar(255, 0, 0), 2);
		cv::putText(outImg4DrawMatches_flann, float2string(top100_d1[99 - k]), cv::Point(100, 10 + (k * 20)), 1, 1, cv::Scalar(255, 0, 0), 2);
		//2
		cv::putText(outImg4DrawMatches_flann, float2string(top100M_d2[99 - k]), cv::Point(10 + img1_w, 10 + (k * 20)), 1, 1, cv::Scalar(255, 0, 0), 2);
		cv::putText(outImg4DrawMatches_flann, float2string(top100_d2[99 - k]), cv::Point(100 + img1_w, 10 + (k * 20)), 1, 1, cv::Scalar(255, 0, 0), 2);
	}

	cv::putText(outImg4DrawMatches_flann, uint2string(descs_1.GetNumOfDescriptors()), cv::Point(img1_w - 50, 10), 1, 1, cv::Scalar(0, 0, 2000), 2);
	cv::putText(outImg4DrawMatches_flann, uint2string(descs_2.GetNumOfDescriptors()), cv::Point(img1_w - 50, 25), 1, 1, cv::Scalar(0, 0, 2000), 2);
	cv::putText(outImg4DrawMatches_flann, int2string(fl_good_matches.size()), cv::Point(img1_w - 50, 50), 1, 1, cv::Scalar(100, 100, 2000), 2);
	cv::imshow("Draw matches flann", outImg4DrawMatches_flann);


	/*cv::Mat outImg4DrawMatches_bf;
	cv::drawMatches(descs_1.GetImageMat(), descs_1.GetOpencvKeypoints(), descs_2.GetImageMat(),
	descs_2.GetOpencvKeypoints(), bf_matches, outImg4DrawMatches_bf);
	cv::imshow("Draw matches bf", outImg4DrawMatches_bf);*/

	cv::waitKey();

	return outImg4DrawMatches;
}

void swapBlocks(std::vector<std::string> & outputLine, std::vector<float> & scoresPP, int m)
{
	std::string block[4];
	block[0] = outputLine[3 + (m * 4)];
	block[1] = outputLine[3 + (m * 4) + 1];
	block[2] = outputLine[3 + (m * 4) + 2];
	block[3] = outputLine[3 + (m * 4) + 3];

	outputLine.at(3 + (m * 4))     = outputLine[7 + (m * 4)];
	outputLine.at(3 + (m * 4) + 1) = outputLine[7 + (m * 4) + 1];
	outputLine.at(3 + (m * 4) + 2) = outputLine[7 + (m * 4) + 2];
	outputLine.at(3 + (m * 4) + 3) = outputLine[7 + (m * 4) + 3];

	outputLine.at(7 + (m * 4))     = block[0];
	outputLine.at(7 + (m * 4) + 1) = block[1];
	outputLine.at(7 + (m * 4) + 2) = block[2];
	outputLine.at(7 + (m * 4) + 3) = block[3];


	float temp		= scoresPP[m];
	scoresPP[m]		= scoresPP[m + 1];
	scoresPP[m + 1] = temp;
}

void blockSorting(std::vector<float> & scoresPP, std::vector<std::string> & outputLine)
{
	printf("\nSORTING:nok ___ ");
	for (int p = 1; p < QUERY_RETURN_SIZE; p++)
	{
		float max = scoresPP[QUERY_RETURN_SIZE - 1];
		for (int m = QUERY_RETURN_SIZE - 2; m >= p - 1; m--)
		{
			if (max>scoresPP[m])
				swapBlocks(outputLine, scoresPP, m);
			else
				max = scoresPP[m];
		}
	}
	printf("SORTING:ok \n");
}

void downsamplingAnalyzer(std::string imgPath, double size1, double size2, std::vector<DMatch >& good_matches)
{
	auto img1 = imread(imgPath.c_str());
	auto img2 = imread(imgPath.c_str());

	uchar_descriptors::resizeImage(&img1, size1);
	uchar_descriptors::resizeImage(&img2, size2);

	uchar_descriptors descriptor_1(imgPath.c_str(), img1,	"", AKAZE_FEATS);
	uchar_descriptors descriptor_2(	imgPath.c_str(), img2,	"", AKAZE_FEATS);

	descriptor_1.ExtractAKAZE();
	descriptor_2.ExtractAKAZE();

	findMatches(descriptor_1, descriptor_2, good_matches);
	/*******************************/
	cv::Mat OutImage4Matcher;
	cv::drawMatches(img1, descriptor_1.GetOpencvKeypoints(), img2,
		descriptor_2.GetOpencvKeypoints(), good_matches, OutImage4Matcher);

	cv::putText(OutImage4Matcher, uint2string(descriptor_1.GetNumOfDescriptors()), cv::Point(descriptor_1.GetImageWidth() - 50, 10), 1, 1, cv::Scalar(0, 0, 200), 2);
	cv::putText(OutImage4Matcher, uint2string(descriptor_2.GetNumOfDescriptors()), cv::Point(descriptor_1.GetImageWidth() - 50, 25), 1, 1, cv::Scalar(0, 0, 200), 2);
	cv::putText(OutImage4Matcher, int2string(good_matches.size()), cv::Point(descriptor_1.GetImageWidth() - 50, 50), 1, 1, cv::Scalar(100, 100, 200), 2);
	cv::imshow("Draw matches", OutImage4Matcher);
	waitKey();
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
