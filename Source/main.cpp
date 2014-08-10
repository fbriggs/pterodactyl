// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <iostream>

#include "wav_def.h"

#include "pParameters.h"

#include "pStringUtils.h"
#include "pFileUtils.h"
#include "pTimer.h"
#include "pTextFile.h"
#include "pVectorUtils.h"

#include "pImageBytes.h"
#include "pImageUtils.h"
#include "pFileIO_BMP.h"
#include "pFastGaussianBlur.h"
#include "pDrawPrimitives.h"

#include "pWavData.h"
#include "pSpectrogram.h"
#include "pNoiseReducer.h"

#include "pExample.h"
#include "pRandomForest_pthread.h"
#include "pKMeansCluster.h"
#include "pFarthestFirstTraversal.h"
#include "pMetaData_Rescale.h"

#include "pExample_MultiLabel.h"
#include "pMultiLabelErrorMeasure.h"
#include "pBinaryRelevance.h"
#include "pMultiLabelEnsembleOfClassifierChains.h"
#include "pRandomForest_OutOfBag.h"

#include "p2DRFSegmentation.h"
#include "p2DSegmentExtraction.h"
#include "p2DSegmentFeatures.h"
#include "pMetaData.h"

#include "pMultiLabelClassificationExperiment.h"

#include "pDatasetSpecific_HJAndrews.h"
#include "pDatasetSpecific_MLSPCompetition.h"
#include "pMIMLClassDiscoveryExperiment.h"

using namespace std;
using namespace multi_label_classification_experiment;

void setupFolders()
{
	system("mkdir src_wavs");
	system("mkdir chunks");
	system("mkdir spectrograms");
	system("mkdir filtered_spectrograms");
	system("mkdir segmentation_examples");
	system("mkdir segmentation_mask");
	system("mkdir segmentation_outlines");
	system("mkdir segmentation_probabilities");
	system("mkdir segments");
	system("mkdir segment_masks");
	system("mkdir segment_box_visualization");
	system("mkdir clustered_segments");
	system("mkdir segmentation_examples_support");		// after running through the process to get a set of segmentation examples, I was asked to give them to someone else. To make them useful, they need to be bundled with the WAV files they came from. those WAVs will be written here.
	system("mkdir unambiguous");
}


// splits wav files into smaller chunks, and does stereo to mono conversion at the same time 
void split(int durationSec, string srcFolder, string destFolder)
{
	cout << "splitting input waves into chunks of " << durationSec << " seconds" << endl;
	cout << "source folder: " << srcFolder << endl;
	cout << "dest folder: " << destFolder << endl;
	
	vector<string> sourceFilenames;
	pFileUtils::getFilesInDirectory(srcFolder, sourceFilenames);
	
	for(int i = 0; i < sourceFilenames.size(); ++i)
	{
		string filename = sourceFilenames[i];
		string prefix = pStringUtils::firstPartOfSplit(filename, '.');
		
		cout << filename << " " << i << " / " << sourceFilenames.size() << endl;
		
		pWavData wav(srcFolder+"/"+filename); 
		
		
		vector<pWavData> chunks;
		wav.splitIntoChunks(durationSec, chunks);
		
		for(int j = 0; j < chunks.size(); ++j)
		{
			string chunkFilename = destFolder + "/" + prefix +  "_" + pStringUtils::intToString(j) + ".wav";
			
			chunks[j].writeWAV(chunkFilename);
		}
	}
}

/// parallel version of split ///


pthread_mutex_t parallelSplitMutex = PTHREAD_MUTEX_INITIALIZER;

struct split_pthread_args
{
	int threadNum_;
	vector<string>* filenames_;
	int* currFile_; // this will keep track of the current file being processed
	PTIMER_MARK* parallelTimer_;
	string srcFolder_;
	string destFolder_;
	int durationSec_; // size of chunk to make
};

void* split_pthread_run(void* ptr)
{
	split_pthread_args& args = *((split_pthread_args*) ptr);
	
	while(true)
	{
		// look for an unclaimed file. this will be a locked operation (but its fast, so it shouldn't slow things down much)
		pthread_mutex_lock(&parallelSplitMutex);
		
		int fileToClaim = (*args.currFile_)++; // save the value, then increment it so the next thread gets the next file
		
		string filename = (*args.filenames_)[fileToClaim];
		string prefix = pStringUtils::firstPartOfSplit(filename, '.');
		float runtimeSec = pTimer::instance()->timeSinceMark(*args.parallelTimer_) / 1000.0f;
		cout << filename << " - " << fileToClaim << " / " << (*args.filenames_).size() << " total time = " << runtimeSec << " time/spectrogram=" << runtimeSec / float(fileToClaim) << " sec" << endl;
		
		pthread_mutex_unlock(&parallelSplitMutex);
		
		pWavData wav(args.srcFolder_+"/"+filename); 
		vector<pWavData> chunks;
		wav.splitIntoChunks(args.durationSec_, chunks);
		
		for(int j = 0; j < chunks.size(); ++j)
		{
			string chunkFilename = args.destFolder_ + "/" + prefix +  "_" + pStringUtils::intToString(j) + ".wav";
			
			chunks[j].writeWAV(chunkFilename);
		}

		if( (*args.currFile_) >= (*args.filenames_).size() ) return NULL;
	}
}

void split_multithread(int durationSec, string srcFolder, string destFolder)
{
	cout << "(parallel) splitting input waves into chunks of " << durationSec << " seconds" << endl;
	cout << "source folder: " << srcFolder << endl;
	cout << "dest folder: " << destFolder << endl;
	
	vector<string> sourceFilenames;
	pFileUtils::getFilesInDirectory(srcFolder, sourceFilenames);


	int numThreads = sysconf( _SC_NPROCESSORS_ONLN ); // set the number of threads to the number of processors
	cout << "# threads = " << numThreads << endl;
	
	vector<pthread_t> threads(numThreads);
	vector<split_pthread_args> threadArgs(numThreads);
	int currFileIndex = 0; // how far along we are
	
	PTIMER_MARK parallelTimer = pTimer::instance()->setMark(); // this timer simply tracks how long the whole parallel computation has run.
	
	for(int t = 0; t < numThreads; ++t)
	{
		split_pthread_args& args = threadArgs[t];
		
		args.threadNum_ = t;
		args.filenames_ = &sourceFilenames;
		args.currFile_ = &currFileIndex;
		args.parallelTimer_ = &parallelTimer;
		args.srcFolder_ = srcFolder;
		args.destFolder_ = destFolder;
		args.durationSec_ = durationSec;
		
		if( t!= 0 ) // the 0'th thread will be this thread
			pthread_create(&threads[t], NULL, split_pthread_run, (void*)&args);
	}
	
	cout << "starting work in main thread" << endl;
	
	// before waiting for other threads to finish, make some progress in this thread as well
	split_pthread_run(&threadArgs[0]);
	
	cout << "wating for threads to finish" << endl;
	
	for(int t = 0; t < numThreads; ++t)
	{
		cout << "waiting for thread " << t << endl;
		pthread_join(threads[t], NULL);
		cout << "thread " << t << " finished" << endl;
	}
	
	cout << "(all threads finished)" << endl;
}

void makeSpectrograms(string srcFolder, string destFolder)
{
	cout << "making spectrograms from wavs" << endl;
	cout << "source folder: " << srcFolder << endl;
	cout << "dest folder: " << destFolder << endl; 
	
	vector<string> sourceFilenames;
	pFileUtils::getFilesInDirectory(srcFolder, sourceFilenames);
	
	for(int i = 0; i < sourceFilenames.size(); ++i)
	{
		string filename = sourceFilenames[i];
		string prefix = pStringUtils::firstPartOfSplit(filename, '.');
		
		cout << filename << endl;
		
		pWavData wav(srcFolder+"/"+filename); 
		
		pSpectrogram spec(wav, pParameters::windowSize, pParameters::windowStep, pParameters::bandPassLow, pParameters::bandPassHigh);
		spec.normalizeSpectraTo01();
		pImageBytes img = spec.makeGreyscaleImageBoostContrast();
		pFileIO_BMP::write(img, destFolder + "/" + prefix + ".bmp");
	}
}

void noiseReductionFilter(string srcFolder, string destFolder)
{
	cout << "applying noise-reduction filters" << endl;
	cout << "source folder: " << srcFolder << endl;
	cout << "dest folder: " << destFolder << endl; 
	
	vector<string> sourceFilenames;
	pFileUtils::getFilesInDirectory(srcFolder, sourceFilenames);
	
	for(int i = 0; i < sourceFilenames.size(); ++i)
	{
		string filename = sourceFilenames[i];
		cout << filename << endl;
		
		pImageBytes spectrogram = pFileIO_BMP::read(srcFolder + "/" + filename);
		pImageBytes filteredSpectrogram = pNoiseReducer::whiteningFilter(spectrogram);
		//pImageBytes renormalized = pImageUtils::matrix2image(pImageUtils::image2matrix(filteredSpectrogram)); // this will force the image to be normalized to have a max value of 1
		pFileIO_BMP::write(filteredSpectrogram, destFolder + "/" + filename);
	}
}

/// parallel spectrogram and noise filter calculation (these are done simulatenously to reduce redundant disk io ///

pthread_mutex_t parallelSpectrogramMutex = PTHREAD_MUTEX_INITIALIZER;

struct spectrogram_pthread_args
{
	int threadNum_;
	vector<string>* filenames_;
	int* currFile_; // this will keep track of the current file being processed
	PTIMER_MARK* parallelTimer_;
	string srcFolder_;
	string destSpectrogramsFolder_;
	string destFilteredFolder_;
};

void* spectrogram_pthread_run(void* ptr)
{
	spectrogram_pthread_args& args = *((spectrogram_pthread_args*) ptr);
	
	while(true)
	{
		// look for an unclaimed file. this will be a locked operation (but its fast, so it shouldn't slow things down much)
		pthread_mutex_lock(&parallelSpectrogramMutex);
		
		int fileToClaim = (*args.currFile_)++; // save the value, then increment it so the next thread gets the next file
		string filename = (*args.filenames_)[fileToClaim];
		string prefix = pStringUtils::firstPartOfSplit(filename, '.');
		float runtimeSec = pTimer::instance()->timeSinceMark(*args.parallelTimer_) / 1000.0f;
		cout << filename << " - " << fileToClaim << " / " << (*args.filenames_).size() << " total time = " << runtimeSec << " time/spectrogram=" << runtimeSec / float(fileToClaim) << " sec" << endl;
		
		pthread_mutex_unlock(&parallelSpectrogramMutex);

		pWavData wav(args.srcFolder_+"/"+filename); 
		pSpectrogram spec(wav, pParameters::windowSize, pParameters::windowStep, pParameters::bandPassLow, pParameters::bandPassHigh);
		spec.normalizeSpectraTo01();
		pImageBytes spectrogram = spec.makeGreyscaleImageBoostContrast();

		pFileIO_BMP::write(spectrogram, args.destSpectrogramsFolder_ + "/" + prefix + ".bmp");
		
		if( args.destFilteredFolder_ != "" ) // this is to allow skipping the noise filter if desired
		{
			pImageBytes filteredSpectrogram = pNoiseReducer::whiteningFilter(spectrogram);
			pFileIO_BMP::write(filteredSpectrogram, args.destFilteredFolder_ + "/" + prefix + ".bmp");
		}
		
		if( (*args.currFile_) >= (*args.filenames_).size() ) return NULL;
	}
}

// note: pass the empty string for destFilteredFolder if you don't want to run the noise filter as well
void makeSpectrogramsAndApplyNoiseFilter_multithread(string srcWavsFolder, string destSpectrogramsFolder, string destFilteredFolder)
{
	cout << "(parallel) make spectrograms and (maybe) apply noise filter" << endl;
	cout << "source folder: " << srcWavsFolder << endl;
	cout << "dest spectrograms folder: " << destSpectrogramsFolder << endl;
	cout << "dest filtered folder: " << destFilteredFolder << endl;

	vector<string> sourceFilenames;
	pFileUtils::getFilesInDirectory(srcWavsFolder, sourceFilenames);


	int numThreads = sysconf( _SC_NPROCESSORS_ONLN ); // set the number of threads to the number of processors
	cout << "# threads = " << numThreads << endl;
	
	vector<pthread_t> threads(numThreads);
	vector<spectrogram_pthread_args> threadArgs(numThreads);
	int currFileIndex = 0; // how far along we are
	
	PTIMER_MARK parallelTimer = pTimer::instance()->setMark(); // this timer simply tracks how long the whole parallel computation has run.
	
	for(int t = 0; t < numThreads; ++t)
	{
		spectrogram_pthread_args& args = threadArgs[t];
		
		args.threadNum_ = t;
		args.filenames_ = &sourceFilenames;
		args.currFile_ = &currFileIndex;
		args.parallelTimer_ = &parallelTimer;
		args.srcFolder_ = srcWavsFolder;
		args.destSpectrogramsFolder_ = destSpectrogramsFolder;
		args.destFilteredFolder_ = destFilteredFolder;
		
		if( t!= 0 ) // the 0'th thread will be this thread
			pthread_create(&threads[t], NULL, spectrogram_pthread_run, (void*)&args);
	}
	
	cout << "starting work in main thread" << endl;
	
	// before waiting for other threads to finish, make some progress in this thread as well
	spectrogram_pthread_run(&threadArgs[0]);
	
	cout << "wating for threads to finish" << endl;
	
	for(int t = 0; t < numThreads; ++t)
	{
		cout << "waiting for thread " << t << endl;
		pthread_join(threads[t], NULL);
		cout << "thread " << t << " finished" << endl;
	}
	
	cout << "(all threads finished)" << endl;
}

void trainRFSegmentation(string segmentationExamplesFolder, string spectrogramsFolder, string outputRFModelFilename)
{
	cout << "training random forest segmentation model" << endl;
	cout << "segmentation examples folder folder: " << segmentationExamplesFolder << endl;
	cout << "spectrograms folder: " << spectrogramsFolder << endl; 
	cout << "output random forest model filename: " << outputRFModelFilename << endl;
	
	
	vector<string> exampleFilenames;
	pFileUtils::getFilesInDirectory(segmentationExamplesFolder, exampleFilenames);
	
	
	vector<pExample> exs; // training examples for the random forest
	
	for(int i = 0; i < exampleFilenames.size(); ++i)
	{
		string filename = exampleFilenames[i];
		cout << filename << " - " << i << " / " << exampleFilenames.size() << endl;
		
		pImageBytes annotation = pFileIO_BMP::read(segmentationExamplesFolder + "/" + filename);
		pImageBytes spectrogram = pFileIO_BMP::read(spectrogramsFolder + "/" + filename);
		
		assert(annotation.w_ == spectrogram.w_ && annotation.h_ == spectrogram.h_ && annotation.channels_ == spectrogram.channels_ && "annotation doesn't match spectrogram; probably one doesn't exist");
		
		p2DRFSegmentation::generateExamplesFromAnnotatedSpectrogram(annotation, spectrogram, exs);
		
		cout << "total training examples: " << exs.size() << endl;
	}
	
	cout << "feature dimension = " << exs[0].featureVector_.size() << endl;
	
	const int numClasses = 2;
	bool storeHistogramInLeaf = true;
	//const int maxTreeDepth = 20;
	const int maxTreeDepth = 1000000; // inifinity for all purposes
	
	pRandomForest_pthread rf;
	//rf.train(pParameters::segmentationNumRandomForestTrees, numClasses, exs, storeHistogramInLeaf, maxTreeDepth);
	rf.train_multithread(pParameters::segmentationNumRandomForestTrees, numClasses, exs, storeHistogramInLeaf, maxTreeDepth);
	rf.save(outputRFModelFilename);
}

void runRFSegmentation(string rfModelFilename, string srcFolder)
{
	cout << "running random forest segmentation" << endl;
	
	system("cd ./segmentation_probabilities; rm *");
	system("cd ./segmentation_mask; rm *");
	system("cd ./segmentation_outlines; rm *");
	
	pRandomForest_pthread rf;
	rf.load(rfModelFilename);
	
	vector<string> sourceFilenames;
	pFileUtils::getFilesInDirectory(srcFolder, sourceFilenames);
	
	for(int i = 0; i < sourceFilenames.size(); ++i)
	{
		PTIMER_MARK timerThisFile = pTimer::instance()->setMark();
		
		string filename = sourceFilenames[i];
		cout << filename << " - " << i << " / " << sourceFilenames.size() << endl;
		
		pImageBytes spectrogram = pFileIO_BMP::read(srcFolder + "/" + filename);
		pImageBytes probMap = p2DRFSegmentation::segmentationProbabilityMap(spectrogram, rf);
		
		pImageBytes blurred = pFastGaussianBlur::blur(probMap, pParameters::probmapBlurRadius);
		pImageBytes mask = pImageUtils::threshold(blurred, pParameters::segmentationProbabilityThreshold);
		
		pImageBytes outlines = pImageUtils::colorize(spectrogram);
		pImageUtils::drawMaskOutline(outlines, mask, pRGB(1, 0, 0));
		
		pFileIO_BMP::write(probMap,		"segmentation_probabilities/" + filename);
		pFileIO_BMP::write(mask,		"segmentation_mask/" + filename);
		pFileIO_BMP::write(outlines,	"segmentation_outlines/" + filename);
		
		cout << "time to segment this file = " << pTimer::instance()->timeSinceMark(timerThisFile) / 1000.0f << "sec" << endl;
	}
}

/// stuff related to running segmentation in parallel.  ///

pthread_mutex_t parallelSegmentationMutex = PTHREAD_MUTEX_INITIALIZER;

struct pSegmentation_pthread_args
{
	int threadNum_;
	vector<string>* filenames_;
	pRandomForest_pthread* rf_;
	//vector<bool>* fileIsClaimed_;
	
	int* currFile_; // this will keep track of the current file being processed
	
	PTIMER_MARK* parallelTimer_;
	string srcFolder_;
};

void* runRFSegmentation_pthread_run(void* ptr)
{
	pSegmentation_pthread_args& args = *((pSegmentation_pthread_args*) ptr);
	
	while(true)
	{
		// look for an unclaimed file. this will be a locked operation (but its fast, so it shouldn't slow things down much)
		pthread_mutex_lock(&parallelSegmentationMutex);
		
		int fileToClaim = (*args.currFile_)++; // save the value, then increment it so the next thread gets the next file
		
		
		
		string filename = (*args.filenames_)[fileToClaim];
		float runtimeSec = pTimer::instance()->timeSinceMark(*args.parallelTimer_) / 1000.0f;
		cout << filename << " - " << fileToClaim << " / " << (*args.filenames_).size() << " total time = " << runtimeSec << " time/spectrogram=" << runtimeSec / float(fileToClaim) << " sec" << endl;
		
		pthread_mutex_unlock(&parallelSegmentationMutex);
		
		pImageBytes spectrogram = pFileIO_BMP::read(args.srcFolder_ + "/" + filename);
		pImageBytes probMap = p2DRFSegmentation::segmentationProbabilityMap(spectrogram, *args.rf_);
		
		pImageBytes blurred = pFastGaussianBlur::blur(probMap, pParameters::probmapBlurRadius);
		pImageBytes mask = pImageUtils::threshold(blurred, pParameters::segmentationProbabilityThreshold);
		
		pImageBytes outlines = pImageUtils::colorize(spectrogram);
		pImageUtils::drawMaskOutline(outlines, mask, pRGB(1, 0, 0));
		
		pFileIO_BMP::write(probMap,		"segmentation_probabilities/" + filename);
		pFileIO_BMP::write(mask,		"segmentation_mask/" + filename);
		pFileIO_BMP::write(outlines,	"segmentation_outlines/" + filename);
		
		if( (*args.currFile_) >= (*args.filenames_).size() ) return NULL;
	}
}

// same as runRFSegmentation, but does it multi-threaded
void runRFSegmentation_multithread(string rfModelFilename, string srcFolder)
{
	cout << "running random forest segmentation (in parallel)" << endl;
	
	system("cd ./segmentation_probabilities; rm *");
	system("cd ./segmentation_mask; rm *");
	system("cd ./segmentation_outlines; rm *");
	
	pRandomForest_pthread rf;
	rf.load(rfModelFilename);
	
	vector<string> sourceFilenames;
	pFileUtils::getFilesInDirectory(srcFolder, sourceFilenames);
	
	int numThreads = sysconf( _SC_NPROCESSORS_ONLN ); // set the number of threads to the number of processors
	assert(numThreads >= 1);
	cout << "# threads = " << numThreads << endl;
	
	vector<pthread_t> threads(numThreads);
	vector<pSegmentation_pthread_args> threadArgs(numThreads);
	int currFileIndex = 0; // how far along we are
	
	PTIMER_MARK parallelTimer = pTimer::instance()->setMark(); // this timer simply tracks how long the whole parallel computation has run.
	
	for(int t = 0; t < numThreads; ++t)
	{
		pSegmentation_pthread_args& args = threadArgs[t];
		
		args.threadNum_ = t;
		args.filenames_ = &sourceFilenames;
		args.rf_ = &rf;
		args.currFile_ = &currFileIndex;
		args.parallelTimer_ = &parallelTimer;
		args.srcFolder_ = srcFolder;
		
		if( t!= 0 ) // the 0'th thread will be this thread
			pthread_create(&threads[t], NULL, runRFSegmentation_pthread_run, (void*)&args);
	}
	
	cout << "starting work in main thread" << endl;
	
	// before waiting for other threads to finish, make some progress in this thread as well
	runRFSegmentation_pthread_run(&threadArgs[0]);
	
	cout << "wating for threads to finish" << endl;
	
	for(int t = 0; t < numThreads; ++t)
	{
		cout << "waiting for thread " << t << endl;
		pthread_join(threads[t], NULL);
		cout << "thread " << t << " finished" << endl;
	}
	
	
	cout << "(all threads finished)" << endl;
}

/// end of multi-threaded segmentation ///

void processSegments(string spectrogramsFolder, string masksFolder, string inputBagIdFilename)
{
	PTIMER_MARK processSegmentsTimer = pTimer::instance()->setMark();

	system("cd segments; rm *");
	system("cd segment_masks; rm *");
	
	bool bagIdSupplied = (inputBagIdFilename != ""); // if a non-empty string is passed for inputBagIdFilename, it will be loaded and used, rather than generated as we go
	
	if( bagIdSupplied ) cout << "analyzing segments with pre-supplied bag_id file" << endl;
	
	pTextFile fBagID(bagIdSupplied ? inputBagIdFilename : "bag_id2filename.txt", bagIdSupplied ? PFILE_READ : PFILE_WRITE);
	pTextFile fSegmentFeatures("segment_features.txt", PFILE_WRITE);
	pTextFile fSegmentRectangles("segment_rectangles.txt", PFILE_WRITE);
	
	if( !bagIdSupplied) fBagID << "bag_id,filename" << "\n";
	fSegmentFeatures << "bag_id,segment_id,[feature vector]" << "\n";
	fSegmentRectangles << "bag_id,segment_id,min_x,max_x,min_y,max_y" << "\n";
	
	// either get the source filenames from the supplied bag_id file, or read the contents of the masksFolder and use those filenames instead
	vector<string> srcFilenames;
	if( bagIdSupplied ) // parse the supplied bag_id file. this will give an order to process the files in
	{
		assert(fBagID.readLine() == "bag_id,filename"); // header check
		
		while(!fBagID.eof())
		{
			vector<string> parts = pStringUtils::split(fBagID.readLine(), ",");
			int bag_id = pStringUtils::stringToInt(parts[0]);
			string filename = parts[1];
			
			srcFilenames.push_back(filename);
			
			if( !(srcFilenames.size()-1 == bag_id) )
			{
				cout << bag_id << " " << srcFilenames.size() -1 << " " << filename << endl;
			}
			
			assert(srcFilenames.size()-1 == bag_id);
		}
	}
	else 
		pFileUtils::getFilesInDirectory(masksFolder, srcFilenames);
		
	// iterate over all files in the masks folder (and also load the corresponding spectrogram)
	for(int i = 0; i < srcFilenames.size(); ++i)
	{
		string filename = srcFilenames[i];
		string prefix = pStringUtils::firstPartOfSplit(filename, '.');
		
		if( !bagIdSupplied) fBagID << i << "," << filename << "\n";
		
		pImageBytes spectrogram = pFileIO_BMP::read(spectrogramsFolder + "/" + filename);
		pImageBytes mask = pFileIO_BMP::read(masksFolder + "/" + filename);
		
		vector<pImageBytes> croppedSegments;
		vector<pImageBytes> croppedMasks;
		vector<pRectangle> segmentRects;	// the rectangle coodinates of each segment are stored here; may be useful for those that prefer to treat segments as rectangles instead of a binary mask
		p2DSegmentExtraction::extractSegments(spectrogram, mask, croppedSegments, croppedMasks, segmentRects);
		
		float runtimeSec = pTimer::instance()->timeSinceMark(processSegmentsTimer) / 1000.0f;
		
		cout << filename << ", " << i << " / " << srcFilenames.size() << " #segments=" << croppedSegments.size() << " total time = " << runtimeSec << " sec, time / spectrogram =  " << runtimeSec / float(i) << " sec" << endl;
		
		// draw rectangles for each segment on the spectrogram as one more kind of visualization
		//pImageBytes spectrogramWithBoxes = pImageUtils::colorize(spectrogram);
		
		
		// now we have all the segments. the next step is to further process each one individually
		for(int j = 0; j < croppedSegments.size(); ++j)
		{
			pImageBytes packed = p2DSegmentExtraction::packVertical(croppedSegments[j], croppedMasks[j]);
			pFileIO_BMP::write(packed,		"segments/"			+ pStringUtils::intToString(i) + "," + pStringUtils::intToString(j) + ".bmp");
			
			// to avoid generating unecessary extra data, packed masks are not generated/saved, but this can be done be un-commenting the code below
			//pImageBytes packedMask = p2DSegmentExtraction::packVertical(croppedMasks[j], croppedMasks[j]);
			//pFileIO_BMP::write(packedMask,	"segment_masks/"	+ pStringUtils::intToString(i) + "," + pStringUtils::intToString(j) + ".bmp");
			
			vector<float> profile_features	= p2DSegmentFeatures::profileFeatures(croppedSegments[j], croppedMasks[j]);
			vector<float> shape_features	= p2DSegmentFeatures::maskShapeFeatures(croppedMasks[j], segmentRects[j]);
			vector<float> hog_features		= p2DSegmentFeatures::hogFeatures(croppedSegments[j], croppedMasks[j]);
			
			vector<float> all_features;
			pVectorUtils::append(profile_features,	all_features);
			pVectorUtils::append(shape_features,	all_features);
			pVectorUtils::append(hog_features,		all_features);
			
			fSegmentFeatures << i << "," << j;
			for(int k = 0; k < all_features.size(); ++k)
				fSegmentFeatures << "," << all_features[k];
			fSegmentFeatures << "\n";
			
			fSegmentRectangles << i << "," << j << "," << segmentRects[j].minX_ << "," << segmentRects[j].maxX_ << "," << segmentRects[j].minY_ << "," << segmentRects[j].maxY_ << "," << "\n";
			
			//pDrawPrimitives::boxOutline(spectrogramWithBoxes, segmentRects[j].minX_, segmentRects[j].maxX_, segmentRects[j].minY_, segmentRects[j].maxY_, pRGB(1,0,0));
		}
		
		//pFileIO_BMP::write(spectrogramWithBoxes, "segment_box_visualization/"+ filename);
	}
	
	fBagID.close();
	fSegmentFeatures.close();
	fSegmentRectangles.close();
}

void generateSegmentMosaic()
{	
	vector<string> segFilenames;
	pFileUtils::getFilesInDirectory("segments", segFilenames);
	
	vector<pImageBytes*> segImages;
	
	for(int i = 0; i < segFilenames.size(); ++i)
	{
		string filename = segFilenames[i];
		cout << filename << ", " << i << " / " << segFilenames.size() << endl;
		segImages.push_back(new pImageBytes(pFileIO_BMP::read("segments/" + filename)));
	}
	
	cout << "building mosaic" << endl;
	
	pImageBytes mosaic = pImageUtils::irregularImageMosaic(segImages, pParameters::segmentMosaicImageWidth);
	
	pFileIO_BMP::write(mosaic, "segment_mosaic.bmp");

	for(int i = 0; i < segImages.size(); ++i) delete segImages[i];
}

void clusterSegments(int k)
{
	system("cd clustered_segments; rm *");

	vector<pAudioFile_MetaData> metaData;
	pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
	pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
	
	cout << "rescaling features" << endl;
	const int featureDim = pMetaData_Parser::getSegmentFeatureDim(metaData);
	pMetaData_Rescale::rescaleTo01(metaData, featureDim);
	
	// put the data into a vector< vector<float>> because that is what the kmeans implementation wants, then apply some feature rescaling
	vector< vector<float> > fvs;
	for(int i = 0; i < metaData.size(); ++i)
		for(int j = 0; j < metaData[i].segments_.size(); ++j)
			fvs.push_back(metaData[i].segments_[j].featureVector_);
	
	//const int k = pParameters::clusterSegmentsVisualizationK;
	
	cout << "running k-means++ with k = " << k << " d = " << featureDim << " n = " << fvs.size() << endl;
	pKMeansCluster km(featureDim, k, fvs.size());
	km.randomizeKMeansPlusPlus(fvs);
	km.runToConvergenceOrMaxIterations(fvs, pParameters::kMeansMaxIterations);
	
	vector< vector<pImageBytes*> > segmentsInCluster(k);
	
	pTextFile clusterFile("segment2cluster.txt", PFILE_WRITE);
	clusterFile << "bag_id,segment_id,cluster_assignment" << "\n";
	for(int i = 0; i < metaData.size(); ++i)
	{
		// "clusterfy" means assign each segment to a cluster (the k-means clustering acts like a classifier in that it maps a feature vector to a cluster label).
		cout << "clusterfying " << metaData[i].filename_ << " " << i << " / " << metaData.size() << endl;
		 
		for(int j = 0; j < metaData[i].segments_.size(); ++j)
		{
			// load the corresponding segment
			pImageBytes* segmentImg = new pImageBytes(pImageUtils::greyscale(pFileIO_BMP::read("segments/" + pStringUtils::intToString(i) + "," + pStringUtils::intToString(j) + ".bmp")));
			
			int clusterLabel = km.getLabelForExample(metaData[i].segments_[j].featureVector_);
			
			clusterFile << i << "," << j << "," << clusterLabel << "\n";
			
			segmentsInCluster[clusterLabel].push_back(segmentImg);
		}	
	}		
	clusterFile.close();
	
	pImageBytes paddingImg(pParameters::segmentMosaicImageWidth, 1, 1);
	pImageUtils::fillChannel(paddingImg, 0, 0.1);
	
	vector<pImageBytes> allClusterImages;
	
	for(int c = 0; c < k; ++c)
	{
		if( segmentsInCluster[c].size() != 0 )
		{
			cout << "building mosaic for cluster " << c << endl;
			pImageBytes mosaic = pImageUtils::irregularImageMosaic(segmentsInCluster[c], pParameters::segmentMosaicImageWidth);
			pFileIO_BMP::write(mosaic, "clustered_segments/cluster" + pStringUtils::intToString(c) + ".bmp");
			
			allClusterImages.push_back(mosaic);
			allClusterImages.push_back(paddingImg);
		}
		
		for(int j = 0; j < segmentsInCluster[c].size(); ++j)
			delete segmentsInCluster[c][j];
	}
	
	pImageBytes stackedClustersImage = pImageUtils::stackVertical(allClusterImages);
	pFileIO_BMP::write(stackedClustersImage, "segment_clusters.bmp");
}


///////////////// species discovery methods ////////////////////

// the idea of these functions for species discovery is to a select a diverse/interesting subset of the data for further inspection.
// we would like to discover as many species as possible for a fixed amount of human effort.

/// one approach is to traverse segments. these few functions do that ///

// select k random segments
void speciesDiscovery_segment_random(int k)
{
	vector<pAudioFile_MetaData> metaData;
	pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
	pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
	
	vector< pair<int, int> > segmentLookup;
	for(int i = 0; i < metaData.size(); ++i)
		for(int j = 0; j < metaData[i].segments_.size(); ++j)
			segmentLookup.push_back(pair<int, int>(i, j));
	
	pRandom::templatedShuffle(segmentLookup);
	
	vector<pImageBytes*> selectedSegmentImages;
	for(int i = 0; i < k; ++i)
	{
		int bag_id = segmentLookup[i].first;
		int segment_id = segmentLookup[i].second;
		pImageBytes* segmentImg = new pImageBytes(pImageUtils::greyscale(pFileIO_BMP::read("segments/" + pStringUtils::intToString(bag_id) + "," + pStringUtils::intToString(segment_id) + ".bmp")));
		selectedSegmentImages.push_back(segmentImg);
	}
	
	pImageBytes mosaic = pImageUtils::irregularImageMosaic(selectedSegmentImages, pParameters::segmentMosaicImageWidth);
	pFileIO_BMP::write(mosaic, "speciesDiscovery_segment_random.bmp");
	
	for(int i = 0; i < selectedSegmentImages.size(); ++i)
		delete selectedSegmentImages[i];
}

// use farthest-first traversal to select k segments (for species discovery).
// the purpose of this is to discover as many species as possible while examining only k segments
void speciesDiscovery_segment_farthest(int k)
{
	vector<pAudioFile_MetaData> metaData;
	pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
	pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
	
	cout << "rescaling features" << endl;
	const int featureDim = pMetaData_Parser::getSegmentFeatureDim(metaData);
	pMetaData_Rescale::rescaleTo01(metaData, featureDim);
	
	// put the data into a vector< vector<float>> because that is what the kmeans implementation wants, then apply some feature rescaling
	vector< vector<float> > fvs;
	vector< pair<int, int> > segmentLookup;
	for(int i = 0; i < metaData.size(); ++i)
	{
		for(int j = 0; j < metaData[i].segments_.size(); ++j)
		{
			fvs.push_back(metaData[i].segments_[j].featureVector_);
			segmentLookup.push_back(pair<int, int>(i, j));
		}
	}
	
	vector<int> traversal = pFarthestFirstTraversal::traverseMaxMin(fvs, k);
	
	vector<pImageBytes*> selectedSegmentImages;
	for(int i = 0; i < k; ++i)
	{
		int bag_id = segmentLookup[traversal[i]].first;
		int segment_id = segmentLookup[traversal[i]].second;
		pImageBytes* segmentImg = new pImageBytes(pImageUtils::greyscale(pFileIO_BMP::read("segments/" + pStringUtils::intToString(bag_id) + "," + pStringUtils::intToString(segment_id) + ".bmp")));
		selectedSegmentImages.push_back(segmentImg);
	}
	
	pImageBytes mosaic = pImageUtils::irregularImageMosaic(selectedSegmentImages, pParameters::segmentMosaicImageWidth);
	pFileIO_BMP::write(mosaic, "speciesDiscovery_segment_farthest.bmp");
	
	for(int i = 0; i < selectedSegmentImages.size(); ++i)
		delete selectedSegmentImages[i];
}

// apply k-means clustering and select one example from each cluster.
void speciesDiscovery_segment_cluster(int k)
{
	vector<pAudioFile_MetaData> metaData;
	pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
	pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
	
	cout << "rescaling features" << endl;
	const int featureDim = pMetaData_Parser::getSegmentFeatureDim(metaData);
	pMetaData_Rescale::rescaleTo01(metaData, featureDim);
		
	// put the data into a vector< vector<float>> because that is what the kmeans implementation wants, then apply some feature rescaling
	vector< vector<float> > fvs;
	for(int i = 0; i < metaData.size(); ++i)
		for(int j = 0; j < metaData[i].segments_.size(); ++j)
			fvs.push_back(metaData[i].segments_[j].featureVector_);
	
	int dim = fvs[0].size();
	cout << "running k-means++ with k = " << k << " d = " << dim << " n = " << fvs.size() << endl;
	pKMeansCluster km(dim, k, fvs.size());
	km.randomizeKMeansPlusPlus(fvs);
	km.runToConvergenceOrMaxIterations(fvs, pParameters::kMeansMaxIterations);
	
	// map each segment to a cluster
	vector< vector< pair<int, int> > > segmentIndicesInCluster(k);
	for(int i = 0; i < metaData.size(); ++i)
	{
		// "clusterfy" means assign each segment to a cluster (the k-means clustering acts like a classifier in that it maps a feature vector to a cluster label).
		cout << "clusterfying " << metaData[i].filename_ << " " << i << " / " << metaData.size() << endl;
		for(int j = 0; j < metaData[i].segments_.size(); ++j)
		{
			int clusterLabel = km.getLabelForExample(metaData[i].segments_[j].featureVector_);
			segmentIndicesInCluster[clusterLabel].push_back(pair<int, int>(i,j));
		}	
	}		
	
	// find the closest example to each cluster center
	vector<pImageBytes*> selectedSegmentImages;
	for(int cluster = 0; cluster < k; ++cluster)
	{
		pair<int, int> closestSegment;
		float closestDist = FLT_MAX;
		
		for(int i = 0; i < segmentIndicesInCluster[cluster].size(); ++i)
		{
			int bag_id = segmentIndicesInCluster[cluster][i].first;
			int segment_id = segmentIndicesInCluster[cluster][i].second;
			
			float dist = km.distanceToCenter(cluster, metaData[bag_id].segments_[segment_id].featureVector_);
			
			if( dist < closestDist )
			{
				closestDist = dist;
				closestSegment = segmentIndicesInCluster[cluster][i];
			}
		}
		
		pImageBytes* segmentImg = new pImageBytes(pImageUtils::greyscale(pFileIO_BMP::read("segments/" + pStringUtils::intToString(closestSegment.first) + "," + pStringUtils::intToString(closestSegment.second) + ".bmp")));
		selectedSegmentImages.push_back(segmentImg);
	}
	
	pImageBytes mosaic = pImageUtils::irregularImageMosaic(selectedSegmentImages, pParameters::segmentMosaicImageWidth);
	pFileIO_BMP::write(mosaic, "speciesDiscovery_segment_cluster.bmp");
	
	for(int i = 0; i < selectedSegmentImages.size(); ++i)
		delete selectedSegmentImages[i];
}

/// another approach is to traverse chunks/recordings. these next functions do that ///

void speciesDiscovery_chunk_random(int k)
{
	// load the dataset	
	vector<pAudioFile_MetaData> metaData;
	pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
	pMetaData_Parser::parseBagLevelFeatures("histogram_of_segments.txt", metaData);
	
	pRandom::templatedShuffle(metaData);
	
	vector<pImageBytes> selectedSpectrograms;
	for(int i = 0; i < k; ++i)
		selectedSpectrograms.push_back(pImageUtils::greyscale(pFileIO_BMP::read("spectrograms/" + metaData[i].filename_)));
	
	pImageBytes stackedImages = pImageUtils::stackVertical(selectedSpectrograms);
	pFileIO_BMP::write(stackedImages, "speciesDiscovery_chunk_random.bmp");
}

void speciesDiscovery_chunk_farthest(int k)
{
		// load the dataset	
	vector<pAudioFile_MetaData> metaData;
	pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
	pMetaData_Parser::parseBagLevelFeatures("histogram_of_segments.txt", metaData);

	// move the histogram of segments features a separate vector for clustering
	vector< vector<float> > fvs;
	for(int i = 0; i < metaData.size(); ++i)
		fvs.push_back(metaData[i].bagLevelFeatureVector_);
		
	vector<int> traversal = pFarthestFirstTraversal::traverseMaxAvgDist(fvs, k);
	
	vector<pImageBytes> selectedSpectrograms;
	for(int i = 0; i < traversal.size(); ++i)
		selectedSpectrograms.push_back(pImageUtils::greyscale(pFileIO_BMP::read("spectrograms/" + metaData[traversal[i]].filename_)));
	
	pImageBytes stackedImages = pImageUtils::stackVertical(selectedSpectrograms);
	pFileIO_BMP::write(stackedImages, "speciesDiscovery_chunk_farthest.bmp");
}

void speciesDiscovery_chunk_cluster(int k)
{
	// load the dataset	
	vector<pAudioFile_MetaData> metaData;
	pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
	pMetaData_Parser::parseBagLevelFeatures("histogram_of_segments.txt", metaData);

	// move the histogram of segments features a separate vector for clustering
	vector< vector<float> > fvs;
	for(int i = 0; i < metaData.size(); ++i)
		fvs.push_back(metaData[i].bagLevelFeatureVector_);
	
	int dim = fvs[0].size();
	cout << "running k-means++ with k = " << k << " d = " << dim << " n = " << fvs.size() << endl;
	pKMeansCluster km(dim, k, fvs.size());
	km.randomizeKMeansPlusPlus(fvs);
	km.runToConvergenceOrMaxIterations(fvs, pParameters::kMeansMaxIterations);
	
	// map each chunk to a cluster
	vector< vector<int> > chunksInCluster(k);
	for(int i = 0; i < metaData.size(); ++i)
	{
		int clusterLabel = km.getLabelForExample(metaData[i].bagLevelFeatureVector_);
		chunksInCluster[clusterLabel].push_back(i);
	}
	
	// find the closest chunk to each cluster center
	vector<pImageBytes> selectedSpectrograms;
	for(int i = 0; i < k; ++i)
	{
		int closestIndex;
		float closestDist = FLT_MAX;
		
		for(int j = 0; j < chunksInCluster[i].size(); ++j)
		{
			float dist = km.distanceToCenter(i, metaData[chunksInCluster[i][j]].bagLevelFeatureVector_);
			
			if( dist < closestDist )
			{
				closestDist = dist;
				closestIndex = chunksInCluster[i][j];
			}
		}
		selectedSpectrograms.push_back(pImageUtils::greyscale(pFileIO_BMP::read("spectrograms/" + metaData[closestIndex].filename_)));
	}
	
	pImageBytes stackedImages = pImageUtils::stackVertical(selectedSpectrograms);
	pFileIO_BMP::write(stackedImages, "speciesDiscovery_chunk_cluster.bmp");
}

///////////////// end of species discovery methods ////////////////////

// make "histogram of segment" (HOS) bag-level feature vectors from segmetn feature vectors 
void histogramOfSegmentFeatures()
{
	vector<pAudioFile_MetaData> metaData;
	pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
	pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
	
	cout << "rescaling features" << endl;
	const int featureDim = pMetaData_Parser::getSegmentFeatureDim(metaData);
	pMetaData_Rescale::rescaleTo01(metaData, featureDim);
	
	// put the data into a vector< vector<float>> because that is what the kmeans implementation wants
	vector< vector<float> > fvs;
	for(int i = 0; i < metaData.size(); ++i)
		for(int j = 0; j < metaData[i].segments_.size(); ++j)
			fvs.push_back(metaData[i].segments_[j].featureVector_);

	// cluster the segments
	const int k = pParameters::histogramOfSegmentsClusterK;
	int dim = fvs[0].size();
	cout << "running k-means++ with k = " << k << " d = " << dim << " n = " << fvs.size() << endl;
	pKMeansCluster km(dim, k, fvs.size());
	km.randomizeKMeansPlusPlus(fvs);
	km.runToConvergenceOrMaxIterations(fvs, pParameters::kMeansMaxIterations);
	
	pTextFile fBagFeatures("histogram_of_segments.txt", PFILE_WRITE);
	fBagFeatures << "bag_id,[histogram of segment features]" << "\n";
	
	// for each bag/recording, build a histogram of segments using the clustering as a codebook
	for(int i = 0; i < metaData.size(); ++i)
	{
		
		// histogram feature
		vector<float> histo(k);
		for(int j = 0; j < metaData[i].segments_.size(); ++j)
		{
			int segmentCluster = km.getLabelForExample(metaData[i].segments_[j].featureVector_);
			histo[segmentCluster] += 1;
		}
		pVectorUtils::normalizeVectorToPDF(histo); // comment this out for raw-count
		
		fBagFeatures << i;
		for(int j = 0; j < histo.size(); ++j) fBagFeatures << "," << histo[j];
		
		fBagFeatures << "\n";
	}
	
	fBagFeatures.close();
}

// histogram of segments is one kind of bag-level feature. however, we can construct additional
// bag-level features, e.g. the # of segments, and the mean/stdev of each instance-level feature.
// this function assumes that the histogram of segments features have already been computed, because 
// they will be included in the final feature vector.
void makeAdditionalBagLevelFeatures()
{
	vector<pAudioFile_MetaData> metaData;
	pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
	pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
	pMetaData_Parser::parseBagLevelFeatures("histogram_of_segments.txt", metaData);
	
	cout << "rescaling features" << endl;
	const int segmentFeatureDim = pMetaData_Parser::getSegmentFeatureDim(metaData);
	pMetaData_Rescale::rescaleTo01(metaData, segmentFeatureDim);
	
	for(int i = 0; i < metaData.size(); ++i)
	{
		if( i%1000 == 0 ) cout << "*"; cout.flush();
		
		// compute the average of the segment features
		vector<float> meanFeatures(segmentFeatureDim, 0);
		for(int j = 0; j < metaData[i].segments_.size(); ++j)
			for(int k = 0; k < segmentFeatureDim; ++k)
				meanFeatures[k] += metaData[i].segments_[j].featureVector_[k];
		
		if( metaData[i].segments_.size() != 0 )
			for(int k = 0; k < segmentFeatureDim; ++k)
				meanFeatures[k] /= float(metaData[i].segments_.size());
		
		// compute the stdev of the segment features
		vector<float> stdevFeatures(segmentFeatureDim, 0);
		for(int j = 0; j < metaData[i].segments_.size(); ++j)
			for(int k = 0; k < segmentFeatureDim; ++k)
				stdevFeatures[k] += pow(meanFeatures[k] - metaData[i].segments_[j].featureVector_[k], 2);
		
		if( metaData[i].segments_.size() != 0 )
			for(int k = 0; k < segmentFeatureDim; ++k)
				stdevFeatures[k] = sqrt(stdevFeatures[k] / float(metaData[i].segments_.size())); 
		
		// append all of the new features to the existing HOS features
		metaData[i].bagLevelFeatureVector_.push_back(metaData[i].segments_.size()); // first new feature: # of segments 
		
		for(int k = 0; k < segmentFeatureDim; ++k)
		{
			metaData[i].bagLevelFeatureVector_.push_back(meanFeatures[k]);
			metaData[i].bagLevelFeatureVector_.push_back(stdevFeatures[k]);
		}
	}
	
	pTextFile fBagFeatures("additional_bag_features.txt", PFILE_WRITE);
	fBagFeatures << "bag_id,[features]" << "\n";
	
	// for each bag/recording, build a histogram of segments using the clustering as a codebook
	for(int i = 0; i < metaData.size(); ++i)
	{
		fBagFeatures << i;
		for(int j = 0; j < metaData[i].bagLevelFeatureVector_.size(); ++j) fBagFeatures << "," << metaData[i].bagLevelFeatureVector_[j];
		fBagFeatures << "\n";
	}
	
	fBagFeatures.close();
}

// this function classifies each spectrogram as rain or non-rain. several files and folders are assumed to exist already:
// folder: /examples_rain_negative
// folder: /examples_rain_positive
// ^-- these folders should contain spectrograms for positive and negative examples of rain. note: the spectrograms will not be used directly.
// all that matters is the filenames in these folders, which are used to get a list of positive and negative examples.
// it is very fast to set these folders up by following this procedure: 
// 1. make the folders /examples_rain_negative and /examples_rain_positive
// 2. run copy_random_subset_reseed to copy some # of spectrograms (e.g. 1000) to /examples_rain_negative
// 3. look at all of the files in /examples_rain_negative, and move any that look like rain to /examples_rain_positive
//
// file: additional_bag_features.txt - these are the features generated by makeAdditionalBagLevelFeatures() 
//
// results generated: 
// a folder /probably_rain containing the spectrograms which got a rain probability > 0.5 (this is just for visualization purposes)
// a file rain_probabilities.txt which gives the probability for each recording to be rain
void runRainClassifier()
{
	vector<pAudioFile_MetaData> metaData;
	pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
	pMetaData_Parser::parseBagLevelFeatures("additional_bag_features.txt", metaData);
	
	map<string, int> filenameToIndex;
	for(int i = 0; i < metaData.size(); ++i)
		filenameToIndex[metaData[i].filename_] = i;
	
	vector<string> positiveFilenames;
	vector<string> negativeFilenames;
	pFileUtils::getFilesInDirectory("examples_rain_positive", positiveFilenames);
	pFileUtils::getFilesInDirectory("examples_rain_negative", negativeFilenames);
	
	vector<pExample> binaryExamples;
	
	for(int i = 0; i < positiveFilenames.size(); ++i)
		binaryExamples.push_back(pExample(metaData[filenameToIndex[positiveFilenames[i]]].bagLevelFeatureVector_, 1));
		
	for(int i = 0; i < negativeFilenames.size(); ++i)
		binaryExamples.push_back(pExample(metaData[filenameToIndex[negativeFilenames[i]]].bagLevelFeatureVector_, 0));
		
	cout << "training rain classifier" << endl;
	pRandomForest_pthread rf;
	const int numTrees = 1000;
	rf.train_multithread(numTrees, 2, binaryExamples);
	
	cout << "applying rain classifier" << endl;
	
	pTextFile resultFile("rain_probabilities.txt", PFILE_WRITE);
	
	resultFile << "bag_id,rain_probability" << "\n";
	
	vector< pair<float, string> > rainProbAndFilename;
	for(int i = 0; i < metaData.size(); ++i)
	{
		if( i % 1000 == 0 ) cout << "*"; cout.flush();
		
		float rainProbability = rf.estimateClassProbabilities(metaData[i].bagLevelFeatureVector_)[1];
		
		rainProbAndFilename.push_back(pair<float, string>(rainProbability, metaData[i].filename_));
		
		resultFile << metaData[i].bag_id_ << "," << rainProbability << "\n";
	}
	resultFile.close();
	
	system("mkdir ./probably_rain");
	system("cd ./probably_rain; rm *");
	cout << "\n----\n" << endl;
	
	sort(rainProbAndFilename.rbegin(), rainProbAndFilename.rend());
	
	for(int i = 0; i < rainProbAndFilename.size(); ++i)
	{
		cout << rainProbAndFilename[i].first << "\t" << rainProbAndFilename[i].second << endl;
		
		if( rainProbAndFilename[i].first >= 0.5 )
			system(string("cp spectrograms/" + rainProbAndFilename[i].second + " probably_rain").c_str());
	}
}

// generate a single text file with everything needed to run a multi-label classification experiment;
// the number of classes and their names, feature vectors, label sets, and information about which fold of cross validation (or train/val/test) that each instance belongs to
void exportSelfContainedMLCExperiment(string outFilename, int numFolds, bool useExistingCVFile)
{
	// load the dataset
	vector<string> speciesCodes;
	vector<string> speciesNames;
	pMetaData_Parser::parseSpeciesList("species_list.txt", speciesCodes, speciesNames);
	int numClasses = speciesCodes.size();
	
	vector<pAudioFile_MetaData> metaData;
	pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
	pMetaData_Parser::parseBagLabelSets("bag_labels.txt", metaData);
	pMetaData_Parser::parseBagLevelFeatures("histogram_of_segments.txt", metaData);
	
	if( useExistingCVFile )
		pMetaData_Parser::parseCrossValidationFolds("CVfolds_"+pStringUtils::intToString(numFolds)+".txt", metaData);
	
	cout << "***** dataset summary *****" << endl;
	int numBags = metaData.size();
	int featureDim = metaData[0].bagLevelFeatureVector_.size();
	cout << "# of species = " << numClasses << endl;
	for(int i = 0; i < numClasses; ++i)
		cout << "\t" << i << "," << speciesCodes[i] << "," << speciesNames[i] << endl;
	cout << "# of audio recordings = " << numBags << endl;
	cout << "***************************" << endl;
	
	// assign each recording to a fold for cross-validation
	vector<int> folds;
	for(int i = 0; i < numBags; ++i) folds.push_back(i % numFolds);
	pRandom::shuffle(folds);
	
	pTextFile f(outFilename, PFILE_WRITE);
	
	f << "# you get 2 lines to write a comment about this file. exactly 2." << "\n";
	f << "# second line of comment" << "\n";
	
	f << "num classes = " << numClasses << "\n";
	f << "feature dim = " << featureDim << "\n";
	f << "partition = k-fold, " << numFolds << "\n";
	
	f << "---" << "\n";
	
	for(int j = 0; j < speciesNames.size(); ++j)
		f << j << "," << speciesNames[j] << "\n";
	
	f << "***" << "\n";
	
	f << "index;partition;labelset;featurevector" << "\n";
	
	for(int i = 0; i < numBags; ++i)
	{	
		if( useExistingCVFile )
			f << i << ";" << metaData[i].crossValidationFold_ << ";"; 
		else 
			f << i << ";" << folds[i] << ";";
		
		if( metaData[i].speciesList_.size() == 0 ) f << ";"; // handle case for empty set
		for(int j = 0; j < metaData[i].speciesList_.size(); ++j)
			f << metaData[i].speciesList_[j]  << (j == metaData[i].speciesList_.size() - 1 ? ";" : ",");
		
		for(int j = 0; j < metaData[i].bagLevelFeatureVector_.size(); ++j)
			f << metaData[i].bagLevelFeatureVector_[j]  << (j == metaData[i].bagLevelFeatureVector_.size() - 1 ? "\n" : ",");
	}
	
	f.close();
}


// sometimes it is helpful to pick a random subset of the dataset, for example,
// to be used as training examples for segmentation. This function copies
// numToCopy (without replacement) items from srcFolder to destFolder (whatever type of file is in those... could be wav, or spectrogram)
void copyRandomSubset(string srcFolder, string destFolder, int numToCopy, bool newSeed = false)
{
	if( newSeed ) srand(time(NULL));

	cout << "reading file list..." << endl;
	vector<string> srcFiles;
	pFileUtils::getFilesInDirectory(srcFolder, srcFiles);
	
	cout << "shuffling " << srcFiles.size() << " filenames" << endl;

	pRandom::templatedShuffle<string>(srcFiles);
	
	for(int i = 0; i < numToCopy; ++i)
	{
		cout << "copying " << i << " / " << numToCopy << endl;
		system(string("cp " + srcFolder + "/" + srcFiles[i] + " " + destFolder + "/" + srcFiles[i]).c_str());
	}
}

// this program generates a very large number of files and puts them in several folders. some times, it is necessary to 
// move files around from one folder to another, but there are so many files that it would crash the OS to do it through the
// standard GUI. also, the linux command 'mv' generates an error "-bash: /bin/mv: Argument list too long" when given 60k+ files.
// this function moves all of the files in srcFolder to destFolder, 1 at a time.  
void moveLargeNumberOfFiles(string srcFolder, string destFolder)
{
	vector<string> srcFiles;
	pFileUtils::getFilesInDirectory(srcFolder, srcFiles);
	
	for(int i = 0; i < srcFiles.size(); ++i)
	{
		string command = "mv " + srcFolder + "/" + srcFiles[i] + " " + destFolder;
		cout << i << " / " << srcFiles.size() << ": " << command << endl;
		system(command.c_str());
	}
}

// after selected a random subset of wavs to be used as segmentation examples, it is necessary to have all of the 
// corresponding wav files (and not every other file), in order to share the segmentation examples with someone else.
// this function examines the contents of the /segmentation_examples folder and writes anything else that would be 
// relevant to /segmentation_examples_support. this function assumes that the WAVs of interest are in /chunks (not in /src_wavs).
// if you want to copy data from src_wavs, you will will have to modify this function. 
// this is not a function I expect to call a lot, so it is not exposed through the command line.
void copySegmentationExamplesSupportingFiles()
{
	vector<string> srcFiles;
	pFileUtils::getFilesInDirectory("segmentation_examples", srcFiles);
	
	for(int i = 0; i < srcFiles.size(); ++i)
	{
		string prefix = pStringUtils::firstPartOfSplit(srcFiles[i], ".");
		string command = "cp chunks/" + prefix + ".wav segmentation_examples_support/" + prefix + ".wav"; 
		cout << command << endl;
		system(command.c_str());
	}
}

// some of the bags may be labeled with only 0 or 1 labels. such bags are
// unambiguous examples of a particular class (or of background noise/uknown classes).
// for various reasons, it is useful to copy such examples to a separate folder
// (e.g. to get a visual guide to the species).
void extractUnambiguousExamples()
{
	// load the dataset
	vector<string> speciesCodes;
	vector<string> speciesNames;
	vector<pAudioFile_MetaData> metaData;
	pMetaData_Parser::parseSpeciesList("species_list.txt", speciesCodes, speciesNames);
	pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
	pMetaData_Parser::parseBagLabelSets("bag_labels.txt", metaData);
	int numClasses = speciesCodes.size();
	
	
	system("mkdir ./unambiguous");
	system("mkdir ./unambiguous/empty");
	system("mkdir ./unambiguous/multibird");
	vector<string> subfolders;
	for(int i = 0; i < numClasses; ++i)
	{
		string subfolder = pStringUtils::intToString(i) + "_" + speciesCodes[i] + "_" + pStringUtils::replaceChar(speciesNames[i], " '", '_');
		subfolders.push_back(subfolder);
		string command = "cd ./unambiguous; mkdir " + subfolder; 
		cout << command << endl;
		system(command.c_str());
	}
	
	int numUnambiguous = 0;
	int numEmpty = 0;
	int numMultibird = 0;
	for(int i = 0; i < metaData.size(); ++i)
	{
		if( metaData[i].speciesList_.size() == 0 )
		{
			// empty set
			string command = "cp ./spectrograms/"+metaData[i].filename_+" ./unambiguous/empty;";
			cout << command << endl;
			system(command.c_str());
			++numEmpty;
		}
		else if( metaData[i].speciesList_.size() == 1 )
		{
			int singleSpeciesLabel = metaData[i].speciesList_[0];
			string command = "cp ./spectrograms/"+metaData[i].filename_+" ./unambiguous/" + subfolders[singleSpeciesLabel] + ";";
			cout << command << endl;
			system(command.c_str());
			++numUnambiguous;
		}
		else
		{
			// multi-bird
			string command = "cp ./spectrograms/"+metaData[i].filename_+" ./unambiguous/multibird;";
			cout << command << endl;
			system(command.c_str());
			++numMultibird;
		}
	}
	
	cout << "found " << numUnambiguous << " unambiguous, " << numEmpty << " empty, " <<numMultibird << " multi-bird examples out of " << metaData.size() << endl;
}

// in some cases, the segmentation algorithm finds no segments, but the manual label set specifies that
// some species is present. this is a segmentation false negative. it is useful to see where these cases occur,
// and how often. this function extracts and counts such examples
void extractSegmentationFalseNegatives()
{
	// load the dataset
	vector<string> speciesCodes;
	vector<string> speciesNames;
	vector<pAudioFile_MetaData> metaData;
	pMetaData_Parser::parseSpeciesList("species_list.txt", speciesCodes, speciesNames);
	pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
	pMetaData_Parser::parseBagLabelSets("bag_labels.txt", metaData);
	pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
	
	system("mkdir ./segmentation_false_negative");
	
	int numSegFalsePositives = 0;
	int numNoSegments = 0;
	int numEmptyLabelSet = 0;
	for(int i = 0; i < metaData.size(); ++i)
	{
		if( metaData[i].speciesList_.size() == 0 ) ++numEmptyLabelSet;
		if( metaData[i].segments_.size() == 0 ) ++numNoSegments;
		
		if( metaData[i].speciesList_.size() != 0 && metaData[i].segments_.size() == 0 )
		{
			++numSegFalsePositives;
			string command = "cp ./spectrograms/"+metaData[i].filename_+" ./segmentation_false_negative/";
			cout << command << endl;
			system(command.c_str());
		}
	}
	
	cout << "# of bags with empty label set = " << numEmptyLabelSet << endl;
	cout << "# of bags with 0 segments = " << numNoSegments << endl;
	cout << "# of segmentation false positives = " << numSegFalsePositives << " / " << metaData.size() << endl;
}

//////////////////////////////////////////////////////////////////////////////////

void runMLCExperiment()
{
	vector<string> datasets, classifiers;
	//datasets.push_back("iphone_mlc.txt");
	//datasets.push_back("hja_bird_song_mlc.txt");
	//datasets.push_back("hja_bird_sdm_mlc.txt");
	//datasets.push_back("moths_mlc.txt");
	datasets.push_back("hja20092010histosegments_mlc.txt");
	
	
	classifiers.push_back("binary_relevance");
	//classifiers.push_back("ecc");
	//classifiers.push_back("ecc_oob");
	
	runFullMultiLabelClassifierExperiment(datasets, classifiers);	
}

// a bunch of different stuff that was done not involving command-line, i.e. just using
// the code as an API
void doNonCLIStuff()
{
	//pDatasetSpecific_HJAndrews::analyzeSegmentsPerHourByCluster(250);
	//pDatasetSpecific_HJAndrews::analyzeSegmentsPerHour();
	//pDatasetSpecific_HJAndrews::analyzeSegmentsPerDay(250);
	//pDatasetSpecific_HJAndrews::prepareLabeled2009_2010data_for_species_discovery();
	
	//pMIMLClassDiscoveryExperiment::runExperiment("random", 300, 100);
	//pMIMLClassDiscoveryExperiment::runExperiment("farthest", 300, 1);
	//pMIMLClassDiscoveryExperiment::runExperiment("coverage_oracle", 300, 1);
	//pMIMLClassDiscoveryExperiment::runExperiment("cluster_coverage", 300, 10);
	//pMIMLClassDiscoveryExperiment::runExperiment("generalized_coverage", 300, 1);
	//pMIMLClassDiscoveryExperiment::runExperiment("instance_coverage", 300, 10);
	
	//pMIMLClassDiscoveryExperiment::runExperiment("random", 100, 100, true);
	//pMIMLClassDiscoveryExperiment::runExperiment("random_non_empty", 100, 100, true);
	
	//pMIMLClassDiscoveryExperiment::runExperiment("farthest", 100, 10, true);
	//pMIMLClassDiscoveryExperiment::runExperiment("coverage_oracle", 100, 100, true);
	//pMIMLClassDiscoveryExperiment::runExperiment("cluster_coverage", 100, 100, true);
	//pMIMLClassDiscoveryExperiment::runExperiment("instance_coverage", 100, 10, true);
	//pMIMLClassDiscoveryExperiment::runExperiment("generalized_coverage", 100, 100, true);
	
	//pMIMLClassDiscoveryExperiment::runExperiment("CCMIFF", 100, 100, true);
	
	//pMIMLClassDiscoveryExperiment::runExperiment("farthest2", 100, 100, true);
	//pMIMLClassDiscoveryExperiment::runExperiment("farthest3", 100, 100, true);
	
	//pMIMLClassDiscoveryExperiment::runExperiment("MIFFk1", 100, 100, true);
	//pMIMLClassDiscoveryExperiment::runExperiment("MIFFk2", 100, 100, true);
	
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("cluster_coverage", 100, 100);
	
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("CCMIFF", 100, 100);
	
	//pMIMLClassDiscoveryExperiment::runExperiment("cluster_centers", 100, 100, true);
	
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("random", "random.txt", 100, 1000);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("random_non_empty", "random_non_empty.txt", 100, 1000);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("coverage_oracle", "coverage_oracle.txt", 100, 1000);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("cluster_coverage", "cc_10.txt", 100, 1000, 10);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("cluster_coverage", "cc_20.txt", 100, 1000, 20);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("cluster_coverage", "cc_30.txt", 100, 1000, 30);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("cluster_coverage", "cc_40.txt", 100, 1000, 40);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("cluster_coverage", "cc_50.txt", 100, 1000, 50);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("cluster_coverage", "cc_100.txt", 100, 1000, 100);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("cluster_coverage", "cc_200.txt", 100, 1000, 200);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("cluster_coverage", "cc_300.txt", 100, 1000, 300);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("cluster_coverage", "cc_400.txt", 100, 1000, 400);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("cluster_coverage", "cc_500.txt", 100, 1000, 500);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("MIFF", "MIFF_k1.txt", 100, 1000, 1);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("MIFF", "MIFF_k2.txt", 100, 1000, 2);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("MIFF", "MIFF_k3.txt", 100, 1000, 3);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("CCMIFF", "CCMIFF_clusters=50,k=1.txt", 100, 1000, 50, 1);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("CCMIFF", "CCMIFF_clusters=50,k=2.txt", 100, 1000, 50, 2);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("CCMIFF", "CCMIFF_clusters=100,k=1.txt", 100, 1000, 100, 1);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("CCMIFF", "CCMIFF_clusters=100,k=2.txt", 100, 1000, 100, 2);
	
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("cluster_coverage_tiebreak", "cctb_50.txt", 100, 1000, 50);
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("cluster_coverage_tiebreak", "cctb_500.txt", 100, 1000, 500);
	
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("unified_CCMIFF", "unified_CCMIFF_50_2.txt", 100, 1000, 50, 2);
	
	//pMIMLClassDiscoveryExperiment::runExperiment_multithread("cluster_centers", "cluster_centers_100.txt", 100, 1000);
	
	//pDatasetSpecific_MLSPCompetition::baselineRFBRHisto();
	//pDatasetSpecific_MLSPCompetition::makeBagLabelsTestHidden();
	//pDatasetSpecific_MLSPCompetition::makeSampleSubmissionFile();
	//pDatasetSpecific_MLSPCompetition::setupTrainTest();
	//pDatasetSpecific_MLSPCompetition::makeTrainAndTestSets();
	//pDatasetSpecific_MLSPCompetition::makeSolutionFileForKaggle();
	//pDatasetSpecific_HJAndrews::visualizeRecordingCoverage();
	
	//pDatasetSpecific_MLSPCompetition::visualizeRecordingCoverage();
	
	//pDatasetSpecific_HJAndrews::speciesDiscovery_randomDawn(100);
	//pDatasetSpecific_HJAndrews::speciesDiscovery_unifiedCCMIFF(100, 100, 2, "ccmiff_c=100_k=2");
	//pDatasetSpecific_HJAndrews::speciesDiscovery_unifiedCCMIFF(100, 1000, 2, "ccmiff_c=1000_k=2");
	//pDatasetSpecific_HJAndrews::speciesDiscovery_MIFF(100, 2, "miff_k=2");
	//pDatasetSpecific_HJAndrews::speciesDiscovery_clusterCenters(100, "cluster_centers");
	//pDatasetSpecific_HJAndrews::discoveryAnalysis_reorderDawn();
	
	//pDatasetSpecific_HJAndrews::speciesDiscovery_rainfilter_unifiedCCMIFF(100, 1000, 2, "rainfilter_ccmiff_c=1000_k=2", 0.5);
	//pDatasetSpecific_HJAndrews::speciesDiscovery_rainfilter_unifiedCCMIFF(100, 100, 2, "rainfilter_ccmiff_c=100_k=2", 0.5);
	//pDatasetSpecific_HJAndrews::speciesDiscovery_rainfilter_unifiedCCMIFF(100, 100, 2, "rainfilter_ccmiff_c=100_k=2_rainthreshold=.1", 0.1);
	//pDatasetSpecific_HJAndrews::speciesDiscovery_rainfilter_unifiedCCMIFF(100, 100, 2, "rainfilter_ccmiff_c=100_k=2_rainthreshold=.01", 0.01);
	//pDatasetSpecific_HJAndrews::speciesDiscovery_rainfilter_unifiedCCMIFF(100, 1000, 2, "rainfilter_ccmiff_c=1000_k=2_rainthreshold=.1", 0.1);
	//pDatasetSpecific_HJAndrews::speciesDiscovery_rainfilter_unifiedCCMIFF(100, 1000, 2, "rainfilter_ccmiff_c=1000_k=2_rainthreshold=.01", 0.01);
	
	//pDatasetSpecific_HJAndrews::speciesDiscovery_rainfilter_MIFF(100, 2, "rainfilter_miff_k=2_rainthreshold=.01", 0.01);
	//pDatasetSpecific_HJAndrews::speciesDiscovery_rainfilter_MIFF(100, 2, "rainfilter_miff_k=2_rainthreshold=.1", 0.1);
	
	//pDatasetSpecific_HJAndrews::speciesDiscovery_rainfilter_clusterCenters(100, "cluster_centers_rainthreshold=.1", .1);
	//pDatasetSpecific_HJAndrews::speciesDiscovery_rainfilter_clusterCenters(100, "cluster_centers_rainthreshold=.01", .01);
	
	//pDatasetSpecific_HJAndrews::discoveryAnalysis_countResults("dawn/labels_dawn_processed.csv");
	//pDatasetSpecific_HJAndrews::discoveryAnalysis_countResults("ccmiff_c=1000_k=2/labels_ccmiff_c=1000_k=2_processed.csv");
	//pDatasetSpecific_HJAndrews::discoveryAnalysis_countResults("ccmiff_c=100_k=2/labels_ccmiff_c=100_k=2_processed.csv"); // TODO
	//pDatasetSpecific_HJAndrews::discoveryAnalysis_countResults("cluster_centers/labels_clustercenters_processed.csv");
	//pDatasetSpecific_HJAndrews::discoveryAnalysis_countResults("miff_k=2/labels_processed.csv");
	
	//pDatasetSpecific_HJAndrews::discoveryAnalysis_countResults("rain-filtered/cluster_centers_rainthreshold=.01/labels_processed.csv");
	//pDatasetSpecific_HJAndrews::discoveryAnalysis_countResults("rain-filtered/cluster_centers_rainthreshold=.1/labels_processed.csv");
	
	//pDatasetSpecific_HJAndrews::discoveryAnalysis_countResults("rain-filtered/rainfilter_ccmiff_c=1000_k=2_rainthreshold=.01/labels_processed.csv");
	//pDatasetSpecific_HJAndrews::discoveryAnalysis_countResults("rain-filtered/rainfilter_ccmiff_c=1000_k=2_rainthreshold=.1/labels_processed.csv");
	
	//pDatasetSpecific_HJAndrews::discoveryAnalysis_countResults("rain-filtered/rainfilter_miff_k=2_rainthreshold=.01/labels_processed.csv");
	//pDatasetSpecific_HJAndrews::discoveryAnalysis_countResults("rain-filtered/rainfilter_miff_k=2_rainthreshold=.1/labels_processed.csv");
	
	//pDatasetSpecific_HJAndrews::discoveryAnalysis_joinSelectionHistograms();
	//pDatasetSpecific_HJAndrews::discoveryAnalysis_makeSelectionHistoGraph();
}

int main (int argc, char * const argv[]) 
{	
	srand(0xDEADBEEF); // out of superstition, I prefer to seed with DEADBEEF instead of 0.
	
	string useageText = 
		"\n\n\n------------------------------------------------------------------" "\n"
		"Pterodactyl 0.1" "\n"
		"Tools for bioacoustics and machine learning" "\n"
		"By Forrest Briggs" "\n"
		"MIT Software License" "\n"
		"\nSee notes below for usage:" "\n\n"
		"==================================================================\n\n\n" 
		
		"Create the folders required for data" "\n"
		"Pterodactyl -setup" "\n" 
		
		"\n------------------------------------------------------------------\n\n"
		
		"Split wav files into smaller chunks:" "\n"
		"Pterodactyl -split [duration_sec] [source_folder] [dest_folder]" "\n"
		"Example:" "\n"
		"Pterodactyl -split 10 src_wavs chunks" "\n"
		
		"\n------------------------------------------------------------------\n\n"
		
		"Generate spectrograms from WAV audio files:" "\n"
		"Pterodactyl -spectrogram [source_folder] [dest_folder]" "\n"
		"Example:" "\n"
		"Pterodactyl -spectrogram chunks spectrograms" "\n"
		
		"\n------------------------------------------------------------------\n\n"
		
		"Apply a noise-reduction filter:" "\n"
		"Pterodactyl -noisefilter [source_folder] [dest_folder]" "\n"
		"Example:"
		"Pterodactyl -noisefilter spectrograms filtered_spectrograms" "\n"
		
		"\n------------------------------------------------------------------\n\n" 
		
		
		"Copy a random subset of files from one folder to another:" "\n"
		"Pterodactyl -copy_random_subset [source_folder] [dest_folder] [numToCopy]" "\n"
		"Notes:"
		"The main use of this function is to copy some random spectrograms to /segmentation_examples for the purpose of labeling. There might be some other uses too."
		"\n------------------------------------------------------------------\n\n" 
		
		
		"Train 2D-Random Forest segmentation:" "\n"
		"Pterodactyl -2drfseg_train [spectrograms_folder] [examples_folder]" "\n"
		"Example:" "\n"
		"Pterodactyl -2drfseg_train filtered_spectrograms segmentation_examples" "\n"
		"Notes:" "\n"
		"Requires some folders to exist / have files in them already:" "\n"
		"	./[spectrograms_folder]			-- should contain all spectrograms that we want to process, after applying the noise filter." "\n"
		"	./[examples_folder]			-- should contain a subset of the spectrograms, marked with red (positive examples) or blue (negative examples)." "\n"
		"Results will be saved to segmentation_RF_model.txt"  "\n"
		
		"\n------------------------------------------------------------------\n\n"
		
		"Run 2D-Random Forest segmentation:" "\n"
		"Pterodactyl -2drfseg_run [spectrograms_folder]" "\n"
		"Example:" "\n"
		"Pterodactyl -2drfseg_run filtered_spectrograms" "\n"
		"Notes:" "\n"
		"Requires some folders to exist / have files in them already:" "\n"
		"	./[spectrograms_folder			-- should contain all spectrograms that we want to process, after applying the noise filter." "\n"
		"	./segmentation_RF_model.txt		-- saved random forest model from training should exist already" "\n"
		"Results will be written to several folders (which must already exist):" "\n"
		"	./segmentation_probabilities		-- greyscale probability maps for segmentation (before thresholding)" "\n"
		"	./segmentation_mask			-- binary masks of segments" "\n"
		"	./segmentation_outlines			-- spectrograms with outlined segments" "\n"
		
		"\n------------------------------------------------------------------\n\n"
		
		"Run segment analysis tools (extract segments from the spectrogram in the form of smaller images, get segment features, etc):" "\n"
		"Pterodactyl -analyze_segments [spectrograms_folder] [masks_folder] [[bagidfilename]]" "\n"
		"Example:" "\n"
		"Pterodactyl -analyze_segments filtered_spectrograms segmentation_mask" "\n"
		"Notes:" "\n"
		"The argument [[bagidfilename]] is optional. If it is supplied, the filename given will be loaded as a proxy for bag_id2filename.txt." "\n"
		"If the argument is not supplied, bag_id2filename.txt will be generated by this function, and assigns bag_id's to filenames in whatever order they are read by the OS directory reading functions"
		"Writes results to the following folderes:" "\n"
		"	./segments				-- images of each individual segment with everything else removed" "\n"
		"	./segment_masks			-- binary masks of each individual segment" "\n"
		
		"\n------------------------------------------------------------------\n\n"
		
		"Visualization: generate a mosaic of segments" "\n"
		"Pterodactyl -segment_mosaic" "\n"
		"Notes:" "\n"
		"Requires segments images to already be in the /segments folder. Run -analyze_segments before running this to setup the required data."
		"Writes results to segment_mosaic.bmp" "\n"

		"\n------------------------------------------------------------------\n\n"
		
		"Visualization: apply k-means++ clustering to segment features, then produce a mosaic image of each cluster" "\n"
		"Pterodactyl -cluster_segments [k]" "\n"
		"Notes:" "\n"
		"Requires segments images to already be in the /segments folder. Run -analyze_segments before running this to setup the required data."
		"Writes results to /clustered_segments" "\n"
	
		"\n------------------------------------------------------------------\n\n"
		
		"Generate a feature vector for each audio recording by aggregating the features of its segments, using a histogram of segments (k-means codebook)." "\n"
		"Pterodactyl -histogram_of_segments" "\n"
		"Notes:" "\n"
		"Requires bag_id2filename.txt and segment_features.txt"

		"\n------------------------------------------------------------------\n\n"
		
		"Export a single self-contained file with everything needed to run a multi-label classification experiment." "\n"
		"Pterodactyl -export_mlc_experiment [outfilename] [num-folds]" "\n"
		"Notes: [outfilename] -- the file that is saved. [num-folds] -- k for k-fold cross-validation." "\n"
		"Requires everything that saves a text file above this in the list already to be run."

		"\n------------------------------------------------------------------\n\n"
		
		"Run a multi-label classification experiment, comparing different classifiers on a ground-truth label set using k-fold cross-validation." "\n"
		"Pterodactyl -mlc_experiment" "\n"
		"Notes:" "\n"
		"Requires species_list.txt, bag_id2filename.txt, histogram_of_segments.txt, bag_labels.txt, CV_fold_10.txt (or similar CV_fold_[k].txt) "		

		"\n==================================================================\n\n\n";
		
	PTIMER_MARK timer = pTimer::instance()->setMark();
		
	p2DSegmentFeatures::initGradientCodebook(); // this must be called to setup HOG features
		
	if( argc == 1 )																				{ cout << useageText << endl; return 0; }
	else if( string(argv[1])	== "-setup")													setupFolders();
	else if( string(argv[1])	== "-extract_unambiguous")										extractUnambiguousExamples();
	else if( string(argv[1])	== "-extract_seg_false_negative")								extractSegmentationFalseNegatives();

	else if( string(argv[1])	== "-split"				&& argc == 5 )							split(pStringUtils::stringToInt(argv[2]), argv[3], argv[4]);
	else if( string(argv[1])	== "-split_multithread"	&& argc == 5 )							split_multithread(pStringUtils::stringToInt(argv[2]), argv[3], argv[4]);
	
	else if( string(argv[1])	== "-analyze_segments"	&& argc == 4)							processSegments(argv[2], argv[3], "");
	else if( string(argv[1])	== "-analyze_segments"	&& argc == 5)							processSegments(argv[2], argv[3], argv[4]);
	else if( string(argv[1])	== "-spectrogram"		&& argc == 4)							makeSpectrograms(argv[2], argv[3]);

	else if( string(argv[1])	== "-noisefilter"		&& argc == 4)							noiseReductionFilter(argv[2], argv[3]);

	else if( string(argv[1])	== "-spectrogram+noisefilter_multithread" && argc == 5)			makeSpectrogramsAndApplyNoiseFilter_multithread(argv[2], argv[3], argv[4]);

	else if( string(argv[1])	== "-copy_random_subset" && argc == 5)							copyRandomSubset(argv[2], argv[3], pStringUtils::stringToInt(argv[4]));
	
	else if( string(argv[1])	== "-copy_random_subset_reseed" && argc == 5)					copyRandomSubset(argv[2], argv[3], pStringUtils::stringToInt(argv[4]), true);
	
	else if( string(argv[1])	== "-move_files"		&& argc == 4)							moveLargeNumberOfFiles(argv[2], argv[3]);

	else if( string(argv[1])	== "-2drfseg_train"		&& argc == 4)							trainRFSegmentation(argv[3], argv[2], "segmentation_RF_model.txt");
	else if( string(argv[1])	== "-2drfseg_run"		&& argc == 3)							runRFSegmentation("segmentation_RF_model.txt", argv[2]);
	else if( string(argv[1])	== "-2drfseg_run_multithread" && argc == 3)						runRFSegmentation_multithread("segmentation_RF_model.txt", argv[2]);
	else if( string(argv[1])	== "-segment_mosaic")											generateSegmentMosaic();
	else if( string(argv[1])	== "-cluster_segments"	&& argc == 3)							clusterSegments(pStringUtils::stringToInt(argv[2]));
	
	else if( string(argv[1])	== "-species_discovery_cluster" && argc == 3)					speciesDiscovery_segment_cluster(pStringUtils::stringToInt(argv[2]));
	else if( string(argv[1])	== "-species_discovery_farthest" && argc == 3)					speciesDiscovery_segment_farthest(pStringUtils::stringToInt(argv[2]));
	else if( string(argv[1])	== "-species_discovery_random" && argc == 3)					speciesDiscovery_segment_random(pStringUtils::stringToInt(argv[2]));
	
	else if( string(argv[1])	== "-species_discovery_chunk_cluster" && argc == 3)				speciesDiscovery_chunk_cluster(pStringUtils::stringToInt(argv[2]));
	else if( string(argv[1])	== "-species_discovery_chunk_farthest" && argc == 3)			speciesDiscovery_chunk_farthest(pStringUtils::stringToInt(argv[2]));
	else if( string(argv[1])	== "-species_discovery_chunk_random" && argc == 3)				speciesDiscovery_chunk_random(pStringUtils::stringToInt(argv[2]));
	
	else if( string(argv[1])	== "-histogram_of_segments")									histogramOfSegmentFeatures();
	else if( string(argv[1])	== "-more_bag_features")										makeAdditionalBagLevelFeatures();
	
	else if( string(argv[1])	== "-rain_classifier")											runRainClassifier();

	else if( string(argv[1])	== "-export_mlc_experiment" && argc == 4)						exportSelfContainedMLCExperiment(argv[2], pStringUtils::stringToInt(argv[3]), true);
	else if( string(argv[1])	== "-mlc_experiment")											runMLCExperiment(); // hidden option available only to those who read the source code :P. with great power comes great responsibility.
	else cout << "error: invalid command line arguments, argc = " << argc << endl;

	cout << "*****" << endl << "Total runtime: " << pTimer::instance()->timeSinceMark(timer) / 1000.0f << " seconds" << endl; 

    return 0;
}




