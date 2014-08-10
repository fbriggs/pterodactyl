// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// parameters which affect the tools in this software are generally stored here, e.g. spectrogram parameters, feature vector parameters,
// parameters for machine learning algorithms, etc.
// (a few may be left behind in other places; further revisions of this code should move any additional parameters here)

#pragma once 

class pParameters
{
public: 
	const static int windowSize = 512;												// spectrogram parameters
	const static int windowStep = 128;
		
	const static int bandPassLow = 0;
	const static int bandPassHigh = 0;
	
	//const static int segmentationNumRandomForestTrees = 100;						// this is the parameter from the ICML 2013 workshop paper
	const static int segmentationNumRandomForestTrees = 16;							// settings for full 2009
	
	
	const static int probmapBlurRadius = 3;
	const static int segmentMosaicImageWidth = 7.5 * 150;
	//const static int clusterSegmentsVisualizationK = 100;							// k for kmeans in clustered segment visualization
	const static int histogramOfSegmentsClusterK = 100;
	const static int kMeansMaxIterations = 100;
	
	const static int defaultCrossValidationFolds = 10;
	
	
	const static int binaryRelevanceNumTrees = 25 * 25;								// number of trees used for each class in pBinaryRelevance classifier
	const static int binaryRelevance_numRepetitions = 10;							// number of trees used for each class in pBinaryRelevance classifier
	
	const static int pMultiLabelEnsembleOfClassifierChains_numChains = 25;			// parameters for ECC algorithm for multi-label classification
	const static int pMultiLabelEnsembleOfClassifierChains_treesPerChain = 25;	
	const static int pMultiLabelEnsembleOfClassifierChains_numRepetitions = 10;

	const static float segmentationProbabilityThreshold;							// you can't set the value of const static float the same way as an int in C++. the value is in pParameters.cpp
};

