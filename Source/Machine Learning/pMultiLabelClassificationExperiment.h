// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// this file contains functions related to running a multi-label classification experiment.

#pragma once

#include "pTextFile.h"
#include "pSignalUtils.h"
#include "pExample_MultiLabel.h"
#include "pBinaryRelevance.h"
#include "pMultiLabelEnsembleOfClassifierChains.h"
#include "pTrainValTest.h"

namespace multi_label_classification_experiment
{

// reads all of the data in a the "self-contained multi-label classification experiment" file format.
void parseSelfContainedMultiLabelClassificationExperiment(string filename, int& _numClasses, int& _featureDim, vector<string>& _classNames, vector<pExample_MultiLabel>& _examples, string& _cvMode, int& _cvFolds)
{
	pTextFile f(filename, PFILE_READ);
	
	/// parse the file header ///
	
	string comment = f.readLine() + "\n" + f.readLine();

	vector<string> numClassesParts = pStringUtils::splitNonEmpty(f.readLine(), "=");
	assert(numClassesParts.size() == 2);
	assert(pStringUtils::pack(numClassesParts[0]) == "num classes");
	_numClasses = pStringUtils::stringToInt(numClassesParts[1]);
	
	vector<string> featureDimParts = pStringUtils::splitNonEmpty(f.readLine(), "=");
	assert(featureDimParts.size() == 2);
	assert(pStringUtils::pack(featureDimParts[0]) == "feature dim");
	_featureDim = pStringUtils::stringToInt(featureDimParts[1]);

	vector<string> partitionParts = pStringUtils::splitNonEmpty(f.readLine(), "=, "); 

	if( partitionParts.size() == 3 && partitionParts[0] == "partition" && partitionParts[1] == "k-fold" )
	{
		_cvMode = "k-fold";
		_cvFolds = pStringUtils::stringToInt(partitionParts[2]);
	}
	else if( partitionParts.size() == 2 && partitionParts[0] == "partition" && partitionParts[1] == "train-val-test" )
	{
		_cvMode = "train-val-test";
		_cvFolds = 0;
	}
	else 
		assert(false && "partition format unrecognized");
	
	assert(f.readLine() == "---" && "expected --- after header");
	
	cout << comment << endl;
	cout << "num classes = " << _numClasses << endl;
	cout << "feature dim = " << _featureDim << endl;
	cout << "cv mode = " << _cvMode << endl;
	cout << "cv folds = " << _cvFolds << endl;
	
	/// parse the class-name list ///
	
	for(int i = 0; i < _numClasses; ++i)
	{
		vector<string> parts = pStringUtils::split(f.readLine(), ",");
		assert(parts.size() == 2);
		
		int index = pStringUtils::stringToInt(parts[0]);
		assert(index == i);
		
		_classNames.push_back(parts[1]); 
	}
	
	assert(f.readLine() == "***" && "expected *** after class list");
	
	for(int i = 0; i < _classNames.size(); ++i)
		cout << _classNames[i] << endl;
	
	assert(f.readLine() == "index;partition;labelset;featurevector");
	
	/// parse the main bulk of the dataset ///
	
	int index = 0;
	while(!f.eof())
	{
		// verify the index
		string line = f.readLine();
		vector<string> parts = pStringUtils::split(line, ";");
		if( parts.size() != 4 )
		{
			cout << "error: " << line << endl;
		}
		assert(parts.size() == 4 && "example does not match the pattern index;partition;labelset;featurevector");
		
		assert(pStringUtils::stringToInt(parts[0]) == index && "examples must be listed sequentially in increasing order of index");

		// parse the partition for cross-validation or train/test/val
		int partition;
		if( _cvMode == "k-fold" )
		{
			partition = pStringUtils::stringToInt(parts[1]);
		}
		
		if( _cvMode == "train-val-test" )
		{
			if( parts[1] == "train" )		partition = kTrain;
			else if( parts[1] == "test" )	partition = kTest;
			else if( parts[1] == "val" )	partition = kValidation;
			else assert(false && "unrecognized partion; must be 'train' 'test' or 'val'");
		}	
		
		// parse the label set
		vector<string> labelSetParts = pStringUtils::split(parts[2], ",");
		vector<int> labelSet;
		for(int j = 0; j < labelSetParts.size(); ++j)
			labelSet.push_back(pStringUtils::stringToInt(labelSetParts[j]));
		
		// parse the feature vector
		vector<string> featureVectorParts = pStringUtils::split(parts[3], ",");
		vector<float> fv;
		for(int j = 0; j < featureVectorParts.size(); ++j)
			fv.push_back(pStringUtils::stringToFloat(featureVectorParts[j]));
		
		_examples.push_back(pExample_MultiLabel(fv, labelSet, partition));
		
		++index;
	}
	
	cout << "parsed " << _examples.size() << " training examples" << endl;
	

	f.close();
	
	/// output some statistics about the dataset ///
	
	cout << "***" << endl;
	
	// count the number of examples with each class
	vector<int> exampleCounts(_numClasses);
	for(int i = 0; i < _examples.size(); ++i)
	{
		for(int j = 0; j < _examples[i].labels_.size(); ++j)
			++exampleCounts[_examples[i].labels_[j]];
	}
	
	for(int j = 0; j < _numClasses; ++j)
		cout << j << "," << _classNames[j] << "," << exampleCounts[j] << endl;
}

// helper functions for multiLabelClassifierExperiment which allow it to be templated with classifiers that 
// take different parameters to train. define a template sepcialization for each classifier that might be used in multiLabelClassifierExperiment
template <class TClassifier, class TClassifierParams> void trainMultiLabelClassifier(TClassifier& classifier, vector<pExample_MultiLabel>& trainExs, int numClasses, TClassifierParams& params) 
{ 
	assert(false && "template specialization required"); 
}

template <> void trainMultiLabelClassifier<pBinaryRelevance, pBinaryRelevance_Params>
	(pBinaryRelevance& classifier, vector<pExample_MultiLabel>& trainExs, int numClasses, pBinaryRelevance_Params& params)
{ 
	classifier.train(trainExs, numClasses, params.numTrees_); 
}

template <> void trainMultiLabelClassifier<pMultiLabelEnsembleOfClassifierChains, pMultiLabelEnsembleOfClassifierChains_Params>
	(pMultiLabelEnsembleOfClassifierChains& classifier, vector<pExample_MultiLabel>& trainExs, int numClasses, pMultiLabelEnsembleOfClassifierChains_Params& params)
{ 
	classifier.train(trainExs, numClasses, params.numTreesPerChain_, params.numChains_, params.useOutOfBagCorrection_); 
}

// runs an experiment to evaluate the accuracy of a collection of multi-label classifiers
// on the dataset. a recording/bag-level label file is required, in addition to some kind of bag-level features
// e.g., histogram of segments. 
// the template parameters are:
// TClassifier			- the type of a multi-label classifier
// TClassifierParams	- the type of a class which stores parameters for use in the classifer
template <class TClassifier, class TClassifierParams>
void multiLabelClassifierExperiment_kFold(
	string datasetName,
	string classifierName,
	int repetitionNum,
	vector<string>& _resultsSummary, // summaries of the experiment will be written to this vector of strings (so they can all be printed out together at the end)
	int numClasses,
	int featureDim,
	int cvFolds,
	vector<pExample_MultiLabel>& allExamples,
	TClassifierParams params) 
{
	float hammingLoss_avg = 0;
	float rankLoss_avg = 0;
	float set01loss_avg = 0;
	float oneError_avg = 0;
	float coverage_avg = 0;
	int TP = 0, TN = 0, FP = 0, FN = 0; // true/false positives/negatives
	
	for(int fold = 0; fold < cvFolds; ++fold)
	{
		cout << "*** fold " << fold << " / " << cvFolds << " ***" << endl;
		
		vector<pExample_MultiLabel> trainExs;
		vector<pExample_MultiLabel> testExs;
		for(int i = 0; i < allExamples.size(); ++i)
		{
			if( allExamples[i].cvFold_ == fold )	testExs.push_back(allExamples[i]);
			else									trainExs.push_back(allExamples[i]);
		}
		
		TClassifier classifier;
		trainMultiLabelClassifier<TClassifier>(classifier, trainExs, numClasses, params); // the templating here handles different classifiers with different parameters
		
		// apply the classifier to the test set and update result statistics
		for(int i = 0; i < testExs.size(); ++i)
		{
			vector<int> predictedLabelSet = classifier.classify(testExs[i].fv_);
			vector<float> scores = classifier.getClassScores(testExs[i].fv_);
			vector<int>& groundTruthLabelSet = testExs[i].labels_;
					
			hammingLoss_avg += pMultiLabelErrorMeasure::hammingLoss(predictedLabelSet, groundTruthLabelSet, numClasses);
			rankLoss_avg	+= pMultiLabelErrorMeasure::rankLoss(scores, groundTruthLabelSet, numClasses);
			set01loss_avg	+= pMultiLabelErrorMeasure::set01loss(predictedLabelSet, groundTruthLabelSet, numClasses);
			oneError_avg	+= pMultiLabelErrorMeasure::oneError(scores, groundTruthLabelSet);
			coverage_avg	+= pMultiLabelErrorMeasure::coverage(scores, groundTruthLabelSet, numClasses);
			pMultiLabelErrorMeasure::true_false_pos_neg(predictedLabelSet, groundTruthLabelSet, numClasses, TP, TN, FP, FN);
			
			pVectorUtils::cout_vector(groundTruthLabelSet);
			cout << "->";
			pVectorUtils::cout_vector(predictedLabelSet);
			cout << endl;
		}
	}
	
	// finish calculation of result statistics
	hammingLoss_avg		/= allExamples.size();
	rankLoss_avg		/= allExamples.size();
	set01loss_avg		/= allExamples.size();
	oneError_avg		/= allExamples.size();
	coverage_avg		/= allExamples.size();
	float accuracy = 1.0f - hammingLoss_avg;
	
	cout << "**************************" << endl;
	cout << "accuracy = " << accuracy << endl;
	cout << "avg set 01 loss = " << set01loss_avg << endl;
	cout << "TP = " << TP << endl;
	cout << "TN = " << TN << endl;
	cout << "FP = " << FP << endl;
	cout << "FN = " << FN << endl;	
	
	_resultsSummary.push_back(datasetName + "\t" + classifierName + "\t" + pStringUtils::intToString(repetitionNum) + "\t" + pStringUtils::floatToString(hammingLoss_avg) + "\t" + pStringUtils::floatToString(rankLoss_avg) + "\t" + pStringUtils::floatToString(set01loss_avg) + 
		"\t" + pStringUtils::floatToString(oneError_avg) + "\t" + pStringUtils::floatToString(coverage_avg));
//		pStringUtils::intToString(TP) + "\t" + pStringUtils::intToString(TN) + "\t" + pStringUtils::intToString(FP) + "\t" + pStringUtils::intToString(FN));
}

template <class TClassifier, class TClassifierParams>
void multiLabelClassifierExperiment_trainValTest(
	string datasetName,
	string classifierName,
	int repetitionNum,
	vector<string>& _resultsSummary, // summaries of the experiment will be written to this vector of strings (so they can all be printed out together at the end)
	int numClasses,
	int featureDim,
	vector<pExample_MultiLabel>& allExamples,
	TClassifierParams params) 
{
	float hammingLoss_avg = 0;
	float rankLoss_avg = 0;
	float set01loss_avg = 0;
	float oneError_avg = 0;
	float coverage_avg = 0;
	int TP = 0, TN = 0, FP = 0, FN = 0; // true/false positives/negatives
	
	// split the dataset up into train/test/val
	vector<pExample_MultiLabel> trainExs;
	vector<pExample_MultiLabel> testExs;
	vector<pExample_MultiLabel> valExs; // TODO: for now, val is not used at all in the experiment. instead the experiment is just run with no parameter tuning / fixed parameters.
	
	for(int i = 0; i < allExamples.size(); ++i)
	{
		switch(allExamples[i].cvFold_)
		{
		case kTrain:		trainExs.push_back(allExamples[i]); break;
		case kTest:			testExs.push_back(allExamples[i]);  break;
		case kValidation:	valExs.push_back(allExamples[i]);	break;
		default: assert(false);
		}
	}
	
	// train
	TClassifier classifier;
	trainMultiLabelClassifier<TClassifier>(classifier, trainExs, numClasses, params); // the templating here handles different classifiers with different parameters
	
	// apply the classifier to all of the "test" examples
	for(int i = 0; i < testExs.size(); ++i)
	{
		vector<int> predictedLabelSet = classifier.classify(testExs[i].fv_);
		vector<float> scores = classifier.getClassScores(testExs[i].fv_);
		vector<int>& groundTruthLabelSet = testExs[i].labels_;

		hammingLoss_avg += pMultiLabelErrorMeasure::hammingLoss(predictedLabelSet, groundTruthLabelSet, numClasses);
		rankLoss_avg	+= pMultiLabelErrorMeasure::rankLoss(scores, groundTruthLabelSet, numClasses);
		set01loss_avg	+= pMultiLabelErrorMeasure::set01loss(predictedLabelSet, groundTruthLabelSet, numClasses);
		oneError_avg	+= pMultiLabelErrorMeasure::oneError(scores, groundTruthLabelSet);		
		coverage_avg	+= pMultiLabelErrorMeasure::coverage(scores, groundTruthLabelSet, numClasses);
		pMultiLabelErrorMeasure::true_false_pos_neg(predictedLabelSet, groundTruthLabelSet, numClasses, TP, TN, FP, FN);
		
		pVectorUtils::cout_vector(groundTruthLabelSet);
		cout << "->";
		pVectorUtils::cout_vector(predictedLabelSet);
		cout << endl;
	}
	
	// finish calculation of result statistics
	hammingLoss_avg		/= testExs.size();
	rankLoss_avg		/= testExs.size();
	set01loss_avg		/= testExs.size();
	oneError_avg		/= testExs.size();
	coverage_avg		/= allExamples.size();
	float accuracy = 1.0f - hammingLoss_avg;
	
	cout << "**************************" << endl;
	cout << "accuracy = " << accuracy << endl;
	cout << "avg set 01 loss = " << set01loss_avg << endl;
	cout << "TP = " << TP << endl;
	cout << "TN = " << TN << endl;
	cout << "FP = " << FP << endl;
	cout << "FN = " << FN << endl;	
	
	_resultsSummary.push_back(datasetName + "\t" + classifierName + "\t" + pStringUtils::intToString(repetitionNum) + "\t" + pStringUtils::floatToString(hammingLoss_avg) + "\t" + pStringUtils::floatToString(rankLoss_avg) + "\t" + pStringUtils::floatToString(set01loss_avg) + 
		"\t" + pStringUtils::floatToString(oneError_avg) + "\t" + pStringUtils::floatToString(coverage_avg));
	//	pStringUtils::intToString(TP) + "\t" + pStringUtils::intToString(TN) + "\t" + pStringUtils::intToString(FP) + "\t" + pStringUtils::intToString(FN));
}

// dispatch to one of the two functions above based on cvMode = k-fold or train/test/val
template <class TClassifier, class TClassifierParams>
void multiLabelClassifierExperiment_cvModeDispatch(
	string cvMode,
	string datasetName,
	string classifierName,
	int repetitionNum,
	vector<string>& _resultsSummary, // summaries of the experiment will be written to this vector of strings (so they can all be printed out together at the end)
	int numClasses,
	int featureDim,
	int cvFolds,
	vector<pExample_MultiLabel>& allExamples,
	TClassifierParams params)
{
	if( cvMode == "k-fold" )
		multiLabelClassifierExperiment_kFold<TClassifier, TClassifierParams>(datasetName, classifierName, repetitionNum, _resultsSummary, numClasses, featureDim, cvFolds, allExamples, params);
			
	else if( cvMode == "train-val-test" )
		multiLabelClassifierExperiment_trainValTest<TClassifier, TClassifierParams>(datasetName, classifierName, repetitionNum, _resultsSummary, numClasses, featureDim, allExamples, params);
	
	else assert(false && "unrecognized partitioning scheme");
}

// do everything involved in running a full multi-label classification experiment
void runFullMultiLabelClassifierExperiment(vector<string>& datasets, vector<string>& classifiers)
{
	vector<string> resultsSummary;
	
	resultsSummary.push_back("dataset\tclassifier\trepetition\tHamming_Loss\tRank_Loss\tSet_01_Loss\tOne_Error\tCoverage");
	
	for(int ds = 0; ds < datasets.size(); ++ds)
	{
		cout << "***** running experiments for dataset " << datasets[ds] << " ****" << endl;
		
		int numClasses;
		int featureDim;
		string cvMode;
		int cvFolds;
		vector<string> classNames;
		vector<pExample_MultiLabel> examples;
		parseSelfContainedMultiLabelClassificationExperiment("mlc_experiments/"+datasets[ds], numClasses, featureDim, classNames, examples, cvMode, cvFolds);
		
		
		cout << "cvMode = " << cvMode << endl;
		cout << "cvFolds = " << cvFolds << endl;
		
		for(int i = 0; i < classifiers.size(); ++i)
		{
			if( classifiers[i] == "binary_relevance" )
			{
				for(int rep = 0; rep < pParameters::binaryRelevance_numRepetitions; ++rep)
					multiLabelClassifierExperiment_cvModeDispatch<pBinaryRelevance, pBinaryRelevance_Params>(
						cvMode, datasets[ds], classifiers[i], rep, resultsSummary, numClasses, featureDim, cvFolds, examples, 
							pBinaryRelevance_Params(						pParameters::binaryRelevanceNumTrees));
			}
			
			if( classifiers[i] == "ecc" )
			{
				for(int rep = 0; rep < pParameters::pMultiLabelEnsembleOfClassifierChains_numRepetitions; ++rep)
					multiLabelClassifierExperiment_cvModeDispatch<pMultiLabelEnsembleOfClassifierChains, pMultiLabelEnsembleOfClassifierChains_Params>(
						cvMode, datasets[ds], classifiers[i], rep, resultsSummary, numClasses, featureDim, cvFolds, examples, 
							pMultiLabelEnsembleOfClassifierChains_Params(	pParameters::pMultiLabelEnsembleOfClassifierChains_treesPerChain, 
																			pParameters::pMultiLabelEnsembleOfClassifierChains_numChains,
																			false));
			}
																		
																		
			if( classifiers[i] == "ecc_oob" )
			{
				for(int rep = 0; rep < pParameters::pMultiLabelEnsembleOfClassifierChains_numRepetitions; ++rep)
					multiLabelClassifierExperiment_cvModeDispatch<pMultiLabelEnsembleOfClassifierChains, pMultiLabelEnsembleOfClassifierChains_Params>(
						cvMode, datasets[ds], classifiers[i], rep, resultsSummary, numClasses, featureDim, cvFolds, examples, 
							pMultiLabelEnsembleOfClassifierChains_Params(	pParameters::pMultiLabelEnsembleOfClassifierChains_treesPerChain, 
																			pParameters::pMultiLabelEnsembleOfClassifierChains_numChains,
																			true));
			}
			
			resultsSummary.push_back("\n"); // add some blank lines
			resultsSummary.push_back("\n");
			
			cout << "********** partial results summary:" << endl;	
			for(int i = 0; i < resultsSummary.size(); ++i)
				cout << resultsSummary[i] << endl;
			cout << "********" << endl;
		}
	}
	
	cout << "********** all experiments finished. results summary:" << endl;
	for(int i = 0; i < resultsSummary.size(); ++i)
		cout << resultsSummary[i] << endl;
}


};