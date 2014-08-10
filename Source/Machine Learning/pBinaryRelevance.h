// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// binary relevance is the simplest algorithm for multi-label classification- treat each label independently
// this is implemented with a random forest as the base-classifier.

#pragma once

#include "pVectorUtils.h"
#include "pExample.h"
#include "pExample_MultiLabel.h"
#include "pRandomForest_pthread.h"
#include "pRandomForest_OutOfBag.h"

class pBinaryRelevance_Params // this struct stores the parameteres of a pBinary relevance classifier, in case it is necessary to pass them around
{
public:
	int numTrees_;
	
	pBinaryRelevance_Params(int numTrees) : numTrees_(numTrees) {}
};

class pBinaryRelevance
{
public:
	int numClasses_;
	vector<pRandomForest_pthread*> rfs_;
	
	vector<float> thresholds_; // one threshold for each class
	
	void train(vector<pExample_MultiLabel>& exs, int numClasses, int numTrees)
	{
		numClasses_ = numClasses;
		
		for(int j = 0; j < numClasses; ++j)
		{
			cout << "building binary dataset for class " << j << endl;
			
			vector<pExample> trainingExamplesForClassJ;
			for(int i = 0; i < exs.size(); ++i)
			{
				int Yj = pVectorUtils::vector_contains(exs[i].labels_, j) ? 1 : 0;
				trainingExamplesForClassJ.push_back(pExample(exs[i].fv_, Yj));
			}
			
			cout << "training binary classifier for class " << j << endl;
			
			rfs_.push_back(new pRandomForest_pthread());
			bool storeHistogramInLeaf = true;
			int maxDepth = 15;
			bool useClassBalancedBootstrapping = false;
			rfs_[j]->train_multithread(numTrees, 2, trainingExamplesForClassJ, storeHistogramInLeaf, maxDepth, useClassBalancedBootstrapping); // use default RF params
		}
		
		// find the optimal threshold for each class
		for(int j = 0; j < numClasses; ++j)
		{
			// simple threshold
			//thresholds_.push_back(0.5);
			
			cout << "computing threshold for class " << j << endl;
			
			vector<float> candidateThresholds;
			for(float thresh = 0.001; thresh < 1.0f; thresh += 0.001) candidateThresholds.push_back(thresh);
			
			vector<float> oobEstimates;
			for(int i = 0; i < exs.size(); ++i)
				oobEstimates.push_back(pRandomForest_OutOfBag::estimateClassProbabilitiesOOB(*rfs_[j], exs[i].fv_, i)[1]);
			
			// evaluate each candidate threshold to see which one gives the best binary accuracy w.r.t. class j
			int minError = INT_MAX;
			float bestThreshold = 0.5;
			
			for(int k = 0; k < candidateThresholds.size(); ++k)
			{
				// compute the number of examples on the wrong side of the threshold
				int err = 0;
				for(int i = 0; i < exs.size(); ++i)
				{
					int binaryLabel = pVectorUtils::vector_contains(exs[i].labels_, j) ? 1 : 0;
					int predictedLabelWithThreshold = oobEstimates[i] > candidateThresholds[k] ? 1 : 0;
					err += abs(binaryLabel-predictedLabelWithThreshold);
				}
				
				if( err < minError )
				{
					minError = err;
					bestThreshold = candidateThresholds[k];
				}
			}
			
			thresholds_.push_back(bestThreshold);
			
		}
		
		cout << "** thresholds **" << endl;
		for(int i = 0; i < thresholds_.size(); ++i)
			cout << i << "\t" << thresholds_[i] << endl;
	}
	
	~pBinaryRelevance()
	{
		for(int i = 0; i < rfs_.size(); ++i)
			delete rfs_[i];
	}
	
	vector<float> getClassScores(vector<float>& fv)
	{
		vector<float> scores;
		for(int j = 0; j < numClasses_; ++j)
			scores.push_back(rfs_[j]->estimateClassProbabilities(fv)[1]);
		return scores;
	}
	
	// TODO: it is often desirable to get both the scores and the class prediction at the same time (to avoid double calculation)

	vector<int> classify(vector<float>& fv)
	{
		vector<int> Y;
		
		for(int j = 0; j < numClasses_; ++j)
		{
			float pJ = rfs_[j]->estimateClassProbabilities(fv)[1];
			if( pJ > thresholds_[j] ) Y.push_back(j);
		}
		
		return Y;
	}
};