// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// functions for computing error/performance measures of a multi-label classifier,
// e.g., Hamming loss, rank-loss, one-error, coverage, set-level 0/1 loss, etc. (not all of these are implemented right now, but thats the kind of thing that would be in here). 

#pragma once

#include <vector>
#include "pVectorUtils.h"

using namespace std;

class pMultiLabelErrorMeasure
{
public:
	// input Y- a list of classes, output - a binary vector
	static vector<int> labelSetToBinaryVector(vector<int>& Y, int numClasses)
	{
		vector<int> V(numClasses, 0);
		
		for(int i = 0; i < Y.size(); ++i)
			V[Y[i]] = 1;
			
		return V;
	}

	// compares two label sets Y and Z. returns the number of differences in the binary representations of Y and Z, divided by the number of classes
	static float hammingLoss(vector<int>& Y, vector<int>& Z, int numClasses)
	{	
		vector<int> vY = labelSetToBinaryVector(Y, numClasses);
		vector<int> vZ = labelSetToBinaryVector(Z, numClasses);
		
		float err = 0;
		for(int i = 0; i < numClasses; ++i)
			err += abs(vY[i] - vZ[i]);
			
		return err/numClasses;
	}
	
	// returns 0 if Y = Z and 1 otherwise
	static float set01loss(vector<int>& Y, vector<int>& Z, int numClasses)
	{	
		vector<int> vY = labelSetToBinaryVector(Y, numClasses);
		vector<int> vZ = labelSetToBinaryVector(Z, numClasses);
		
		for(int i = 0; i < numClasses; ++i)
			if( vZ[i] != vY[i] ) return 1;
		return 0;
	}
	
	// updates _TP, _TN, _FP, _FN with true/false positives/negatives. 
	static void true_false_pos_neg(vector<int>& prediction, vector<int>& groundTruth, int numClasses, int& _TP, int& _TN, int& _FP, int& _FN)
	{
		vector<int> vGroundTruth = labelSetToBinaryVector(groundTruth, numClasses);
		vector<int> vPrediction = labelSetToBinaryVector(prediction, numClasses);
		
		for(int j = 0; j < numClasses; ++j)
		{
			if( vGroundTruth[j] == 1 && vPrediction[j] == 1 ) ++_TP;
			if( vGroundTruth[j] == 1 && vPrediction[j] == 0 ) ++_FN;
			if( vGroundTruth[j] == 0 && vPrediction[j] == 1 ) ++_FP;
			if( vGroundTruth[j] == 0 && vPrediction[j] == 0 ) ++_TN;
		}
	}
	
	// same version of rank-loss defined in "Acoustic classification of multiple simultaneous bird species: A multi-instance multi-label approach"
	static float rankLoss(vector<float>& scores, vector<int>& Y, int numClasses)
	{
		float loss = 0;
		
		for(int j = 0; j < numClasses; ++j)
		for(int k = 0; k < numClasses; ++k)
		{
			if( pVectorUtils::vector_contains(Y, j) && !pVectorUtils::vector_contains(Y, k) )
			{
				if( scores[j] <= scores[k] )
					loss += 1;
			}
		}
		return loss / (Y.size() * (numClasses - Y.size()));;
	}
	
	// one error is the fraction of bags for which the top-scoring label is not in the ground-truth label set
	static float oneError(vector<float>& scores, vector<int>& Y)
	{
		int topScoringLabel = pVectorUtils::argmax(scores);
	
		return pVectorUtils::vector_contains(Y, topScoringLabel) ? 0 : 1;
	}
	
	// the class-scores can be ranked from highest to lowest. coverage measures how far down on that
	// ranking we have to go before getting all of the true classes
	static float coverage(vector<float>& scores, vector<int>& Y, int numClasses)
	{
		vector< pair<float,int> > scoreClassPairs;
		for(int j = 0; j < numClasses; ++j)
			scoreClassPairs.push_back(pair<float,int>(scores[j], j));
			
		sort(scoreClassPairs.begin(), scoreClassPairs.end());
		reverse(scoreClassPairs.begin(), scoreClassPairs.end());
		
		// go down the list of classes ordered by score (descending)
		vector<bool> covered(numClasses);
		for(int j = 0; j < numClasses; ++j)
		{
			covered[scoreClassPairs[j].second] = true;
			
			// check if all of the classes in the label set are covered
			bool allCovered = true;
			for(int k = 0; k < Y.size(); ++k)
			{
				if( !covered[Y[k]] )
					allCovered = false;
			}
			
			if( allCovered )
				return j;
		}
		
		assert(false); // this is just to stop a warning; the logic of this function should prevent it from getting here
	}
};