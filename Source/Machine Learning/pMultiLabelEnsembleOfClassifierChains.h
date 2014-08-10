// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// this is a basic implementation of the 'ensemble of probabalistic classifier chains' algorithm, using random forest as base-classifier
// this means it is an ensemble of classifier chains, where each base-classifier outputs a probability. each chain in the ensemble estimates the
// probability of each label to be present given an input feature vector. the probability estimates of all chains in the ensemble are averaged
// and compared to a threshold of 0.5 to produce the final predicted label set.
#pragma once

#include "pRandom.h"
#include "pExample.h"
#include "pExample_MultiLabel.h"
#include "pRandomForest_pthread.h"
#include "pRandomForest_OutOfBag.h"

// this struct stores the parameteres of a pMultiLabelEnsembleOfClassifierChains relevance classifier, in case it is necessary to pass them around
class pMultiLabelEnsembleOfClassifierChains_Params 
{
public:
	int numTreesPerChain_;
	int numChains_;
	bool useOutOfBagCorrection_;
	
	pMultiLabelEnsembleOfClassifierChains_Params(int numTreesPerChain, int numChains, bool useOutOfBagCorrection)
		: numTreesPerChain_(numTreesPerChain), numChains_(numChains), useOutOfBagCorrection_(useOutOfBagCorrection) {}
};

// a chain a sequence of c classifiers, where c is the number of classes
struct pMLC_Chain
{
	vector<int> pi_;								// pi is short for permutation. pi_[j] is the j'th class in this chain's order
	vector<pRandomForest_pthread*> rfs_;			// rfs_[j] is a binary classifier which makes a decision about the class pi_[j]
};

// this is the ensemble. it consists of a collection of pMLC_Chain's
class pMultiLabelEnsembleOfClassifierChains
{
public:
	int c_;
	int numChains_;
	vector<pMLC_Chain> chains_;
	vector<float> thresholds_; // one threshold for each class
	
	void train(vector<pExample_MultiLabel>& exs, int numClasses, int numTreesPerChain, int numChains, bool useOutOfBagCorrection)
	{
		cout << "pMultiLabelEnsembleOfClassifierChains" << endl;
		cout << "n = " << exs.size() << endl;
		cout << "c = " << numClasses << endl;
		cout << "T/chain = " << numTreesPerChain << endl;
		cout << "# chains =" << numChains << endl;
		cout << "OOB = " << (useOutOfBagCorrection ? "true" : "false") << endl; 
		
		c_ = numClasses;
		numChains_ = numChains;
		
		for(int cn = 0; cn < numChains; ++cn)
		{
			cout << "building chain " << cn << endl;
			
			chains_.push_back(pMLC_Chain());
			pMLC_Chain& chain = chains_[cn];
			
			for(int i = 0; i < c_; ++i) chain.pi_.push_back(i);
			pRandom::shuffle(chain.pi_);
			
			// the first classifier in the chain gets a SISL dataset consisting of the feature vector paired
			// with the label that is first in the chain's order.
			vector<pExample> sislExamples;
			for(int i = 0; i < exs.size(); ++i) 
				sislExamples.push_back(pExample(exs[i].fv_)); 
			
			// in each iteration of the loop below, we will maintain the invariant that sislExamples stores the feature vectors augmented with the correct number of bits of the (maybe estimated) label set
			for(int j = 0; j < c_; ++j)
			{
				if( j > 0 ) // if this is not the first classifier in the chain, add one more element to its feature vector, which is a bit of the label set (or an OOB prediction of it)
				{
					for(int i = 0; i < exs.size(); ++i)
					{
						float psuedoLabel;
						if( useOutOfBagCorrection )	psuedoLabel = pVectorUtils::vector_contains(exs[i].labels_, chain.pi_[j-1]) ? 1 : 0;
						else						psuedoLabel = pRandomForest_OutOfBag::estimateClassProbabilitiesOOB(*chain.rfs_[j-1], sislExamples[i].featureVector_, i)[1];
						
						//cout << "psuedo label = " << (pVectorUtils::vector_contains(exs[i].labels_, chain.pi_[j-1]) ? 1 : 0) << " oob psuedo label = " << (pRandomForest_OutOfBag::estimateClassProbabilitiesOOB(*chain.rfs_[j-1], sislExamples[i].featureVector_, i)[1]) << endl;
						
						sislExamples[i].featureVector_.push_back(psuedoLabel);
					}
				}
				
				// update the target label for the SISL examples for the current class j
				for(int i = 0; i < exs.size(); ++i) 
					sislExamples[i].classLabel_ = pVectorUtils::vector_contains(exs[i].labels_, chain.pi_[j]) ? 1 : 0;
				
				chain.rfs_.push_back(new pRandomForest_pthread);
				
				bool storeHistogramInLeaf = true;
				int maxDepth = 15;
				bool useClassBalancedBootstrapping = false;
				chain.rfs_[j]->train_multithread(numTreesPerChain, 2, sislExamples, storeHistogramInLeaf, maxDepth, useClassBalancedBootstrapping);
			}
		}
		
		/// find optimal thresholds for each class //
		cout << "computing OOB score vectors for threshold calibration" << endl;
		vector< vector<float> > oobScoreVectors;
		for(int i = 0; i < exs.size(); ++i)
			oobScoreVectors.push_back(scoreOOB(exs[i].fv_, i));
			
		for(int j = 0; j < numClasses; ++j)
		{
			cout << "computing threshold for class " << j << endl;
			
			vector<float> candidateThresholds;
			for(float thresh = 0.001; thresh < 1.0f; thresh += 0.001) candidateThresholds.push_back(thresh);
			
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
					int predictedLabelWithThreshold = oobScoreVectors[i][j] > candidateThresholds[k] ? 1 : 0;
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
	
	// produce a vector of class-scores using OOB estimate for the indexed instance
	vector<float> scoreOOB(vector<float>& fv, int index)
	{
		vector<float> votes(c_, 0);
		
		for(int cn = 0; cn < numChains_; ++cn)
		{
			vector<float> fvThisChain(fv);
			
			for(int j = 0; j < c_; ++j)
			{
				float p = pRandomForest_OutOfBag::estimateClassProbabilitiesOOB(*chains_[cn].rfs_[j], fvThisChain, index)[1];
				votes[chains_[cn].pi_[j]] += p;
				fvThisChain.push_back(p);
			}
		}
		
		for(int j = 0; j < c_; ++j) votes[j] /= float(numChains_);
			
		return votes;
	}
	
		
	vector<float> getClassScores(vector<float>& fv)
	{
		vector<float> votes(c_, 0);
		
		for(int cn = 0; cn < numChains_; ++cn)
		{
			vector<float> fvThisChain(fv);
			for(int j = 0; j < c_; ++j)
			{
				float p = chains_[cn].rfs_[j]->estimateClassProbabilities(fvThisChain)[1];
				votes[chains_[cn].pi_[j]] += p;
				fvThisChain.push_back(p);
			}
		}
		
		for(int j = 0; j < c_; ++j) votes[j] /= float(numChains_);
		
		return votes;
	}
	
	vector<int> classify(vector<float>& fv)
	{
		vector<float> votes = getClassScores(fv);
		
		vector<int> Y;
		for(int j = 0; j < c_; ++j) if( votes[j] > thresholds_[j] ) Y.push_back(j);
		return Y;
	}
	
	~pMultiLabelEnsembleOfClassifierChains() 
	{
		for(int cn = 0; cn < chains_.size(); ++cn)
			for(int j = 0; j < chains_[cn].rfs_.size(); ++j) delete chains_[cn].rfs_[j]; 
	}
};