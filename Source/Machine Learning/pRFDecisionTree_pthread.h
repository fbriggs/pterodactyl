// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// the decision tree used in a random forest
#pragma once

#include "pExample.h"
#include "pVectorUtils.h"
#include "pRandom.h"
#include "pStringUtils.h"

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

// used to sort pExamples on a particular feature
class pExampleComparator_pRFDecisionTree_pthread
{
public:
	pExampleComparator_pRFDecisionTree_pthread(int indexToCompare) : indexToCompare_(indexToCompare) {}
	
	int indexToCompare_;
	
	bool operator ()(const pExample* a, const pExample* b) 
	{
		return a->featureVector_[indexToCompare_] < b->featureVector_[indexToCompare_];
	}
};

// this class implements the CART algorithm for decision tree learning, with minor modifications for use with Random Forest
// multiple classes are supported, and all features are assumed to be real/float valued
class pRFDecisionTree_pthread
{
public:
	bool isLeaf_;								// nodes are either internal (i.e. test a condition such as x10 < 5.3), or leaves (which are labeled with a class)

	vector<float> labelHistogram_;				// if this node is a leaf, it stores a histogram over classes
	
	int testFeature_;							// if this node is internal, this is the feature being tested
	float testThreshold_;						// if this node is internal, this is the threshold against which the tested feature is compared
	pRFDecisionTree_pthread* lessSubtree_;		// if this node is internal, this is the subtree for examples with attribute less than the test threshold
	pRFDecisionTree_pthread* moreSubtree_;		// if this node is internal, this is the subtree for examples with attribute greater than or equal to the test threshold

	// grow a tree from a labeled training data set, following the CART algorithm, with Gini impurity, and the following modification for Random Forest:
	// rather than chosing the best split variable amongst all of them, select a random subset. numFeaturesToTry is the number of features in this subset.
	// extending the basic random forest slightly, this implementation has an option to put a max depth on the tree (in which case it is much more likely that leaves will not be pure).
	// basic RF stores the single majority class label in leaves. this implementation has the option to store histograms of labels in the leaves instead, which some papers suggest
	// improves performance (controlled by storeHistogramInLeaf).
	pRFDecisionTree_pthread(vector<pExample*>& examples, int numClasses, int featureDim, int numFeaturesToTry, bool storeHistogramInLeaf, int depth, int maxDepth)
	{		
		// this will be used as the threshold difference for determining equality
		// between floats. without this measure, bogus results are possible due to a failed
		// equality test betweeen numbers that are basically equal
		float epsilon = 10E-7; 
		
		assert( examples.size() != 0 );
		assert( numFeaturesToTry <= featureDim );
		
		lessSubtree_ = NULL;
		moreSubtree_ = NULL;
		
		// if all examples are the same class, make this node a leaf with their label
		if( depth >= maxDepth || allSameClass(examples) )
		{
			makeLeaf(examples, numClasses, storeHistogramInLeaf);
			return;
		}
		
		// we know its not a leaf now... (unless all of the features we select randomly end up having identical values. this case is handled at the end)
		isLeaf_ = false;
		
		// pick a random subset of features to evaluate
		vector<int> allFeatures;
		for(int i = 0; i < featureDim; ++i) allFeatures.push_back(i);
		
		pRandom::shuffle(allFeatures);
		vector<int> subsetOfFeatures;
		for(int i = 0; i < numFeaturesToTry; ++i) subsetOfFeatures.push_back(allFeatures[i]);
		
		// keep track of the best feature and split threshold evaluated so far
		float lowestGini = FLT_MAX;
		bool foundASplit = false; // it is possible that due to redundant values, there are no possible splits... 
		int splitIndexThatWasFound = -1;
		
		int numExamples = examples.size();
		
		// find the optimal feature and value on which to split
		for(int i = 0; i < subsetOfFeatures.size(); ++i)
		{
			int f = subsetOfFeatures[i]; // the current feature to try
			 
			sort(examples.begin(), examples.end(), pExampleComparator_pRFDecisionTree_pthread(f));
			
			vector<int> classCountLess(numClasses); // for the current split index, this is the number of examples that fall in the 'less' partion, for each class
			vector<int> classCountMore(numClasses); // same as above, but for the 'more' partition
			
			
			int numLess = 0;
			int numMore = numExamples;

			// initially, only thet first example is on the less side of the partition
			//++classCountLess[examples[0]->classLabel_];
			
			// before we start evaluating splits, count the number of examples on the 'more' side of the partition (which is everything)
			for(int ex = 0; ex < numExamples; ++ex)
				++classCountMore[examples[ex]->classLabel_];
			
			// consider all of the possible values on which to split for this variable
			// there are at most # of examples - 1 such values. we only need to test
			// at places where the class label switches, since
			// the optimal cut point must be the midpoint between two such examples
			for(int ex = 1; ex < numExamples; ++ex)
			{
				// shift the example on the more side of the split we just tested over to the less side
				--classCountMore[examples[ex-1]->classLabel_];
				--numMore;
				++classCountLess[examples[ex-1]->classLabel_];
				++numLess;
				
				// check if the classes differ between this example and the previous one.
				// if not, we don't need to consider a split between them.
				// if two consective (sorted) examples have the same feature value, there is no point in trying to split
				// between them either.
		
				if( examples[ex - 1]->classLabel_ == examples[ex]->classLabel_ ) continue;
				if( fabs(examples[ex - 1]->featureVector_[f] - examples[ex]->featureVector_[f]) < epsilon ) continue;
				
				// otherwise, find the midpoint between the two examples for the current variable.
								
				// evaluate the quality of this split based on gini index
				float giniLess = gini(classCountLess, numLess, numClasses);
				float giniMore = gini(classCountMore, numMore, numClasses);

				float fracLess = (float)numLess / (float) numExamples;
				float fracMore = (float)numMore / (float) numExamples;
				
				float weightedGini = giniLess * fracLess + giniMore * fracMore;
				
				foundASplit = true;
				
				// if this new split is the best so far, keep it
				if( weightedGini < lowestGini )
				{
					splitIndexThatWasFound = ex;
					
					
					lowestGini = weightedGini;
					testFeature_ = f;
					testThreshold_ = (examples[ex - 1]->featureVector_[f] + examples[ex]->featureVector_[f]) / 2.0f;
					
					assert(numLess != 0);
					assert(numMore != 0);
				}	
				
				
			}
		}
		
		// it is possible that we don't find a split that makes any improvement (i.e. all instances go in one child and none go in the other).
		// this can happen if the randomly chosen features to split on are all identical. in this case, we will declare this node a leaf
		// and set its class to the majority from the examples (with a random tie break).
		if( !foundASplit )
		{
			makeLeaf(examples, numClasses, storeHistogramInLeaf);
			return;
		}

		// now that we have the variable on which to split, recursively build the subtrees
		vector<pExample*> examplesLess;
		vector<pExample*> examplesMore;
		splitExamples(examples, testFeature_, testThreshold_, examplesLess, examplesMore);

		lessSubtree_ = new pRFDecisionTree_pthread(examplesLess, numClasses, featureDim, numFeaturesToTry, storeHistogramInLeaf, depth + 1, maxDepth);
		moreSubtree_ = new pRFDecisionTree_pthread(examplesMore, numClasses, featureDim, numFeaturesToTry, storeHistogramInLeaf, depth + 1, maxDepth);
	}
	
	void makeLeaf(vector<pExample*>& examples, int numClasses, bool storeHistogramInLeaf)
	{
		isLeaf_ = true;
		labelHistogram_ = classHistogram(examples, numClasses);
		if( storeHistogramInLeaf) makeHistogramIntoMajority(labelHistogram_);
	}
	
	// test if all of the labels passed in are the same
	static bool allSameClass(vector<pExample*>& examples)
	{
		for(int i = 0; i < examples.size(); ++i)
			if( examples[i]->classLabel_ != examples[0]->classLabel_ ) return false;
		return true;
	}
	
	// split the examples into outExamplesLess and outExamplesGreaterEq by comparing them on splitFeature against splitThreshold
	static void splitExamples(vector<pExample*>& examples, int splitFeature, float splitThreshold, vector<pExample*>& outExamplesLess, vector<pExample*>& outExamplesGreaterEq)
	{
		int numExamples = examples.size();
		for(int i = 0; i < numExamples; ++i)
		{
			if( examples[i]->featureVector_[splitFeature] < splitThreshold )
				outExamplesLess.push_back(examples[i]);
			else
				outExamplesGreaterEq.push_back(examples[i]);
		}
	}
	
	// measure the gini impurity in the class labels for a collection of examples
	// profiling indicates that a lot of time is spent in this function, so it is optimized
	static float gini(vector<int>& classCounts, int numExamples, int numClasses)
	{
		int i;
		float sum = 1;
		float oneOverNumExamples = 1.0 / (float) numExamples;
		float oneOverNumExamplesSquared = oneOverNumExamples * oneOverNumExamples;
		for(i = 0; i < numClasses; ++i) sum -= classCounts[i]  * classCounts[i] * oneOverNumExamplesSquared;
		return sum;
	}
	
	static vector<float> classHistogram(vector<pExample*>& examples, int numClasses)
	{
		vector<float> classCounts(numClasses);
		for(int i = 0; i < examples.size(); ++i)
			classCounts[examples[i]->classLabel_] += 1; // possible improvement: add a small random number here to cause ties to be broken randomly for majority
		
		pVectorUtils::normalizeVectorToPDF(classCounts);
		
		return classCounts;
	}
		
	// return the class label for feature vector as determined by this decision tree
	int classify(vector<float>& featureVector)
	{
		if( isLeaf_ ) return pVectorUtils::argmax(labelHistogram_);
		
		if( featureVector[testFeature_] < testThreshold_ )
			return lessSubtree_->classify(featureVector);
		else
			return moreSubtree_->classify(featureVector);
	}
	
	// given an input feature vector, find the corresponding leaf
	// and return the class histogram it stores
	vector<float> getClassHistogramForInput(vector<float>& featureVector)
	{
		if( isLeaf_ ) return labelHistogram_;
		
		if( featureVector[testFeature_] < testThreshold_ )
			return lessSubtree_->getClassHistogramForInput(featureVector);
		else
			return moreSubtree_->getClassHistogramForInput(featureVector);
	}
	
	// take a histogram, find the majority class, make that probability 1 and all others 0
	void makeHistogramIntoMajority(vector<float>& hist)
	{
		int majority = pVectorUtils::argmax(hist);
		for(int i = 0; i < hist.size(); ++i)
			hist[i] = (i == majority) ? 1 : 0;
	}
	
	// deconstructor recusrively deletes children
	~pRFDecisionTree_pthread()
	{
		if( lessSubtree_ != NULL ) delete lessSubtree_;
		if( moreSubtree_ != NULL ) delete moreSubtree_;
	}

	//// file io ////
	
	// return a string representation of this tree
	string toString()
	{
		if( isLeaf_ )
		{
			//return "[" + pStringUtils::intToString(classLabel_) + "]";
			
			string s = "[";
			for(int i = 0; i < labelHistogram_.size(); ++i)
				s += pStringUtils::floatToString(labelHistogram_[i]) + (i == labelHistogram_.size() - 1 ? "]" : ",");
			
			return s;
		}
		else 
		{
			return "(" + pStringUtils::intToString(testFeature_) + "," + pStringUtils::floatToString(testThreshold_) + "," + lessSubtree_->toString() + "," + moreSubtree_->toString() + ")";
		}
	}
	
	// parse a decision tree from a saved string
	pRFDecisionTree_pthread(string src, int numClasses)
	{
		if( src[0] == '[' ) // leaf
		{
			vector<string> parts = pStringUtils::splitNonEmpty(src, "[],");
			
			if( parts.size() != numClasses )
			{
				cout << "PARSER ERROR" << endl;
				cout << "src=" << src << endl;
				pVectorUtils::cout_vector(parts, "^^|^^");
				cout << endl;
				exit(0);
			}
			
			assert(parts.size() == numClasses);
			
			for(int i = 0; i < parts.size(); ++i)
				labelHistogram_.push_back(pStringUtils::stringToFloat(parts[i]));
			
			isLeaf_ = true;
			lessSubtree_ = NULL;
			moreSubtree_ = NULL;
		}
		else if( src[0] == '(' ) // non-leaf
		{
			vector<string> parts = pStringUtils::splitNonEmpty(src, "(),");
			testFeature_ = pStringUtils::stringToInt(parts[0]);
			testThreshold_ = pStringUtils::stringToFloat(parts[1]);
			
			// find the second comma, which is where the lessSubtree starts
			int secondCommaIndex = 0;
			int commasFound = 0;
			while( secondCommaIndex < src.size() && commasFound < 2)
			{
				if( src[secondCommaIndex] == ',' )
					commasFound++;
				
				++secondCommaIndex;
			}
			
			// find the next comma at the same paren level, which is the end of the
			// first subtree and the start of the second
			int thirdCommaIndex = secondCommaIndex;
			int indentLevel = 0;
			for(; thirdCommaIndex < src.size(); ++thirdCommaIndex)
			{
				if( src[thirdCommaIndex] == ',' && indentLevel == 0 )
					break;
				
				if( src[thirdCommaIndex] == '(' || src[thirdCommaIndex] == '[' )
					++indentLevel;
				
				if( src[thirdCommaIndex] == ')' || src[thirdCommaIndex] == ']' )
					--indentLevel;
			}
			
			string lessStr = src.substr(secondCommaIndex, thirdCommaIndex - secondCommaIndex);
			string moreStr = src.substr(thirdCommaIndex + 1, src.size() - thirdCommaIndex - 2);
			
			isLeaf_ = false;
			lessSubtree_ = new pRFDecisionTree_pthread(lessStr, numClasses);
			moreSubtree_ = new pRFDecisionTree_pthread(moreStr, numClasses);
		}
		else assert(false); // parse error
	}
};