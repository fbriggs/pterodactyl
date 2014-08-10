// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// extensions of the basic random forest related to out-of-bag (OOB) estimation. 
// OOB estimation can be used to get an estimate of the classifier's accuracy without cross-validation,
// and to make predictions about data in the training set without overfitting

#pragma once

#include "pRandomForest_pthread.h"

class pRandomForest_OutOfBag
{
public:
	// use the random forest to classify one of the examples in its training set 
	// (taking votes only from trees for which the example is out-of-bag).
	// feature - the feature vector for the example
	// index - the index of the example in the training dataset
	// returns a histogram of class probabilities
	static vector<float> estimateClassProbabilitiesOOB(pRandomForest_pthread& rf, vector<float>& feature, int index)
	{
		assert(feature.size() == rf.featureDim_ ); 
		
		vector<float> outClassProbabilities(rf.numClasses_);
		
		for(int i = 0; i < rf.trees_.size(); ++i)
		{
			if( rf.inBag_[i][index] ) continue;
			
			vector<float> histo = rf.trees_[i]->getClassHistogramForInput(feature);
			
			for(int c = 0; c < rf.numClasses_; ++c)
				outClassProbabilities[c] += histo[c];
		}
		
		pVectorUtils::normalizeVectorToPDF(outClassProbabilities);
		
		return outClassProbabilities;
	}
};