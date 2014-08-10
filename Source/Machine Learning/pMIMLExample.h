// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// this class defines a bag of instances with multiple labels (i.e the object of interest in a MIML dataset), as well as some 
// possible other meta data such as its fold in cross-validation, or a bag-level feature vector

#pragma once

#include <vector>
using namespace std;

class pMIMLExample
{
public:
	vector< vector<float> > instances_; // instance features vectors. first index is instance, second is feature
	vector<int> labels_;				// this is the set of labels for the bag
	vector<float> bagFeature_;			// this is scratch space for use by MIML algorithms that associate features with bags
	vector<int> instanceLabels_;		// given n instances, this is a vector of n class labels. won't necessarily be available for all datasets
	
	int bag_id_;						// this mostly here to verify the correctness of the parser
	int fold_;							// this is used in cross-validation
	
	pMIMLExample() {}
	
	// create a bag with unkown instance labels
	pMIMLExample(vector< vector<float> >& instances, vector<int>& labels)
		: instances_(instances), labels_(labels) {}
	
	// create a bag with known instance labels
	pMIMLExample(vector< vector<float> >& instances, vector<int>& labels, vector<int>& instanceLabels)
		: instances_(instances), labels_(labels), instanceLabels_(instanceLabels) {}
};


