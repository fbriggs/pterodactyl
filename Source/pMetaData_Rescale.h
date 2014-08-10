// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// functions to rescale feature vectors for use with distance-based methods (e.g. k-means, SVM)

#pragma once

#include "pMetaData.h"
#include <vector>

using namespace std;

class pMetaData_Rescale
{
public:
	// takes a reference to the audio-file meta data with segment features already loaded.
	// rescales each feature independently to the range [0,1]
	static void rescaleTo01(vector<pAudioFile_MetaData>& afmd, const int featureDim)
	{		
		cout << "feature dim = " << featureDim << endl;
		
		for(int f = 0; f < featureDim; ++f)
		{
			float minVal = FLT_MAX;
			float maxVal = -FLT_MAX;
			
			for(int i = 0; i < afmd.size(); ++i)
			for(int j = 0; j < afmd[i].segments_.size(); ++j)
			{
				minVal = min(minVal, afmd[i].segments_[j].featureVector_[f]);
				maxVal = max(maxVal, afmd[i].segments_[j].featureVector_[f]);
			}
			
			assert(maxVal != minVal); // this will cause divide by 0, and indicates a feature with all the same value which is useless
			
			for(int i = 0; i < afmd.size(); ++i)
			for(int j = 0; j < afmd[i].segments_.size(); ++j)
				afmd[i].segments_[j].featureVector_[f] = (afmd[i].segments_[j].featureVector_[f] - minVal) / (maxVal - minVal);
		}
	}
};