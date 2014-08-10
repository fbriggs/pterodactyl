// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// functions for re-scaling feature vectors for use with distance/kernel-based methods. 
#pragma once

class pFeatureRescale
{
public:
	static void rescaleForKMeans(vector< vector<float> >& fvs)
	{
		int featureDim = fvs[0].size();
		int numExamples = fvs.size();
		
		for(int j = 0; j < featureDim; ++j)
		{
			float fMin = FLT_MAX;
			float fMax = -FLT_MAX;
			
			for(int i = 0; i < numExamples; ++i)
			{
				fMin = min(fMin, fvs[i][j]);
				fMax = max(fMax, fvs[i][j]);
			}
			
			cout << "feature " << j << ": " << fMin << ", " << fMax << endl;
			
			float diff = fMax - fMin;
			
			assert(diff != 0); // a diff of zero here indicates a bad feature (i.e. all zero for every instance)
			
			for(int i = 0; i < numExamples; ++i)
				fvs[i][j] = (fvs[i][j] - fMin) / diff;
		}
	}
};