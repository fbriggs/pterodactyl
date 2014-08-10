// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#pragma once

#include "pImageBytes.h"

class pNoiseReducer
{
public:
	static pImageBytes whiteningFilter(pImageBytes& src)
	{
        pImageBytes newImage(src.w_, src.h_, 1);
        
		// enrgSum  - vector of float vectors containing the energy envelope for each column
		// totalCol - the number of pixels in each column
		vector< vector<float> > enrgSum;
		float totalCol = src.h_;
        
		// for each column create and push onto engrSum a float vector that contains:
		// (0) sumCol - the sum of the column divided by the # of pixels in the column
		// (1) i - the column that the sum is from
		for(int i = 0; i < src.w_; ++i)
		{
			float sum = 0;
			vector<float> sumCol;
			for(int j = 0; j < src.h_; ++j)
				sum += src.getPixelF(i,j,0)*src.getPixelF(i, j, 0);
			sumCol.push_back(sum/totalCol);
			sumCol.push_back(i);
			enrgSum.push_back(sumCol);
		}
		// sort enrgSum by the sumCol
		sort(enrgSum.begin(), enrgSum.end(), sortFloatVector);
		for(int i = 0; i < src.h_; ++i)
		{
			// for each row take the sum of the pixels in the lowest enrgPrcnt(20%) of sumCol
			float sumRow = 0;
			float enrgPrcnt = .2;
			for(int j = 0; j < enrgSum.size()*enrgPrcnt; ++j)
				sumRow += src.getPixelF((int)enrgSum[j][1], i, 0)*src.getPixelF((int)enrgSum[j][1], i, 0);
			// add to keep from dividing by zero later
			sumRow += .00000000001;
            
			// for each pixel set the value to the original pixel divided by the sumRow found above
			for(int k = 0; k < src.w_; ++k)
			{
				float pixel = src.getPixelF(k, i, 0);
				pixel /= sqrt(sumRow);
				pixel = min(pixel, 1.0f);
				newImage.setPixelF(k, i, 0, pixel);
			}
		}
		// return the new Image with all of the new pixels set
		return newImage;
	}
    
    //sorts a vector of floats by there 0ith element from low to high
	static bool sortFloatVector(vector<float> i, vector<float> j)
	{
		return (i[0] < j[0]);
	}
	
};