// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once

#include "pImageBytes.h"
#include "pVectorUtils.h"

class pFastGaussianBlur
{
public:

	// a faster gaussian blur that uses the separability of the gaussian kernel
	// assumes a 1 channel image
	// applies a convolution over a window of pixels that is 2 * radius + 1 in width
	// with sigma = radius
	static pImageBytes blur(pImageBytes& src, int radius)
	{
		//assert(src.channels_ == 1);
		
		float sigma = radius;
		vector<float> kernel;
		for(int u = -radius; u <= radius; ++u)
			kernel.push_back(exp(-0.5 * float(u * u) / float(sigma * sigma)));
		
		assert(kernel.size() == 2 * radius + 1);
		
		pVectorUtils::normalizeVectorToPDF(kernel);
	
		pImageBytes blurred(src.w_, src.h_, src.channels_);
		pImageBytes blurred2(src.w_, src.h_, src.channels_);
		
		for(int c = 0; c < src.channels_; ++c)
		{
			for(int x = 0; x < src.w_; ++x)
			{
				for(int y = 0; y < src.h_; ++y)
				{
					float conv = 0;
					
					for(int v = y - radius, k = 0; v <= y + radius; ++v, ++k)
						conv += kernel[k] * src.getPixelF_extend(x, v, c);
			
					blurred.setPixelF(x, y, c, conv);
				}
			}
			
			
			for(int x = 0; x < src.w_; ++x)
			{
				for(int y = 0; y < src.h_; ++y)
				{
					float conv = 0;
					
					for(int u = x - radius, k = 0; u <= x + radius; ++u, ++k)
						conv += kernel[k] * blurred.getPixelF_extend(u, y, c);
					
					blurred2.setPixelF(x, y, c, conv);
				}
			}
		}
		
		
		return blurred2;
	}
};