// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// this file defines image convolution and correlation functions
#pragma once

#include "pImageBytes.h"
#include "pImageUtils.h"
#include "pMatrix.h"

// this class defines a method for computing correlation between images
class pCorrelation
{
public:
	// compute the correlation between two images. the output
	// will be an image with the same dimensions as srcA
	static pImageBytes twoImageCorrelation(pImageBytes& srcA, pImageBytes& srcB)
	{
		assert(srcA.channels_ == 1 && srcB.channels_ == 1); // only for 1-channel images
		
		int halfBWidth  = srcB.w_ / 2;
		int halfBHeight = srcB.h_ / 2;
		
		// make matrix to hold the non-normalized values
		pMatrix corr(srcA.w_, srcA.h_);
		
		// compute the correlation at each pixel
		for(int x = 0; x < srcA.w_; ++x)
		{
			cout << "*"; cout.flush(); // progress bar
			for(int y = 0; y < srcA.h_; ++y)
			{
				for(int i = -halfBWidth;  i <= halfBWidth;  ++i)
				for(int j = -halfBHeight; j <= halfBHeight; ++j)
					corr.data_[x][y] += (srcA.getPixel_extend(x + i, y + j, 0) / 255.0f)  * 
										(srcB.getPixel(i + halfBWidth, j + halfBHeight, 0) / 255.0f);
			}
		}
		
		// normalize the results
		corr.normalizeTo01();
		
		pImageBytes output(srcA.w_, srcA.h_, 1);
		for(int x = 0; x < srcA.w_; ++x)
		for(int y = 0; y < srcA.h_; ++y)
			output.setPixelF(x, y, 0, corr(x, y));
		
		return output;
	}
};

// this class defines functions for image convolutions
class pConvolution
{
public:
	static pMatrix kernel_avg(int size)
	{
		pMatrix k(size, size);
		for(int i = 0; i < size; ++i)
		for(int j = 0; j < size; ++j)
			k.data_[i][j] = 1.0 / (float(size) * float(size));
		return k;
	}
	
	// returns the variation of the laplacian kernel with 8 at the center
	static pMatrix kernel_laplacian8()
	{
		return pMatrix(	-1, -1, -1,
						-1,  8, -1,
						-1, -1, -1);
	}
	
	// returns the variation of the laplacian kernel with -8 at the center
	static pMatrix kernel_laplacianNeg8()
	{
		return pMatrix(	1,  1, 1,
					    1, -8, 1,
					    1,  1, 1);
	}
	
	// returns the variation of the laplacian kernel with 4 at the center
	static pMatrix kernel_laplacian4()
	{
		return pMatrix(	0, -1,  0,
					   -1,  4, -1,
					    0, -1,  0);
	}
	
	// compute N(x,y) = exp(-(x^2 + y^2)/2 sigma^2)
	static float normal(float x, float y, float sigma2) { return exp(- (x*x + y*y) / (2 * sigma2)); }
	
	// return a 3x3 gaussian kernel with parameter sigma
	static pMatrix kernel_gaussian(float sigma2)
	{
		return pMatrix(	normal(-1,	-1, sigma2),		normal(0,	-1,	sigma2),	normal(1,	-1,	sigma2), 
						normal(-1,	0,	sigma2),		normal(0,	0,	sigma2),	normal(1,	0,	sigma2), 
						normal(-1,	1,	sigma2),		normal(0,	1,	sigma2),	normal(1,	1,	sigma2)).normalizeSumTo1();
	}
	
	static pMatrix kernel_gaussian(int size, float sigma2)
	{
		int halfWidth = size / 2;
		int halfHeight = size / 2;
		
		pMatrix k(size, size);
		for(int i = -halfWidth; i <= halfWidth; ++i)
		for(int j = -halfHeight; j <= halfHeight; ++j)
			k.data_[i + halfWidth][j + halfHeight] = normal(i,	j,	sigma2);
		
		k.normalizeSumTo1();
		
		return k;
	}
	
	// apply a convolution to img at point (x,y) in channel c, using the supplied kernel 
	static float applyConvolutionToPixel(pImageBytes& img, pMatrix& kernel, int x, int y, int c = 0)
	{
		float sum = 0;
		int halfWidth = kernel.cols_ / 2;
		int halfHeight = kernel.rows_ / 2;
		
		for(int i = -halfWidth; i <= halfWidth; ++i)
		for(int j = -halfHeight; j <= halfHeight; ++j)
			sum += kernel(i + halfWidth, j + halfHeight) * img.getPixel_extend(x - i, y - j, c);
		
		return pImageUtils::clamp255(sum);
	}
	
	// apply the supplied convolution kernel to all pixels in img and return a new image
	static pImageBytes applyConvolutionToAllPixels(pImageBytes& img, pMatrix& kernel)
	{
		pImageBytes output(img.w_, img.h_, img.channels_);
		
		for(int c = 0; c < img.channels_; ++c)
		for(int x = 0; x < img.w_; ++x)
		for(int y = 0; y < img.h_; ++y)
			output.setPixel(x,y,c, applyConvolutionToPixel(img, kernel, x, y, c));
		
		return output;
	}
};