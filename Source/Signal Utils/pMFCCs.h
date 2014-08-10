// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


// this class extracts MFCCs from a signal frame
// the implementation follows from the description in "Comparative Evaulation of Various MFCC implementations on the Speaker Verification Task";
// of the several methods described in this paper, the one implemented here is the MFCC-FB24, based on the Hidden Markov Model Toolkit (HTK)

#pragma once

#include <iostream>
#include <math.h>
#include "pImageBytes.h"
#include "pFileIO_BMP.h"
#include <vector>
#include <assert.h>
using namespace std;

class pMFCCs
{
public:
	static const int numFilters = 24;
	static const int numMFCCs = 12;
	
	static vector< vector<float> > filterBank_;
	static int maxNonZeroFilterElement_;
	
	static void makeFilterBank(
		int numSpectrumBins,					// the number of bins in the magnitude spectrum
		float samplingFrequency,				// the nuber of samples per second of audio
		float fLow,								// minimum frequency included by filters
		float fHigh								// maximum frequency included by filters
	)	
	{		
		cout << "*** making mfcc filterbank ***" << endl;
		for(int i = 0; i < numFilters; ++i )
		{
			cout << "making filter " << i << endl;
			vector<float> filter;
			for(int k = 0; k < numSpectrumBins; ++k)
				filter.push_back(integrateTriangleFilterOverBin(i, k, numSpectrumBins, samplingFrequency, fLow, fHigh));
			filterBank_.push_back(filter);
		}
		
		// find the index of the highest frequency spectrum bin that will be given non-zero weight by any
		// filter. we can save a lot of time by only going this far in the calculation of the filters
		for(int i = 0; i < numSpectrumBins; ++i)
			if( filterBank_[numFilters - 1][i] > 0 )
				maxNonZeroFilterElement_ = i;
	
		cout << "max non-zero filter element = " << maxNonZeroFilterElement_ << endl;
	}
	
	// evaluate triangle filter i at bin index k.
	// NOTE: k is a continuous value, and the filter is not constant
	// over the range of a bin. it is necessary to integrate this function
	// to get a filter that isn't badly aliased. there is some redundant 
	// calculation here, but it doesn't matter because this is done at initialization, only once.
	static float triangleFilter(int i, float k, int numSpectrumBins, float samplingFrequency, float fLow, float fHigh)
	{
		float freqToBinCoversion = (float)numSpectrumBins / samplingFrequency;
		float fbi =		freqToBinCoversion * filterBoundary(i, fLow, fHigh);
		float fbiM1 =	freqToBinCoversion * filterBoundary(i-1, fLow, fHigh);
		float fbiP1 =	freqToBinCoversion * filterBoundary(i+1, fLow, fHigh);
		
		if( k < fbiM1 ) return 0;
		else if( fbiM1 <= k && k <= fbi ) return  ((float)k - fbiM1)/(fbi - fbiM1);
		else if( fbi <= k && k <= fbiP1 ) return  (fbiP1 - (float)k)/(fbiP1 - fbi);
		return 0;
	}
	
	// numerically integrate triangle filter i over the range k to k+1.
	// this returns the appropriate value for filter i, bin k to use in
	// the dot product with the spectrum to compute MFCs.
	static float integrateTriangleFilterOverBin(int i, float k, int numSpectrumBins, float samplingFrequency, float fLow, float fHigh)
	{
		float sum = 0;
		float deltaK = 0.01;
		for(float j = k; j < k + 1; j += deltaK )
			sum += deltaK * triangleFilter(i, j, numSpectrumBins, samplingFrequency, fLow, fHigh);
		return sum;
	}
	
	
	static void makeFilterbankImage()
	{
		pImageBytes img(numFilters, filterBank_[0].size(), 1);
		
		
		for(int i = 0; i < numFilters; ++i)
			for(int j = 0; j < filterBank_[0].size(); ++j)
				img.setPixelF(i, j, 0, filterBank_[i][j]);
		
		pFileIO_BMP::write(img, "misc_output/mfcc_filter_bank.bmp");
	}
	
	// convert a frequency in the linear scale to the mel scale
	static float fLinToMel(float fLin) { return 2595.0 * log10(1 + fLin / 700); }
	
	// convert a frequency in the mel scale to the linear scale
	static float fMelToLin(float fMel) { return 700.0 * (pow(10, fMel / 2595.0) -1.0); }
	
	// compute the value f_bi described in eqn (7)
	static float filterBoundary(int i, float fLow, float fHigh)
	{
		float deltaFMel = (fLinToMel(fHigh) - fLinToMel(fLow)) / ((float)numFilters + 0);
		float fciMel = fLinToMel(fLow) + i * deltaFMel;
		return fMelToLin(fciMel);
	}
	
	static void getFrameMFCCs(vector<float>& powerSpectrum, vector<float>& outMFCCs, vector<float>& outMFCs)
	{
		for(int i = 0; i < numFilters; ++i)
		{
			float filterResponse = computeFilterResponse(filterBank_[i], powerSpectrum, maxNonZeroFilterElement_);
			
			assert(!isinf(filterResponse));
			assert(!isnan(filterResponse));
			
			if( filterResponse < 0.001 ) filterResponse = 0.001; // THIS IS NECESSARY TO PREVENT log(0) = -inf from messing everything up
			float mfc = log10( filterResponse );
			
			assert(!isinf(mfc));
			assert(!isnan(mfc));
			
			outMFCs.push_back( mfc ); // this is the calculation in eqn (11)
		}
		
		vector<float> fullListOfMFCCs;
		discreteCosineTransform(outMFCs, fullListOfMFCCs);
		
		// take only the first several MFCCs, not all of them
		for(int i = 0; i < numMFCCs; ++i)
			outMFCCs.push_back(fullListOfMFCCs[i]);
	}
	
	static float computeFilterResponse(vector<float>& filter, vector<float>& powerSpectrum, int maxNonZeroFilterElement)
	{
		assert(powerSpectrum.size() >= maxNonZeroFilterElement);
		float sum = 0;
		for(int i = 0; i < maxNonZeroFilterElement; ++i)
		{
			sum += filter[i] * powerSpectrum[i];
		}	
		return sum;
	}
	
	// compute the DCT of the vector x, and write the results to outCoefs
	// this is O(n^2), but it is very fast for n=24 (faster than a nlogn implementation in lib FFTW)
	static void discreteCosineTransform(vector<float>& x, vector<float>& outCoefs)
	{
		int M = x.size();
		for(int j = 0; j < M; ++j)
		{
			float Cj = 0;
			for(int i = 0; i < M; ++i)
				Cj += x[i] * cos((float)(j+1) * ((float)(i+1) - .5) * M_PI / (float)M);
			
			outCoefs.push_back(Cj);
		}
	}
};