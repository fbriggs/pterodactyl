// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// segment a spectrogram by using a random forest to classify each pixel
#pragma once

#include "pExample.h"
#include "pImageBytes.h"
#include "pRandom.h"
#include "pRandomForest_pthread.h"

#include <vector>
using namespace std;

class p2DRFSegmentation
{
public:
	const static int featureBoxRadius = 8; // this is the setting for the full 2009 dataset
		
	// appends results to exs, an existing vector of examples
	static void generateExamplesFromAnnotatedSpectrogram(pImageBytes& annotation, pImageBytes& spectrogram, vector<pExample>& exs)
	{
		 
		// these are the parameters for the full 2009 experiment 
		float bgSampleProb = 0.01; // fraction of background / negative pixels that will be sampled
		float positiveSampleProb = 0.1; // fraction of positive / bird sound pixels that will be sampled
		float negSampleProb = 0.1; // fraction of blue / explicitly negative pixels that will be sampled
		
		// TODO: don't use pixels on the edge as training examples?
		for(int x = 0; x < annotation.w_; ++x)
		for(int y = 0; y < annotation.h_; ++y)
		{
			int label = 0;
			if( annotation.getPixel(x, y, 0) == 255 && annotation.getPixel(x, y, 1) == 0 && annotation.getPixel(x, y, 2) == 0 )
				label = 1;
			if( annotation.getPixel(x, y, 0) == 0 && annotation.getPixel(x, y, 1) == 0 && annotation.getPixel(x, y, 2) == 255 )
				label = 2;
			
			if( (label == 0 && pRandom::randInRange(0, 1) <= bgSampleProb) || 
				(label == 1 && pRandom::randInRange(0,1) <= positiveSampleProb ) ||
				(label == 2 && pRandom::randInRange(0,1) <= negSampleProb ) )
			{
				if( label == 2 ) label = 0;
				vector<float> fv = getFeatureVector(spectrogram, x, y);
				exs.push_back(pExample(fv, label));
			}
		}
	}
	
	// returns a probability map for segmentation given an input spectrogram and pre-trained random forest
	static pImageBytes segmentationProbabilityMap(pImageBytes& spectrogram, pRandomForest_pthread& rf)
	{
		pImageBytes probMap(spectrogram.w_, spectrogram.h_, 1);
		
		for(int x = 0; x < spectrogram.w_; ++x)
		for(int y = 0; y < spectrogram.h_; ++y)
		{
			vector<float> fv = getFeatureVector(spectrogram, x, y);
			float p = rf.estimateClassProbabilities(fv)[1];
			probMap.setPixelF(x, y, 0, p);
		}
		
		return probMap;
	}
	
	static vector<float> getFeatureVector(pImageBytes& spectrogram, int x, int y)
	{
		vector<float> fv;
		
		float sum = 0;
		int numPixelsInWindow = 0;
		
		for(int u = x - featureBoxRadius; u <= x + featureBoxRadius; ++u)
		for(int v = y - featureBoxRadius; v <= y + featureBoxRadius; ++v)
		{
			float p = spectrogram.getPixelF_(u, v, 0); // uses protected acces to return 0 for pixels off the edge
			fv.push_back(p);
			
			sum += p;
			
			++numPixelsInWindow;
		}
		
		fv.push_back(sum/float(numPixelsInWindow)); // average intensity in the window
		fv.push_back(y); // the center "frequency" of the window (in units of bins) 
		
		
		// this feature was originally added for the ICML 2013 bird challenge to help with segmentation in rain
		// it seems to work pretty well, so it is being included from now on. to go one step further, this now includes the frame before and after as well
		for(int j = 0; j < spectrogram.h_; ++j)
		{
			fv.push_back(spectrogram.getPixelF_(x, j));
			//fv.push_back(spectrogram.getPixelF_(x-1, j)); // idea: add frame before/after as well... may dilute the feature
			//fv.push_back(spectrogram.getPixelF_(x+1, j));
		}
		
		
		// TODO: add more features?
		
		return fv;
	}
};