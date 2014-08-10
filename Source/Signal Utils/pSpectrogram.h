// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// this class computes and stores a spectrogram from wav data

#pragma once

#include "pWavData.h"
#include "FFT.h"
#include "pImageBytes.h"
#include "pMFCCs.h"
#include "pSignalUtils.h"
#include <iostream>
#include <math.h>
#include <vector>
using namespace std;


class pSpectrogram
{
public:
	int windowSize_;
	int spectrumBins_;
	int numFrames_;
	int bandPassLow_;
	int bandPassHigh_;
	int samplingFrequency_;
	
	vector< vector<float> > spectra_;
	vector< vector<float> > frameMFCs_;		// mel frequency coefs
	vector< vector<float> > frameMFCCs_;	// mel frequency cepstral coefs
	int numMFCsPerFrame_;
	int numMFCCsPerFrame_;
	
	~pSpectrogram() {spectra_.clear(); } // this should happen automatically, but maybe it isn't... 
	
	// make a spectrogram from wavData
	pSpectrogram(pWavData& wavData, int windowSize, int windowStep, int bandPassLow, int bandPassHigh)
	{
		bandPassLow_ = bandPassLow;
		samplingFrequency_ = wavData.sampleRate_;
		windowSize_ = windowSize;
		
		int specBinsNoBandPass = windowSize / 2;
		
		spectrumBins_ = specBinsNoBandPass - bandPassLow - bandPassHigh;
		
		vector<float>& samples = wavData.samples_;
		float* windowSamples = new float[windowSize];
		float* realCoefs = new float[specBinsNoBandPass];
		float* imagCoefs = new float[specBinsNoBandPass];
		
		int currOffset = 0;
		while(currOffset + windowSize < samples.size())
		{
			// extract the samples in this frame
			for(int i = 0; i < windowSize; ++i)
				windowSamples[i] = samples[i + currOffset];
			
			// apply a hamming window to frame samples
			WindowFunc(2, windowSize, windowSamples); 
			
			// do an FFT from real samples to complex FFT coefs
			RealFFT(windowSize, windowSamples, realCoefs, imagCoefs);
			
			// save the spectrum as the absolute magnitude of the complex fft coefficients
			vector<float> spectrum;
			for(int i = bandPassLow; i < specBinsNoBandPass - bandPassHigh; ++i)
			{
				float coef = sqrt(realCoefs[i] * realCoefs[i] + imagCoefs[i] * imagCoefs[i]);
				assert( !isnan(coef) );
				spectrum.push_back( coef );
			}
			
			spectra_.push_back(spectrum);
			
			currOffset += windowStep;
		}
		
		delete[] windowSamples;
		delete[] realCoefs;
		delete[] imagCoefs;
		
		numFrames_ = spectra_.size();
	}
	
	// set the first n elements (low frequency) in the spectrogram to 0.
	// this is equivelent to a high pass filter. we will keep these
	// elements around even though they are 0 to make the MFCC calculations simpler
	void zeroLowFreqElements(int numElements)
	{
		for(int i = 0; i < spectra_.size(); ++i)
			for(int j = 0; j < numElements; ++j)
				spectra_[i][j] = 0;
	}
	
	// make the spectra for each frame into a probability density function (i.e. elements add to 1). 
	void makeSpectraIntoPDFs()
	{
		// for each frame
		for(int i = 0; i < spectra_.size(); ++i)
		{
			float sum = 0;
			for(int j = 0; j < spectrumBins_; ++j)
				sum += spectra_[i][j];
			
			if( sum != 0 )
				for(int j = 0; j < spectrumBins_; ++j)
					spectra_[i][j] /= sum;
		}
	}
	
	void normalizeSpectraTo01()
	{
		float m = 0;
		for(int i = 0; i < spectra_.size(); ++i)
			for(int j = 0; j < spectrumBins_; ++j)
				m = max(m,spectra_[i][j]);
		
		
		assert(m!=0);
		
		// for each frame
		if( m != 0 )
			for(int i = 0; i < spectra_.size(); ++i)
				for(int j = 0; j < spectrumBins_; ++j)
				{
					spectra_[i][j] /= m;
					assert(!isnan(spectra_[i][j]));
					assert(!isinf(spectra_[i][j]));
				}
	}
	
	// this is a basic de-noising algorithm. we compute the average spectrum over all frames,
	// then subtract it from each frame (and if a result is < 0, we set it to 0).
	void subtractAvgSpectrum()
	{
		// compute the average spectrum
		vector<float> avgSpectrum(spectrumBins_);
		for(int i = 0; i < spectra_.size(); ++i)
			for(int j = 0; j < spectrumBins_; ++j)
				avgSpectrum[j] += spectra_[i][j];
		
		for(int j = 0; j < spectrumBins_; ++j)
			avgSpectrum[j] /= (float)spectra_.size();
		
		// subtract it from each frame
		for(int i = 0; i < spectra_.size(); ++i)
			for(int j = 0; j < spectrumBins_; ++j)
				spectra_[i][j] = max(0.0f, spectra_[i][j] - avgSpectrum[j]);
	}
	
	void computeSmoothedEnergyEnvelope(vector<float>& smoothedEnergy, int blurSize)
	{
		vector<float> energyForFrame;
		for(int i = 0; i < spectra_.size(); ++i)
		{
			float energy = 0;
			for(int j = 0; j < spectrumBins_; ++j) energy += spectra_[i][j];
			energyForFrame.push_back(energy);
		}
		
		pSignalUtils::uniformBlur(energyForFrame, smoothedEnergy, blurSize);
		pSignalUtils::normalizeVectorToMaxValueOf1(smoothedEnergy);
		
	}
	
	void boostContrast()
	{
		for(int i = 0; i < spectra_.size(); ++i)
			for(int j = 0; j < spectrumBins_; ++j)
				spectra_[i][j] = sqrt(spectra_[i][j]);
	}
	
	// convert a fourier coeffient index to its corresponding frequency
	float coefToFrequency(int coef)
	{
		// since we are dropping the first 8 coefs, we need to add 8 to get back to a normal coef index
		float c = coef + bandPassLow_ + 1; // the + 1 here is to go from indices starting at 0 to indices starting at 1... not sure about this
		float freq = samplingFrequency_ * (c-1)/windowSize_;
		return freq;
	}
	
	// make an image from this spectrogram. assumes that makeSpectraIntoPDFs has been called already
	pImageBytes makeGreyscaleImage()
	{
		pImageBytes img(spectra_.size(), spectrumBins_, 1);
		
		for(int i = 0; i < spectra_.size(); ++i)
			for(int j = 0; j < spectrumBins_; ++j)
				img.setPixelF(i, j, 0, spectra_[i][j]);
		
		return img;
	}
	
	// similar to makeGreyscaleImage, but boosts the contrast of the image
	pImageBytes makeGreyscaleImageBoostContrast()
	{
		pImageBytes img(spectra_.size(), spectrumBins_, 1);
		
		for(int i = 0; i < spectra_.size(); ++i)
			for(int j = 0; j < spectrumBins_; ++j)
				img.setPixelF(i, j, 0, sqrt(spectra_[i][j]));
		
		return img;
	}
	
	/*
	 // makes  a color spectrogram image by applying a color map
	 pImage* makeColorImage()
	 {
	 pImage* img = new pImage(spectra_.size(), spectrumBins_, 3);
	 
	 for(int i = 0; i < spectra_.size(); ++i)
	 for(int j = 0; j < spectrumBins_; ++j)
	 {
	 Pixel p;
	 float h = spectra_[i][j] * 360 * 4 + 180;
	 float s = 1;
	 float v = .5 +  .5 * spectra_[i][j];
	 pColorSpaces::HSVtoRGB( h, s, v, p.r, p.g, p.b);
	 img->setPixel(i, j, p);
	 }
	 
	 return img;
	 }*/
	
	// don't normalize the spectrum before calling this
	void computeMFCCs()
	{		
		for(int i = 0; i < spectra_.size(); ++i)
		{
			for(int j = 0; j < spectrumBins_; ++j)
			{
				assert(!isnan(spectra_[i][j]));
				assert(!isinf(spectra_[i][j]));
			}
			
			vector<float> mfcs;
			vector<float> mfccs;
			pMFCCs::getFrameMFCCs(spectra_[i], mfccs, mfcs);
			frameMFCCs_.push_back(mfccs);
			frameMFCs_.push_back(mfcs);
		}
		
		numMFCCsPerFrame_ = frameMFCCs_[0].size();
		numMFCsPerFrame_ = frameMFCs_[0].size();
	}
	
	/*
	 pImage* makeMFCCImage()
	 {
	 int numMFCCsPerFrame = frameMFCCs_[0].size();
	 
	 pImage* img = new pImage(spectra_.size(), numMFCCsPerFrame, 1);
	 
	 // find the max mfcc, so we can normalize the image
	 float maxMFCC = 0;
	 float minMFCC = 0;
	 for(int i = 0; i < spectra_.size(); ++i)
	 for(int j = 0; j < numMFCCsPerFrame; ++j)
	 {
	 maxMFCC = max(maxMFCC, frameMFCCs_[i][j]);
	 minMFCC = min(minMFCC, frameMFCCs_[i][j]);
	 }
	 
	 assert( !isnan( minMFCC) );
	 assert( !isnan( maxMFCC) );
	 
	 float scale = maxMFCC - minMFCC;
	 assert( scale != 0);
	 
	 for(int i = 0; i < spectra_.size(); ++i)
	 for(int j = 0; j < numMFCCsPerFrame; ++j)
	 {
	 float mfccForPixel = frameMFCCs_[i][j];
	 
	 assert( !isnan(mfccForPixel));
	 
	 float p = (mfccForPixel - minMFCC) /scale;
	 
	 //cout << p << " " << mfccForPixel << " " << scale << " " << minMFCC << " " << maxMFCC << endl;
	 
	 assert( !isnan(p));
	 assert( !isinf(p));
	 img->setPixel_(i, j, 0, p);
	 }
	 
	 return img;
	 }
	 
	 pImage* makeMFCImage()
	 {
	 int numMFCsPerFrame = frameMFCs_[0].size();
	 
	 pImage* img = new pImage(spectra_.size(), numMFCsPerFrame, 1);
	 
	 // find the max mfcc, so we can normalize the image
	 float maxMFC = 0;
	 float minMFC = 0;
	 for(int i = 0; i < spectra_.size(); ++i)
	 for(int j = 0; j < numMFCsPerFrame; ++j)
	 {
	 maxMFC = max(maxMFC, frameMFCs_[i][j]);
	 minMFC = min(minMFC, frameMFCs_[i][j]);
	 }
	 
	 assert( !isnan( minMFC) );
	 assert( !isnan( maxMFC) );
	 
	 float scale = maxMFC - minMFC;
	 assert( scale != 0);
	 
	 for(int i = 0; i < spectra_.size(); ++i)
	 for(int j = 0; j < numMFCsPerFrame; ++j)
	 {
	 float mfcForPixel = frameMFCs_[i][j];
	 
	 assert( !isnan(mfcForPixel));
	 
	 float p = (mfcForPixel - minMFC) /scale;
	 
	 assert( !isnan(p));
	 assert( !isinf(p));
	 img->setPixel_(i, j, 0, p);
	 }
	 
	 return img;
	 }
	 */
};