// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#pragma once

#include <string>
#include <vector>
#include <math.h>
using namespace std;

class pSignalUtils
{
public:
	// returns the index of the max value in v
	static int findMax(vector<float>& v)
	{
		float maxVal = -10E50;
		int maxIndex = 0;
		for(int i = 0; i < v.size(); ++i)
			if( v[i] > maxVal )
			{
				maxVal = v[i];
				maxIndex = i;
			}
		
		return maxIndex;
	}
	
	// zero all elements in v in the range a...b inclusive
	static void zeroRange(vector<float>& v, int a, int b)
	{
		for(int i = a; i <= b; ++i) v[i] = 0;
	}
	
	
	// extract a frame of frameSize, starting at sample startOffset from samples. store the output in outFrameSamples
	static void extractFrame(vector<float>& samples, vector<float>& outFrameSamples, int startOffset, int frameSize)
	{
		for(int i = 0; i < frameSize; ++i)
		{
			if( i + startOffset >= samples.size() ) return;
			outFrameSamples.push_back(samples[i + startOffset]);
		}
	}
		
	// do a uniform blur on samples and write the results into outSamples.
	// blur will go a distance of neighborhoodSize away from center (not including center itself). 
	static void uniformBlur(vector<float>& samples, vector<float>& outSamples, int neighborhoodSize)
	{
		for(int i =  0; i < samples.size(); ++i)
			outSamples.push_back(samples[i]);
		
		for(int i = neighborhoodSize; i < samples.size() - neighborhoodSize; ++i)
		{
			for(int j = 1; j <= neighborhoodSize; ++j)
			{
				outSamples[i] += samples[i + j];
				outSamples[i] += samples[i - j];
			}
			outSamples[i] /= (float)(neighborhoodSize * 2 + 1);
		}
	}
	
	// normalize v by finding the max, then dividing all of them by that amount
	// (assumes all positive samples). results in values ranging from 0 to 1.
	static void normalizeVectorToMaxValueOf1(vector<float>& v)
	{
		float maxVal = 0;
		for(int i = 0; i < v.size(); ++i)
			maxVal = max(maxVal, v[i]);
		
		if(maxVal != 0)
			for(int i = 0; i < v.size(); ++i)
				v[i] /= maxVal;
	}	
	
	// normalize the vector v so that its elements sum to 1 (i.e. it is a probability density function)
	static void normalizeVectorToPDF(vector<float>& v)
	{
		float sum = 0;
		for(int i = 0; i < v.size(); ++i)
			sum += v[i];
		
		if( sum == 0 ) return;
		
		for(int i = 0; i < v.size(); ++i)
			v[i] /= sum;
	}
	
	// return the frame corresponding to time t (in seconds)
	static int timeToFrameIndex(float t, int windowSize, int windowStep, int sampleRate)
	{
		return (windowSize / windowStep) * t * (float)sampleRate / (float)windowSize;
	}
	
	// return the time corresponding to frame i
	static float frameIndexToTime(int i, int windowSize, int windowStep, int sampleRate)
	{
		return (float)i / (((float)windowSize / (float)windowStep) * (float)sampleRate / (float)windowSize);
	}
	
	// compute the Kullback-Leibler divergence of P from Q
	static float computeKLDivergence(vector<float>& P, vector<float>& Q)
	{
		float sum = 0;
		for(int i = 0; i < P.size(); ++i)
			sum += P[i] * log(P[i] / Q[i]);
		return sum;
	}
	
	static float computeDistSquared(vector<float>& X, vector<float>& Y)
	{
		float sum = 0;
		for(int i = 0; i < X.size(); ++i)
		{
			float diff = X[i] - Y[i];
			sum += diff * diff;
		}
		return sum;
	}
	
	// make a vector that is a uniform pdf with k elements
	static void makeUniformDistribution(int k, vector<float>& outU)
	{
		for(int i = 0; i < k; ++i) outU.push_back(1.0f / (float)k);
	}
	
	static float coefToFrequency(float coef, float bandPassLow, float samplingFrequency, float windowSize)
	{
		// since we are dropping the first 8 coefs, we need to add 8 to get back to a normal coef index
		float c = coef + bandPassLow + 1; // the + 1 here is to go from indices starting at 0 to indices starting at 1... not sure about this
		float freq = samplingFrequency * (c-1)/windowSize;
		return freq;
	}
};