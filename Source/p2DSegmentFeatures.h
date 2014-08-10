// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// this file defines functions for computing a feature vector to describe a (2D) segment of audio from the spectrogram and segmentation mask
#pragma once

#include "pImageBytes.h"
#include "pImageUtils.h"
#include "pVectorND.h"
#include "pFastGaussianBlur.h"
#include "pVectorUtils.h"
#include "pRectangle.h"

class p2DSegmentFeatures
{
public:
	/// histogram of gradients feature ///
	
	static vector<pVectorND> gradientCodebook_;
	const static int numGradientCodebookElements = 16;
	
	static void initGradientCodebook()
	{
		for(int i = 0; i < numGradientCodebookElements; ++i)
		{
			float theta = ((float)i / (float)numGradientCodebookElements) * 2.0 * M_PI;
			gradientCodebook_.push_back(pVectorND(cos(theta), sin(theta)));
		}
	}
	
	static int findGradientBin(pVectorND& g)
	{
		vector<float> dots;
		for(int i = 0; i < numGradientCodebookElements; ++i)
		{
			dots.push_back(g.dot(gradientCodebook_[i]));
		}
		return pVectorUtils::argmax(dots);
	}
	
	static vector<float> hogFeatures(pImageBytes& croppedSpec, pImageBytes& croppedMask)
	{
		vector<float> histogram(numGradientCodebookElements);
		
		pImageBytes blurred = pFastGaussianBlur::blur(croppedSpec, 3);

		float minMag2 = 0.01;
		
		for(int x = 0; x < croppedSpec.w_; ++x)
		for(int y = 0; y < croppedSpec.h_; ++y)
		if( croppedMask.getPixel(x, y, 0) > 0)
		{
			pVectorND g = pImageUtils::grad(blurred, x, y);
			//pVectorND g = grad(croppedSpec, x, y);
			if( g.lengthSquared() < minMag2 ) continue;
			histogram[findGradientBin(g)] += 1;
		}
		pVectorUtils::normalizeVectorToPDF(histogram);
		
		return histogram;
	}
		
	/// time/frequency profile descriptors & other statistics ///
	
	// take the cropped spectrogram of a segmentand compute properties
	// of its time and frequency profiles, and other statistical properties
	static vector<float> profileFeatures(pImageBytes& croppedSpec, pImageBytes& croppedMask)
	{
		vector<float> freq_profile = frequencyProfile(croppedSpec);		// can be used as a feature by iteslf. equivelent to avg spectrum pdf
		vector<float> time_profile = timeProfile(croppedSpec);			// this cant be used directly because it varies in length
		
		float freq_gini					= pVectorUtils::giniForDistribution(freq_profile);	// "uncertainty" of the profiles
		float time_gini					= pVectorUtils::giniForDistribution(time_profile);
		
		float rel_freq_mean				= relative_mean(freq_profile);
		float rel_time_mean				= relative_mean(time_profile);
		
		float rel_freq_variance			= relative_kth_centeral_moment(freq_profile, rel_freq_mean, 2);
		float rel_time_variance			= relative_kth_centeral_moment(time_profile, rel_time_mean, 2);
		
		float rel_freq_skewness			= relative_kth_centeral_moment(freq_profile, rel_freq_mean, 3);
		float rel_time_skewness			= relative_kth_centeral_moment(time_profile, rel_time_mean, 3);
		
		float rel_freq_kurtosis			= relative_kth_centeral_moment(freq_profile, rel_freq_mean, 4);
		float rel_time_kurtosis			= relative_kth_centeral_moment(time_profile, rel_time_mean, 4);
		
		float rel_freq_maxima			= (float)pVectorUtils::argmax(freq_profile) / (freq_profile.size() - 1);
		float rel_time_maxima			= (float)pVectorUtils::argmax(time_profile) / (time_profile.size() - 1);
		
		float avg_intensity_in_mask		= avgInMaskedRegion(croppedSpec, croppedMask);
		float stdev_intensity_in_mask	= stdDevInMaskedRegion(croppedSpec, croppedMask, avg_intensity_in_mask);
		
		vector<float> features;
		
		features.push_back(freq_gini);
		features.push_back(time_gini);
		features.push_back(rel_freq_mean);
		features.push_back(rel_time_mean);
		features.push_back(rel_freq_variance);
		features.push_back(rel_time_variance);
		features.push_back(rel_freq_skewness);
		features.push_back(rel_time_skewness);
		features.push_back(rel_freq_kurtosis);
		features.push_back(rel_time_kurtosis);
		features.push_back(rel_freq_maxima);
		features.push_back(rel_time_maxima);
		features.push_back(avg_intensity_in_mask);
		features.push_back(stdev_intensity_in_mask);
		
		return features;
	}
	
	// compute the mean of f, but rescale so the range of x-values is [0,1]
	static float relative_mean(vector<float>& f)
	{
		float E = 0;
		for(int i = 0; i < f.size(); ++i)
		{
			float x_i = (float)i / ((float)f.size() - 1);
			E += x_i * f[i];
		}
		return E;
	}
	
	// compute the kth releative moment of f, but rescale so the range of x-values is [0,1]
	// mu_k = E[(x-mu)^k]
	static float relative_kth_centeral_moment(vector<float>& f, float mean, float k)
	{
		float mu_k = 0;
		for(int i = 0; i < f.size(); ++i)
		{
			float x_i = (float)i / ((float)f.size() - 1);
			mu_k += pow(x_i - mean, k) * f[i];
		}
		return mu_k;
	}
	
	static float avgInMaskedRegion(pImageBytes& croppedSpec, pImageBytes& croppedMask)
	{
		float sum = 0;
		int num = 0;
		for(int x = 0; x < croppedSpec.w_; ++x)
		for(int y = 0; y < croppedSpec.h_; ++y)
		if( croppedMask.getPixel(x, y, 0) > 0 )
		{
			sum += croppedSpec.getPixel(x, y, 0);
			++num;
		}
		return sum / (float)num;
	}
	
	static float stdDevInMaskedRegion(pImageBytes& croppedSpec, pImageBytes& croppedMask, float avg)
	{
		float sum = 0;
		int num = 0;
		for(int x = 0; x < croppedSpec.w_; ++x)
		for(int y = 0; y < croppedSpec.h_; ++y)
		if( croppedMask.getPixel(x, y, 0) > 0 )
		{
			sum += (croppedSpec.getPixel(x, y, 0) - avg) * (croppedSpec.getPixel(x, y, 0) - avg);
			++num;
		}
		return sqrt(sum / (float)num);
	}
	
	
	// takes a spectrogram as input and computes the frequency 
	// profile (i.e. sum over rows/time), normalized to a pdf
	static vector<float> frequencyProfile(pImageBytes& src)
	{
		vector<float> profile;
		
		for(int y = 0; y < src.h_; ++y)
		{
			float sum = 0;
			for(int x = 0; x < src.w_; ++x)
				sum += src.getPixel(x, y, 0);
			profile.push_back(sum);
		}
		pVectorUtils::normalizeVectorToPDF(profile);
		return profile;
	}
	
	// takes a spectrogram as input and computes the time 
	// profile (i.e. sum over cols/freq), normalized to a pdf
	static vector<float> timeProfile(pImageBytes& src)
	{
		vector<float> profile;
		
		for(int x = 0; x < src.w_; ++x)		
		{
			float sum = 0;			
			for(int y = 0; y < src.h_; ++y)
				sum += src.getPixel(x, y, 0);
			profile.push_back(sum);
		}
		pVectorUtils::normalizeVectorToPDF(profile);
		return profile;
	}
	
	/// operations on binary masks ///
	
	// compute a vector of features describing the shape of the component mask
	static vector<float> maskShapeFeatures(pImageBytes& src, pRectangle& rect)
	{
		vector<float> features;
		
		int& minX = rect.minX_;
		int& maxX = rect.maxX_;
		int& minY = rect.minY_;
		int& maxY = rect.maxY_;
				
		float min_freq			= minY;									features.push_back(min_freq);
		float max_freq			= maxY;									features.push_back(max_freq);
		float bandwidth			= maxY - minY;							features.push_back(bandwidth);
		float duration			= maxX - minX;							features.push_back(duration);

		float area				= pImageUtils::maskArea(src);			features.push_back(area);
		float perimeter			= pImageUtils::maskPerimeter(src);		features.push_back(perimeter);
		
		float non_compactness	= perimeter * perimeter / area;			features.push_back(non_compactness);
		float rectangularity	= area / (bandwidth * duration);		features.push_back(rectangularity);
		
		return features;
	}
};