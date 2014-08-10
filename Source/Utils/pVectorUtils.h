// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// misc. functions for working with STL vectors...
#pragma once

#include <iostream>
#include <vector>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <assert.h>

using namespace std;

class pVectorUtils
{
public:
	static float sample_mean(vector<float>& x)
	{
		float sum = 0;
		for(int i = 0; i < x.size(); ++i)
			sum += x[i];
		return sum / float(x.size());
	}
	
	static float sample_var(vector<float>& x, float mean)
	{
		float s2 = 0;
		for(int i = 0; i < x.size(); ++i)
			s2 += (x[i] - mean) * (x[i] - mean);
		return s2 / float(x.size() - 1);
	}
	
	template <class T> static bool vector_contains(vector<T>& v, const T& x)
	{
		for(int i = 0; i < v.size(); ++i)
			if( v[i] == x ) return true;
		return false;
	}
	
	// various versions of argmin and argmax
	static int argmax(vector<int>& v)
	{
		int maxSoFar = INT_MIN;
		int maxIndex = 0;
		for(int i = 0; i < v.size(); ++i)
		{
			if( v[i] > maxSoFar )
			{
				maxSoFar = v[i];
				maxIndex = i;
			}
		}
		return maxIndex;
	}
	
	static int argmax(vector<float>& v)
	{
		float maxSoFar = -FLT_MAX;
		int maxIndex = 0;
		for(int i = 0; i < v.size(); ++i)
		{
			if( v[i] > maxSoFar )
			{
				maxSoFar = v[i];
				maxIndex = i;
			}
		}
		return maxIndex;
	}
	
	static int argmin(vector<int>& v)
	{
		int minSoFar = INT_MAX;
		int minIndex = 0;
		for(int i = 0; i < v.size(); ++i)
		{
			if( v[i] < minSoFar )
			{
				minSoFar = v[i];
				minIndex = i;
			}
		}
		return minIndex;
	}
	
	static int argmin(vector<float>& v)
	{
		float minSoFar = FLT_MAX;
		int minIndex = 0;
		for(int i = 0; i < v.size(); ++i)
		{
			if( v[i] < minSoFar )
			{
				minSoFar = v[i];
				minIndex = i;
			}
		}
		return minIndex;
	}
	
	// if v contains x, returns the index of x, otherwise -1
	template <class T> static int index_of(vector<T>& v, const T& x)
	{
		for(int i = 0; i < v.size(); ++i)
			if( v[i] == x ) return i;
		return -1;
	}
	
	// normalize a vector so its elements add up to 1
	static void normalizeVectorToPDF(vector<float>& w)
	{
		float sum = 0;
		for(int i = 0; i < w.size(); ++i) sum += w[i];
		if( sum == 0 ) return;
		
		for(int i = 0; i < w.size(); ++i) w[i] /= sum;
	}
	
	static void normalizeToUnit(vector<float>& w)
	{
		float sum = 0;
		for(int i = 0; i < w.size(); ++i) sum += w[i] * w[i];
		if( sum == 0 ) return;
		
		sum = sqrt(sum);
		
		for(int i = 0; i < w.size(); ++i) w[i] /= sum;
	}
	
	static float giniForDistribution(vector<float>& prob)
	{
		float gini = 1;
		for(int i = 0; i < prob.size(); ++i) gini -= prob[i] * prob[i];
		return gini;
	}
	
	template <class T> static void cout_vector(vector<T>& v, string delim = ",")
	{
		for(int i = 0; i < v.size(); ++i) cout << v[i] << (i == v.size() - 1 ? "" : delim);
		//cout << endl;
	}
	
	template <class T> static void append(vector<T>& src, vector<T>& dst)
	{
		for(int i = 0; i < src.size(); ++i) dst.push_back(src[i]);
	}
	
	// compute the squared distance between two vectors of floats
	// refactor: if more distance metrics end up in the project, they should be moved to pDistanceMetric.h
	static float distanceSquared(const vector<float>& x, const vector<float>& y)
	{
		assert(x.size() == y.size());
		float sum = 0;
		for(int i = 0; i < x.size(); ++i)
			sum += (x[i] - y[i]) * (x[i] - y[i]);
		return sum;
	}
};