// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once

#include <iostream>
#include <cstdlib>
#include <math.h>
#include <vector>

using namespace std;

// a collection of functions for generating random numbers 
class pRandom
{
public:
	// returns a random float between the range min and max
	static float randInRange(float min, float max)
	{
		return min + (max - min) * (float) rand() / (float) RAND_MAX;
	}
	
	// returns a random float from a normal distribution with mean zero
	// and standard deviation sigma (not tested for accuracy, but it seems pretty close).
	static float randGaussian(float sigma)
	{
		float x, y, r2;

		do
		{
		  // choose x, y in uniform square (-1,-1) to (+1,+1) 
		  x = -1 + 2 * randInRange(0, 1);
		  y = -1 + 2 * randInRange(0, 1);

		  // see if it is in the unit circle
		  r2 = x * x + y * y;
		}
		while (r2 > 1.0 || r2 == 0);

		// Box-Muller transform
		return sigma * y * sqrt (-2.0 * log (r2) / r2);
	}
	
	// return a random integer between 0 and pdf.size() - 1, according
	// to the multinomial distribution described by pdf
	static int randomMultinomial(vector<float>& pdf)
	{
		float r = randInRange(0, 1);
		float sumP = 0;
		for(int i = 0; i < pdf.size(); ++i)
		{
			sumP += pdf[i];
			if( r < sumP ) 
				return i;
		}
		
		return pdf.size() - 1;
	}

	// shuffle the vector of ints v to a new random permutation.
	// i think this is called djikstra's shuffling algorithm.
	static void shuffle(vector<int>& v)
	{
		// loop invariant: the v[0] through v[i] are shuffled
		for(int i = 0; i < v.size(); ++i)
		{
			// chose a random element that is not shuffled
			int r = i + (rand() % (v.size() - i));
			
			// swap the ith element with the chosen element
			int temp = v[i];
			v[i] = v[r];
			v[r] = temp;
		}
	}
	
	// shuffle a vector<T>
	template <class T>
	static void templatedShuffle(vector<T>& v)
	{
		// loop invariant: the v[0] through v[i] are shuffled
		for(int i = 0; i < v.size(); ++i)
		{
			// chose a random element that is not shuffled
			int r = i + (rand() % (v.size() - i));
			
			// swap the ith element with the chosen element
			T temp = v[i];
			v[i] = v[r];
			v[r] = temp;
		}
	}
	
	// REFACTOR: shuffle for ints should be a template specialization (or might not need to be separately defined at all)
	
};
