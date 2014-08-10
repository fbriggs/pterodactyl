// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// this file defines functions for computing "distances" between bags

#pragma once

#include "pDistanceMetric.h"
#include "pMIMLExample.h"

class pBagDistance
{
public:
	
	// return the average Hausdorff distance between bags A and B
	static float DHavg(pMIMLExample& A, pMIMLExample&B)
	{
		assert(A.instances_.size() != 0 && B.instances_.size() != 0);
		
		float sum = 0;
		
		// for each a in A
		for(int a = 0; a < A.instances_.size(); ++a)
		{
			// find the distance to its nearest neighbor in B
			float minD = FLT_MAX;
			for(int b = 0; b < B.instances_.size(); ++b)
				minD = min(minD, (float)sqrt(distL2Squared(A.instances_[a], B.instances_[b])));
				//minD = min(minD, (float)distL1(A.instances_[a], B.instances_[b]));
			sum += minD;
		}
		
		// for each b in B
		for(int b = 0; b < B.instances_.size(); ++b)
		{
			// find the distance to its nearest neighbor in A
			float minD = FLT_MAX;
			for(int a = 0; a < A.instances_.size(); ++a)
				minD = min(minD, (float)sqrt(distL2Squared(A.instances_[a], B.instances_[b])));
				//minD = min(minD, (float)distL1(A.instances_[a], B.instances_[b]));
			sum += minD;
		}
		
		return sum / ((float)A.instances_.size() + (float)B.instances_.size());
	}
	
	// compute the distance between bags A and B as the distance between
	// the farthest pair (a,b) of instance a in A and b in B
	static float farthestInstancePair(pMIMLExample& A, pMIMLExample&B)
	{
		float maxDist = 0;
		for(int a = 0; a < A.instances_.size(); ++a)
			for(int b = 0; b < B.instances_.size(); ++b)
				maxDist = max(maxDist, sqrt(distL2Squared(A.instances_[a], B.instances_[b])));
		return maxDist;
	}
	
	static float closestInstancePair(pMIMLExample& A, pMIMLExample&B)
	{
		float minDist = FLT_MAX;
		for(int a = 0; a < A.instances_.size(); ++a)
			for(int b = 0; b < B.instances_.size(); ++b)
				minDist = min(minDist, sqrt(distL2Squared(A.instances_[a], B.instances_[b])));
		return minDist;
	}
};