// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// this file implements a farthest first traversal of a dataset. 
// farthest-first is useful for identifying all of the examples of each class in a datase while inspecting as few examples as possible (the class-discovery problem, e.g. species discovery).

#pragma once

#include <hash_set.h>
#include <vector>
#include <algorithm>
#include "pVectorUtils.h"

using namespace std;

class pFarthestFirstTraversal
{
public:
	// fvs contains an array of feature vectors.  starting with startIndex, 
	// visit numToTraverse points, and return their indices. the order is defined greedily by
	// always selecting the next point that is farthest from all currently selected points.
	static vector<int> traverseMaxMin(vector< vector<float> >& fvs, int numToTraverse)
	{
		// start at an arbitrary point
		int startIndex = 0;
		
		vector<int> traversal;
		hash_set<int> inSet;
		traversal.push_back(startIndex);
		inSet.insert(startIndex);
		
		for(int selectingExample = 1; selectingExample < numToTraverse; ++selectingExample) // select the i'th example
		{
			cout << "selecting example " << selectingExample << endl;
			
			vector<float> candidateDistances(fvs.size(), -1); // the default value is -1. this will never be the max in a collection of valid distances
		
			// consider example j as a candidate to be added to the set
			for(int j = 0; j < fvs.size(); ++j) // no room for optimization here: because k << n, there isn't much gain to be had here by somehow processing (n-k) examples instead of n
			{
				if( j % 1000 == 0 ) cout << "*"; cout.flush();
				
				if( inSet.count(j) != 0 ) continue; // if example j is already in the set, dont try to add it again
				
				// compute the distance from example j to all examples currently in the set, and 
				// take the min over all of those distances
				float minDistToSet = FLT_MAX;
				for(int k = 0; k < traversal.size(); ++k)
					minDistToSet = min(minDistToSet, pVectorUtils::distanceSquared(fvs[j], fvs[traversal[k]]));

				candidateDistances[j] = minDistToSet;
			}
			
			int exampleToSelect = pVectorUtils::argmax(candidateDistances);
			
			cout << "dist to selected example = " << candidateDistances[exampleToSelect] << endl;			
			traversal.push_back(exampleToSelect);
			inSet.insert(exampleToSelect);
		}
		
		return traversal;
	}
	
	// alternative to max-min traversal, the max-average traversal selects the point with the greatest average distance to points in the current set
	static vector<int> traverseMaxAvgDist(vector< vector<float> >& fvs, int numToTraverse)
	{
		// start at an arbitrary point
		int startIndex = 0;
		
		vector<int> traversal;
		hash_set<int> inSet;
		traversal.push_back(startIndex);
		inSet.insert(startIndex);
		
		for(int selectingExample = 1; selectingExample < numToTraverse; ++selectingExample) // select the i'th example
		{
			cout << "selecting example " << selectingExample << endl;
			
			vector<float> candidateDistances(fvs.size(), -1); // the default value is -1. this will never be the max in a collection of valid distances
		
			// consider example j as a candidate to be added to the set
			for(int j = 0; j < fvs.size(); ++j) // no room for optimization here: because k << n, there isn't much gain to be had here by somehow processing (n-k) examples instead of n
			{
				if( j % 1000 == 0 ) cout << "*"; cout.flush();
				
				if( inSet.count(j) != 0 ) continue; // if example j is already in the set, dont try to add it again
				
				// compute the average distance to the set from example j
				float avgDist = 0;
				for(int k = 0; k < traversal.size(); ++k)
					avgDist += sqrt(pVectorUtils::distanceSquared(fvs[j], fvs[traversal[k]]));
				avgDist /= float(traversal.size());	
				
				candidateDistances[j] = avgDist;
			}
			
			int exampleToSelect = pVectorUtils::argmax(candidateDistances);
			
			cout << "dist to selected example = " << candidateDistances[exampleToSelect] << endl;			
			traversal.push_back(exampleToSelect);
			inSet.insert(exampleToSelect);
		}
		
		return traversal;
	}
};