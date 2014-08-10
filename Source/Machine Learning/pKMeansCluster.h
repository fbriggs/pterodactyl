// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// this is a basic implementation of lloyd's hueristic for the kmeans problem.
// 3 methods of seeding are provided: 
// 1 - uniform random in the range (0, 1)
// 2 - seed by selecting from the examples uniformly at random
// 3 - use the kmeans++ seeding algorithm 

#pragma once

#include "pRandom.h"
#include <vector>

using namespace std;

class pKMeansCluster
{
public:
	int d_, k_, n_;
	vector< vector<float> > centers_;	// centers_[i][j] is the j-th coordinate of the i-th centroid 
	vector<int> labels_;				// exampleLabels_[i] stores the label (i.e. which centroid) example i belongs to
	
	// setup a k-means clusterer for a d-dimensional space, and n examples to cluster
	pKMeansCluster(int d, int k, int n) : 
		d_(d), 
		k_(k), 
		n_(n), 
		labels_(n)						// initialize exampleLables to be n-dimensional
	{
		// initialize clusters
		for(int i = 0; i < k; ++i)
		{
			vector<float> center(d);
			centers_.push_back(center);
		}		
	}
	
	// initialize the centroids to random locations.
	// for now, they are uniformly distributed int the range [0, 1]
	void randomizeCentriods()
	{
		for(int i = 0; i < k_; ++i)
			for(int j = 0; j < d_; ++j)
				centers_[i][j] = pRandom::randInRange(0, 1);
	}
	
	// set the centroids to random examples... maybe this is a better
	// starting condition than uniform?
	void randomizeCentroidsFromExamples(vector< vector<float> >& vs)
	{
		for(int i = 0; i < k_; ++i)
		{
			int example = rand()%n_;
			for(int j = 0; j < d_; ++j)
				centers_[i][j] = vs[example][j];
		}
	}
	
	// set the ith center to the vector v
	void setCenterToVector(int i, vector<float>& v)
	{
		centers_[i].clear();
		for(int j = 0; j < v.size(); ++j)
			centers_[i].push_back(v[j]);
	}
	
	// compute the distance squared from the vector v to center i
	float distanceToCenter(int i, vector<float>& v)
	{
		float sum = 0;
		float diff;
		for(int j = 0; j < v.size(); ++j)
		{
			diff = centers_[i][j] -v[j];
			sum += diff * diff;
		}
		
		return sum;
	}
	
	// this is an implementation of the seeding algorithm described in "k-means++: the advantages of careful seeding", Arthur and Vasssilvitskii
	void randomizeKMeansPlusPlus(vector< vector<float> >& vs)
	{
		cout << "randomizing with kmeans++ seeding..." << endl;
		// pick the first center at random from amongst the example points
		setCenterToVector(0, vs[rand()%vs.size()]);
		
		// pick the remianing centers from the examples according to the 
		// probability distribution given by Arthur and Vasssilvitskii
		for(int i = 1; i < k_; ++i) // pick a new center c_i
		{
			cout << "determining center " << i << endl;
			// at this point, centers 0 through i-1 have been chosen already. 
			// now we must compute the probablility to choose each example as the new center
			// as the distance from that example to the closest existing center, normalized by 
			// the sum of those distances
			vector<float> probabilityToChooseExample;
			for(int j = 0; j < vs.size(); ++j) // for each example, find its probability
			{
				float minDist = 10E100; // basically infinity
				for(int k = 0; k < i; ++k) // for each center that has already been chosen, find the distance from the example to that center
					minDist = min(minDist, distanceToCenter(k, vs[j]));
					
				probabilityToChooseExample.push_back(minDist);
			}
			
			int numExamples = vs.size();
			
			// normalize the probabilities so that they collectively for a pdf
			// NOTE: this normalization could be making a lot of numerical error since the number of things 
			// in the distribution is very large. could be fixed by changing the code to sample from
			// the multinomial distribution so that it doesn't require a normalized pdf
			float sumP = 0;
			for(int j = 0; j < numExamples; ++j) sumP += probabilityToChooseExample[j];
			assert(sumP != 0);
			for(int j = 0; j < numExamples; ++j) probabilityToChooseExample[j] /= sumP;
		
			// choose an example to use as the next center from the multinomial distribution
			// described by probabilityToChooseExample
			int newExampleIndex = pRandom::randomMultinomial(probabilityToChooseExample);
			
			setCenterToVector(i, vs[newExampleIndex]);

		}
	}
	
	// returns the index of the centroid that v is closest to
	int getLabelForExample(vector<float>& v)
	{
		float dontCare;
		return getLabelForExampleAndReturnDistance(v, dontCare);
	}
	
	// returns the index of the centroid that v is closest to, and write that
	// distance to outDist
	int getLabelForExampleAndReturnDistance(vector<float>& v, float& outDist)
	{
		int closestK = 0;
		float closestDist = 10E100; // basically infinity
		
		// for each centroid
		for(int i = 0; i < k_; ++i)
		{
			// compute the euclidean distance squared to that centroid
			// TODO: this is now duplicated code.. refactor.
			float dist = 0;
			for(int j = 0; j < d_; ++j)
			{
				float delta = v[j] - centers_[i][j];
				dist += delta * delta;
			}
			
			
			// if it is less than the current best distance, update the best distance
			if( dist < closestDist )
			{
				closestK = i;
				closestDist = dist;
				
			}
		}
		outDist = closestDist;
		return closestK;
	}
	
	// for each vector in vs, label it with the index of its closest centroid
	// returns the average distance from each example to its closest center
	float mapExamplesToCentroids(vector< vector<float> >& vs)
	{
		int n = vs.size();
		//cout << "mapping examples to centroids, # examples: " << n << endl;
		
		float avgDistToClosestCenter = 0;
		for(int i = 0; i < n; ++i)
		{
			float distToClosestCenter;
			labels_[i] = getLabelForExampleAndReturnDistance(vs[i], distToClosestCenter);
			avgDistToClosestCenter += distToClosestCenter / (float) n;
			//if( i % 100 == 0 )
			//{
			//	cout << "*";
			//	cout.flush();
			//}
		}
		//cout << endl;
		//cout << "average distance to closest center: " << avgDistToClosestCenter << endl;
		return avgDistToClosestCenter;
	}

	// move the centroids to the average position of each vector that is labeled
	// as being in that cluster
	void updateCentroids(vector< vector<float> >& vs)
	{
		
		//cout << "updating centroids" << endl;
		
		// for each centroid
		for(int i = 0; i < k_; ++i)
		{
			vector<float> newCenter(d_); // initially filled with 0s
			int examplesInThisCluster = 0;
			
			// for each example that is in this centroid
			for(int j = 0; j < n_; ++j)
			{
				if( labels_[j] != i ) continue; // only consider examples in this centroid
				++examplesInThisCluster;
				
				for(int x = 0; x < d_; ++x)
					newCenter[x] += vs[j][x];
			}
			
			if( examplesInThisCluster != 0 )
			for(int j = 0; j < d_; ++j)
				newCenter[j] /= (float)examplesInThisCluster;
			
			setCenterToVector(i, newCenter); // this wasn't here before... which means the error never decreased!
			
			//if( i % 100 == 0 )
			//{
			//	cout << "*";
			//	cout.flush();
			//}
		}
		//cout << endl;
	}
	
	// run llyod's heuristic until either it has been performed maxIterations times,
	// or convergence is reached
	void runToConvergenceOrMaxIterations(vector< vector<float> >& fvs, int maxIterations)
	{
		float prevAvgDist = -1;
		for(int i = 0; i < maxIterations; ++i)
		{
			cout << "*"; cout.flush();
			
			//cout << "i = " << i << endl;
			float avgDist = mapExamplesToCentroids(fvs);
			if( avgDist == prevAvgDist)
			{
				cout << "reached convergence" << endl;
				break;
			}
			prevAvgDist = avgDist;
			updateCentroids(fvs);
		}
	}
};

