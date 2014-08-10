// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// the MIML class discovery problem is: given a MIML dataset, query bags for their label sets (in a batch of k) to find as many classes as possible

#pragma once

#include "pRandom.h"
#include "pMIMLExample.h"
#include "pFarthestFirstTraversal.h"
#include "pKMeansCluster.h"
#include "pBagDistance.h"
#include <pthread.h>
#include <map>

using namespace std;

// this struct is used to run multiple experiments in parallel with pthreads
struct pMIMLClassDiscoveryExperiment_runExperiment_args
{
	int threadNum_;
	string discoveryMethod_;
	int numToSelect_;
	int numRepetitions_;
	int param1_, param2_; // these are parameters that can be passed to the discovery method
	
	vector<pMIMLExample>* mimlExs_;
	vector< vector<int> >* classDiscoveryCurves_; // results are written here
};

pthread_mutex_t pMIMLClassDiscoveryExperiment_mutex = PTHREAD_MUTEX_INITIALIZER;

class pMIMLClassDiscoveryExperiment
{
public:
	// selects bags randomly (without replacement). returns a vector of selected bag indices 
	static vector<int> discoveryMethod_random(vector<pMIMLExample>& mimlExs, int numToSelect)
	{
		vector<int> allBagIndices;
		for(int i = 0; i < mimlExs.size(); ++i)
			allBagIndices.push_back(i);
		
		pRandom::shuffle(allBagIndices);
		
		vector<int> selectedBags;
		for(int i = 0; i < numToSelect; ++i)
			selectedBags.push_back(allBagIndices[i]);
		
		return selectedBags;
	}
	
	// selects bags randomly (without replacement), but ignores bags with 0 instances.
	// after all non-empty bags have been exhausted, it will then start randomly selecting bags with 0 instances.
	// note: this may miss some bags with non-empty label sets if segmentation false negatives occur.
	static vector<int> discoveryMethod_randomNonEmpty(vector<pMIMLExample>& mimlExs, int numToSelect)
	{
		
		
		vector<int> emptyBags;
		vector<int> nonEmptyBags;
		
		for(int i = 0; i < mimlExs.size(); ++i)
		{
			if( mimlExs[i].instances_.size() == 0 )
				emptyBags.push_back(i);
			else
				nonEmptyBags.push_back(i);
		}
		
		pRandom::shuffle(emptyBags);
		pRandom::shuffle(nonEmptyBags);
		pVectorUtils::append(emptyBags, nonEmptyBags); // append the empty bag list to the end of the non-empty bag list
		
		vector<int> selectedBags;
		for(int i = 0; i < numToSelect; ++i)
			selectedBags.push_back(nonEmptyBags[i]);
		
		return selectedBags;
	}
	
	// this method cheats by using the bag label sets (which aren't available in practice), and uses the
	// greedy algorithm for the max-coverage problem (which is the best possible approximation assuming p != np).
	// hence this can be viewed as an upper-bound on the performance that can be achieved
	static vector<int> discoveryMethod_coverageOracle(vector<pMIMLExample>& mimlExs, int numToSelect)
	{
		vector<int> selectedBags;
		set<int> knownClasses;
		
		for(int i = 0; i < numToSelect; ++i)
		{
			int bestBag = 0;
			int bestNewClasses = 0;
			
			for(int j = 0; j < mimlExs.size(); ++j) // consider bag j as a candidate for selection
			{
				set<int> newKnownClasses(knownClasses);
				
				for(int k = 0; k < mimlExs[j].labels_.size(); ++k)
					newKnownClasses.insert(mimlExs[j].labels_[k]);
				
				int newNumKnownClasses = newKnownClasses.size();
				
				if( newNumKnownClasses > bestNewClasses )
				{
					bestNewClasses = newNumKnownClasses;
					bestBag = j;
				}
			}
			
			selectedBags.push_back(bestBag);
			
			// update the official know classes set with the selected bag
			for(int k = 0; k < mimlExs[bestBag].labels_.size(); ++k)
				knownClasses.insert(mimlExs[bestBag].labels_[k]);
		}
		
		return selectedBags;
	}
	
	// cluster the instances using k-means++, then find the instance closest to each cluster center, and select the corresponding
	// bag. the clusters will be ordered from largest to smallest (in terms of # of instances they contain).
	static vector<int> discoveryMethod_clusterCenters(vector<pMIMLExample>& mimlExs, int numToSelect)
	{
		// put the instances in a datastructure that can be used by k-means
		vector< vector<float> > fvs;
		for(int i = 0; i < mimlExs.size(); ++i)
			for(int j = 0; j < mimlExs[i].instances_.size(); ++j)
				fvs.push_back(mimlExs[i].instances_[j]);
		
		// run kmeans++ on the instances
		int featureDim = fvs[0].size();
		pKMeansCluster km(featureDim, numToSelect, fvs.size());
		km.randomizeKMeansPlusPlus(fvs);
		km.runToConvergenceOrMaxIterations(fvs, 10000000); // the large # here means "no limit on # of iterations to converge"
		
		// figure out how many instances are in each cluster, and find the closest instance to each cluster center
		vector<int> numInstancesInCluster(numToSelect, 0);
		vector<int> closestInstanceToClusterCenterBag(numToSelect, 0); // this stores the index of the bag containing the closest instance to the cluster center
		vector<float> closestInstancesToClusterCenterDist(numToSelect, FLT_MAX);
		
		for(int i = 0; i < mimlExs.size(); ++i)
		{
			for(int j = 0; j < mimlExs[i].instances_.size(); ++j)
			{
				float _distToClusterCenter;
				int clusterAssignment = km.getLabelForExampleAndReturnDistance(mimlExs[i].instances_[j], _distToClusterCenter);
				numInstancesInCluster[clusterAssignment]++;
				
				if( _distToClusterCenter <= closestInstancesToClusterCenterDist[clusterAssignment] )
				{
					closestInstancesToClusterCenterDist[clusterAssignment] = _distToClusterCenter;
					closestInstanceToClusterCenterBag[clusterAssignment] = i;
				}
			}
		}
		
		// sort the clusters descending by size (and carry the corresponding selected bag along with it)
		vector< pair<int, int> > clusterSizeAndSelectedBagPairs;
		for(int i = 0; i < numToSelect; ++i)
			clusterSizeAndSelectedBagPairs.push_back(pair<int, int>(numInstancesInCluster[i], closestInstanceToClusterCenterBag[i]));
		
		sort(clusterSizeAndSelectedBagPairs.rbegin(), clusterSizeAndSelectedBagPairs.rend());
		
		// unzip the selected bags
		vector<int> selectedBags;
		for(int i = 0; i < numToSelect; ++i)
			selectedBags.push_back(clusterSizeAndSelectedBagPairs[i].second);

		return selectedBags;
	}
	
	// similar to the coverage orcale method, but doesn't cheat. instead, it clusters the instaces with kmeans,
	// then applies the greedy max coverage algorithm to clusters instead of classes. one of the issues with this approach is
	// that it is hard to know in advance how many clusters to use, and one hopes that clusters and classes are well aligned.
	static vector<int> discoveryMethod_clusterCoverage(vector<pMIMLExample>& mimlExs, int numToSelect, int k)
	{
		// put the instances in a datastructure that can be used by k-means
		vector< vector<float> > fvs;
		for(int i = 0; i < mimlExs.size(); ++i)
			for(int j = 0; j < mimlExs[i].instances_.size(); ++j)
				fvs.push_back(mimlExs[i].instances_[j]);

		// run kmeans++ on the instances
		int featureDim = fvs[0].size();
		pKMeansCluster km(featureDim, k, fvs.size());
		km.randomizeKMeansPlusPlus(fvs);
		km.runToConvergenceOrMaxIterations(fvs, 10000000); // the large # here means "no limit on # of iterations to converge"
		
		// for each bag, get a list of clusters it contains
		vector< vector<int> > clusterSets;
		for(int i = 0; i < mimlExs.size(); ++i)
		{
			vector<int> clustersInThisBag;
			for(int j = 0; j < mimlExs[i].instances_.size(); ++j)
			{
				int instanceCluster = km.getLabelForExample(mimlExs[i].instances_[j]);
				if( !pVectorUtils::vector_contains(clustersInThisBag, instanceCluster) )
					clustersInThisBag.push_back(instanceCluster);
			}
			clusterSets.push_back(clustersInThisBag);
		}
		
		// select bags by greedy max-coverage of clusters
		vector<int> selectedBags;
		set<int> knownClusters; // these are all of the clusters coverage so far by selected bags
		
		for(int i = 0; i < numToSelect; ++i)
		{
			int bestBag = 0;
			int bestCoverage = 0;
			
			for(int j = 0; j < mimlExs.size(); ++j) // consider bag j as a candidate
			{
				set<int> newKnownClusters(knownClusters);
				
				for(int k = 0; k < clusterSets[j].size(); ++k)
					newKnownClusters.insert(clusterSets[j][k]);
				
				int newCoverage = newKnownClusters.size();
				
				if( newCoverage > bestCoverage )
				{
					bestCoverage = newCoverage;
					bestBag = j;
				}
			}
			
			if( bestCoverage == 0 )
			{
				cout << "best coverage = 0, num selected = " << i << endl;
			}
			
			selectedBags.push_back(bestBag);
			
			// update the set of covered clusters with the newly selected bag
			for(int k = 0; k < clusterSets[bestBag].size(); ++k)
				knownClusters.insert(clusterSets[bestBag][k]);
		}
			
		return selectedBags;
	}
	
	// same as cluster coverage, but breaks ties (multiple bags with the same coverage), using MIFF with k=1.
	// TODO: delete this function. it will be made unecessary after implementing the fully unified version of CC/MIFF
	static vector<int> discoveryMethod_clusterCoverage_tiebreakMIFF(vector<pMIMLExample>& mimlExs, int numToSelect, int numClusters)
	{
		// put the instances in a datastructure that can be used by k-means
		vector< vector<float> > fvs;
		for(int i = 0; i < mimlExs.size(); ++i)
			for(int j = 0; j < mimlExs[i].instances_.size(); ++j)
				fvs.push_back(mimlExs[i].instances_[j]);

		// run kmeans++ on the instances
		int featureDim = fvs[0].size();
		pKMeansCluster km(featureDim, numClusters, fvs.size());
		km.randomizeKMeansPlusPlus(fvs);
		km.runToConvergenceOrMaxIterations(fvs, 10000000); // the large # here means "no limit on # of iterations to converge"
		
		// for each bag, get a list of clusters it contains
		vector< vector<int> > clusterSets;
		for(int i = 0; i < mimlExs.size(); ++i)
		{
			vector<int> clustersInThisBag;
			for(int j = 0; j < mimlExs[i].instances_.size(); ++j)
			{
				int instanceCluster = km.getLabelForExample(mimlExs[i].instances_[j]);
				if( !pVectorUtils::vector_contains(clustersInThisBag, instanceCluster) )
					clustersInThisBag.push_back(instanceCluster);
			}
			clusterSets.push_back(clustersInThisBag);
		}
		
		// select bags by greedy max-coverage of clusters
		vector<int> selectedBags;
		set<int> selectedBagSet;
		set<int> knownClusters; // these are all of the clusters coverage so far by selected bags
		vector< vector<float> > coveredInstances; // these are the instances covered by selected bags. used by MIFF, not cluster coverage
		
		// keep going until we have selected the requested number of bags.
		// NOTE: once all clusters are covered, all of the remaining bags will be tied for coverage improvements, and MIFF will be the selection criteria
		while(selectedBags.size() < numToSelect)
		{
			int bestCoverage = 0;
			vector<int> bagsTiedForBest;
			
			for(int j = 0; j < mimlExs.size(); ++j) // consider bag j as a candidate
			{
				if( selectedBagSet.count(j) != 0 ) continue; // bag j has already been selected. skip it
				
				set<int> newKnownClusters(knownClusters);
				for(int k = 0; k < clusterSets[j].size(); ++k)
					newKnownClusters.insert(clusterSets[j][k]);
					
				int newCoverage = newKnownClusters.size();
				if( newCoverage > bestCoverage )
				{
					bestCoverage = newCoverage;
					bagsTiedForBest.clear();
					bagsTiedForBest.push_back(j);
				}
				else if( newCoverage == bestCoverage )
				{
					bagsTiedForBest.push_back(j);
				}
			}
			
			// find the bag that is tied for best (according to cluster coverage), with the instance that is furthest 
			// from its nearest neighbor in the set of covered instances
			vector<float> bagScores(mimlExs.size(), -1);
			for(int j = 0; j < bagsTiedForBest.size(); ++j)
			{
				int currBag = bagsTiedForBest[j];
				for(int m = 0; m < mimlExs[currBag].instances_.size(); ++m)
					bagScores[currBag] = max(bagScores[currBag], distanceToNearestNeighbor(mimlExs[currBag].instances_[m], coveredInstances));
			} 
			
			int bestBag = pVectorUtils::argmax(bagScores);
			
			// update list and set of selected bags
			selectedBags.push_back(bestBag);
			selectedBagSet.insert(bestBag);
			
			// update the set of covered clusters with the newly selected bag
			for(int k = 0; k < clusterSets[bestBag].size(); ++k)
				knownClusters.insert(clusterSets[bestBag][k]);
				
			// update the instances covered by selected bags
			for(int k = 0; k < mimlExs[bestBag].instances_.size(); ++k)
				coveredInstances.push_back(mimlExs[bestBag].instances_[k]);
		}
		
		return selectedBags;
	}
	
	// applies cluster coverage and breaks ties with MIFF with any k. this is a full unification of cluster coverage and MIFF
	static vector<int> discoveryMethod_unified_clusterCoverage_MIFF(vector<pMIMLExample>& mimlExs, int numToSelect, int numClusters, int miffK)
	{
		// put the instances in a datastructure that can be used by k-means
		vector< vector<float> > fvs;
		for(int i = 0; i < mimlExs.size(); ++i)
			for(int j = 0; j < mimlExs[i].instances_.size(); ++j)
				fvs.push_back(mimlExs[i].instances_[j]);

		// run kmeans++ on the instances
		int featureDim = fvs[0].size();
		pKMeansCluster km(featureDim, numClusters, fvs.size());
		km.randomizeKMeansPlusPlus(fvs);
		km.runToConvergenceOrMaxIterations(fvs, 10000000); // the large # here means "no limit on # of iterations to converge"
		
		// for each bag, get a list of clusters it contains
		vector< vector<int> > clusterSets;
		for(int i = 0; i < mimlExs.size(); ++i)
		{
			vector<int> clustersInThisBag;
			for(int j = 0; j < mimlExs[i].instances_.size(); ++j)
			{
				int instanceCluster = km.getLabelForExample(mimlExs[i].instances_[j]);
				if( !pVectorUtils::vector_contains(clustersInThisBag, instanceCluster) )
					clustersInThisBag.push_back(instanceCluster);
			}
			clusterSets.push_back(clustersInThisBag);
		}
		
		// select bags by greedy max-coverage of clusters
		vector<int> selectedBags;
		set<int> selectedBagSet;
		set<int> knownClusters; // these are all of the clusters coverage so far by selected bags
		vector< vector<float> > coveredInstances; // these are the instances covered by selected bags. used by MIFF, not cluster coverage
		
		// keep going until we have selected the requested number of bags.
		// NOTE: once all clusters are covered, all of the remaining bags will be tied for coverage improvements, and MIFF will be the selection criteria
		while(selectedBags.size() < numToSelect)
		{
			if( knownClusters.size() == numClusters )
				cout << "all clusters covered, bags selected = " << selectedBags.size() << endl;
			else
				cout << knownClusters.size() << " clusters covered, bags selected = " << selectedBags.size() << endl;
				
			int bestCoverage = 0;
			vector<int> bagsTiedForBest;
			
			for(int j = 0; j < mimlExs.size(); ++j) // consider bag j as a candidate
			{
				if( selectedBagSet.count(j) != 0 ) continue; // bag j has already been selected. skip it
				
				set<int> newKnownClusters(knownClusters);
				for(int k = 0; k < clusterSets[j].size(); ++k)
					newKnownClusters.insert(clusterSets[j][k]);
					
				int newCoverage = newKnownClusters.size();
				if( newCoverage > bestCoverage )
				{
					bestCoverage = newCoverage;
					bagsTiedForBest.clear();
					bagsTiedForBest.push_back(j);
				}
				else if( newCoverage == bestCoverage )
				{
					bagsTiedForBest.push_back(j);
				}
			}
			
			cout << "# bags tied for best coverage = " << bagsTiedForBest.size() << endl;
			
			// find the bag that is tied for best (according to cluster coverage), with the instance that is furthest 
			// from its nearest neighbor in the set of covered instances
			vector<float> bagScores(mimlExs.size());
			
			// fill in the bag scores with random negative numbers. at the end of considering all bags,
			// the highest scoring bag will be selected. however, any bag with 0 instances will be skipped.
			// putting random negative numbers in here causes the algorithm to degenerate into "random" if all other
			// options are exhausted.
			for(int j = 0; j < bagScores.size(); ++j) 
				bagScores[j] = -pRandom::randInRange(1,2); // the range [-2,-1] is used here to be totally sure there is no floating point issue near 0
			
			for(int j = 0; j < bagsTiedForBest.size(); ++j)
			{
				int currBag = bagsTiedForBest[j];
				if( mimlExs[currBag].instances_.size() == 0 ) continue; // if the bag has 0 instances, skip it
				
				vector< vector<float> > coveredInstancesInThisBag;
				hash_set<int> coveredInstancesInThisBagSet;
				
				// find the k "best" instances in the bag. the bag's score for selection
				// is the sum of the distances from those best k instances to their nearest covered neighbor
				float bagScore = 0;
				
				for(int l = 0; l < miffK; ++l)
				{
					float maxDistToNearestCoveredNeighbor = 0;
					int bestInstanceToCover = 0;
					
					for(int m = 0; m < mimlExs[currBag].instances_.size(); ++m)
					{
						if( coveredInstancesInThisBagSet.count(m) != 0 ) continue;
						
						// considering instance m as one of the best in the bag, its distance to its nearest covered neighbor
						// is the min of the distance to its nearest covered neighbor in already selected bags, and other instances in the same bag
						float distFromInstanceMToNearestCoveredNeighbor = 
							min(distanceToNearestNeighbor(mimlExs[currBag].instances_[m], coveredInstances), 
								distanceToNearestNeighbor(mimlExs[currBag].instances_[m], coveredInstancesInThisBag));
								
						if( distFromInstanceMToNearestCoveredNeighbor >= maxDistToNearestCoveredNeighbor )
						{
							maxDistToNearestCoveredNeighbor = distFromInstanceMToNearestCoveredNeighbor;
							bestInstanceToCover = m;
						}
					}
					
					coveredInstancesInThisBag.push_back(mimlExs[currBag].instances_[bestInstanceToCover]);
					coveredInstancesInThisBagSet.insert(bestInstanceToCover);
					bagScore += maxDistToNearestCoveredNeighbor;
				}
				
				bagScores[currBag] = bagScore;
			} 
			
			int bestBag = pVectorUtils::argmax(bagScores);
			
			// update list and set of selected bags
			selectedBags.push_back(bestBag);
			selectedBagSet.insert(bestBag);
			
			// update the set of covered clusters with the newly selected bag
			for(int k = 0; k < clusterSets[bestBag].size(); ++k)
				knownClusters.insert(clusterSets[bestBag][k]);
				
			// update the instances covered by selected bags
			for(int k = 0; k < mimlExs[bestBag].instances_.size(); ++k)
				coveredInstances.push_back(mimlExs[bestBag].instances_[k]);
		}
		
		return selectedBags;
	}
	
	// solves a generalized version of the max coverage problem to promote cluster diversity
	static vector<int> discoveryMethod_generalizedMaxCoverage(vector<pMIMLExample>& mimlExs, int numToSelect)
	{
		// put the instances in a datastructure that can be used by k-means
		cout << "flattening intances" << endl;
		vector< vector<float> > fvs;
		for(int i = 0; i < mimlExs.size(); ++i)
			for(int j = 0; j < mimlExs[i].instances_.size(); ++j)
				fvs.push_back(mimlExs[i].instances_[j]);
		
		// run kmeans++ on the instances
		cout << "clustering instances" << endl;
		int featureDim = fvs[0].size();
		int k = 50; // TODO: this is just a placeholder value, in the absence of a better way to pick k
		pKMeansCluster km(featureDim, k, fvs.size());
		km.randomizeKMeansPlusPlus(fvs);
		km.runToConvergenceOrMaxIterations(fvs, 10000000); // the large # here means "no limit on # of iterations to converge"
		
		// for each bag, get a list of clusters it contains
		vector< vector<int> > clusterSets;
		
		// at the same time...
		// compute hausdorff distance between each pair of clusters.
		// the interface for computing hausdorff distance expects a pMIMLExample,
		// so we will construct one for each cluster
		vector<pMIMLExample> clustersAsBags(k);
	
		cout << "clusterfying" << endl;
		for(int i = 0; i < mimlExs.size(); ++i)
		{
			vector<int> clustersInThisBag;
			
			for(int j = 0; j < mimlExs[i].instances_.size(); ++j)
			{
				int instanceCluster = km.getLabelForExample(mimlExs[i].instances_[j]);
				
				clustersAsBags[instanceCluster].instances_.push_back(mimlExs[i].instances_[j]);
				
				if( !pVectorUtils::vector_contains(clustersInThisBag, instanceCluster) )
					clustersInThisBag.push_back(instanceCluster);
			}
			
			clusterSets.push_back(clustersInThisBag);
		}
		
		// compute a distance measure between each pair of clusters
		cout << "computing distances between cluster pairs" << endl;
		vector< vector<float> > clusterPairDistances(k, vector<float>(k, 0));
		for(int i = 0; i < k; ++i)
		for(int j = 0; j < k; ++j)
			//clusterPairDistances[i][j] = pBagDistance::closestInstancePair(clustersAsBags[i], clustersAsBags[j]);
			//clusterPairDistances[i][j] = pBagDistance::farthestInstancePair(clustersAsBags[i], clustersAsBags[j]);
			clusterPairDistances[i][j] = pBagDistance::DHavg(clustersAsBags[i], clustersAsBags[j]);
		
		
		// apply the greedy algorithm for generalized max-coverage
		vector<int> selectedBags;
		set<int> coveredClusters;
		
		// start out with the single bag which achieves the best objective.
		double currObjective = 0;
		int firstBag = 0;
		for(int i = 0; i < mimlExs.size(); ++i) // consider selecting bag i first
		{
			double firstBagObjective = 0;
			set<int> firstCoveredClusters;
			for(int j = 0; j < clusterSets[i].size(); ++j)
				firstCoveredClusters.insert(clusterSets[i][j]);
			
			for(set<int>::iterator p = firstCoveredClusters.begin(); p != firstCoveredClusters.end(); ++p)
			for(set<int>::iterator q = firstCoveredClusters.begin(); q != firstCoveredClusters.end(); ++q)
				firstBagObjective += clusterPairDistances[*p][*q];
			
			if( firstBagObjective >= currObjective )
			{
				currObjective = firstBagObjective;
				firstBag = i;
			}
		}
		
		selectedBags.push_back(firstBag);
		for(int i = 0; i < clusterSets[firstBag].size(); ++i)
			coveredClusters.insert(clusterSets[firstBag][i]);
		
		// select the rest of the bags greedily
		for(int i = 1; i < numToSelect; ++i)
		{
			double bestObjective = 0;
			int bestBag = 0;
			cout << "selecting bag " << i << endl;
			
			cout << "# covered clusters = " << coveredClusters.size() << endl; 
			if( coveredClusters.size() == k ) cout << "warning: covered all clusters" << endl;
			
			for(int j = 0; j < mimlExs.size(); ++j) // consider selecting bag j
			{
				// figure out which clusters are covered as a result of selecting bag i, which were not already covered before
				vector<int> newlyCoveredClusters;
				for(int l = 0; l < clusterSets[j].size(); ++l)
				{
					if( coveredClusters.count(clusterSets[j][l]) == 0 )
						newlyCoveredClusters.push_back(clusterSets[j][l]);
				}
				
				// the objective achieved by selecting bag j is the sum of
				// (a) the objective before selecting bag j
				double newObjective = currObjective;
				
				// (b) the distance between all pairs of clusters (p,q) such that p was already covered and q is newly covered
				for(set<int>::iterator p = coveredClusters.begin(); p != coveredClusters.end(); ++p)
				for(int q = 0; q < newlyCoveredClusters.size(); ++q)
					newObjective += clusterPairDistances[*p][newlyCoveredClusters[q]];
				
				// (c) the distance between all pairs of newly covered clusters
				for(int p = 0; p < newlyCoveredClusters.size(); ++p)
				for(int q = 0; q < newlyCoveredClusters.size(); ++q)
					newObjective += clusterPairDistances[newlyCoveredClusters[p]][newlyCoveredClusters[q]];
				
				if( newObjective >= bestObjective )
				{
					bestObjective = newObjective;
					bestBag = j;
				}
			}
			
			// add bestBag to the selections, and update cluster coverage
			currObjective = bestObjective;
			selectedBags.push_back(bestBag);
			for(int j = 0; j < clusterSets[bestBag].size(); ++j)
				coveredClusters.insert(clusterSets[bestBag][j]);
			
			cout << "objective = " << currObjective << ", log(objective)=" << log(currObjective) << endl;
		}
		
		return selectedBags;
	}
	
	// selects bags to maximize the diversity of the covered instances. this can be viewed as a form
	// of generalized max-coverage. it should (hopefully) be better than discoveryMethod_instance_farthest,
	// because it considers the value of all instances that are covered by each bag rather than the value of a single instance per bag.
	// this method was motivated by the obsevation that cluster coverage methods tend to cover all clusters before finding all classes,
	// leaving some classes undiscovered with no further room for improvement
	static vector<int> discoveryMethod_instance_coverage(vector<pMIMLExample>& mimlExs, int numToSelect)
	{
		vector<int> selectedBags;
		set<int> selectedBagSet;
		vector< vector<float> > coveredInstances;
	
		// start with a random bag. it must not contain 0 instances!
		int randomBag;
		do { randomBag = rand()%mimlExs.size(); } while( mimlExs[randomBag].instances_.size() == 0);
		
		selectedBags.push_back(randomBag);
		selectedBagSet.insert(randomBag);
		for(int j = 0; j < mimlExs[randomBag].instances_.size(); ++j)
			coveredInstances.push_back(mimlExs[randomBag].instances_[j]);
		
		// the objective after selecting one bag is just the sum of pairwise distances between its instances
		double currDistanceSum = 0;
		for(int p = 0; p < mimlExs[randomBag].instances_.size(); ++p)
		for(int q = p + 1; q < mimlExs[randomBag].instances_.size(); ++q)
			currDistanceSum += distL2Squared(mimlExs[randomBag].instances_[p], mimlExs[randomBag].instances_[q]);
		
		//double currObjective = currDistanceSum / float(coveredInstances.size() * (coveredInstances.size()-1));
		double currObjective = currDistanceSum;
		
		cout<< "randomBag=" << randomBag << " # inst=" << mimlExs[randomBag].instances_.size() << endl;
		
		for(int i = 1; i < numToSelect; ++i)
		{
			cout << "selecting bag " << i << ", # instances covered = " << coveredInstances.size() << ", objective=" << currObjective << endl;
			
			int bestBag = 0;
			double bestObjective = 0;
			double bestNewDistanceSum = 0;
			
			for(int j = 0; j < mimlExs.size(); ++j) // consider selecting bag j
			{
				if( selectedBagSet.count(j) > 0 ) continue; // dissalow re-selecting the same bag
				if( mimlExs[j].instances_.size() == 0 ) continue; // don't query empty bags
				
				// TODO: it might be possible to improve the asymptotic complexity using an iterator over a set
				
				// selecting bag j adds all of its instances to the set of covered instances.
				// the new objective with bag j is the sum of:
				// (a) the objective before adding bag j
				double newDistanceSum = currDistanceSum;
				
				// (b) the distance between all pairs of instances in bag j
				for(int p = 0; p < mimlExs[j].instances_.size(); ++p)
				for(int q = p + 1; q < mimlExs[j].instances_.size(); ++q)
					newDistanceSum += distL2Squared(mimlExs[j].instances_[p], mimlExs[j].instances_[q]);
				
				// (c) the distance between all pairs of instances (p,q) such that p is in bag j and q was already covered
				for(int p = 0; p < mimlExs[j].instances_.size(); ++p)
				for(int q = 0; q < coveredInstances.size(); ++q)
					newDistanceSum += distL2Squared(mimlExs[j].instances_[p], coveredInstances[q]);
				
				//double newObjective = newDistanceSum  / double((coveredInstances.size() + mimlExs[j].instances_.size()) * (coveredInstances.size() + mimlExs[j].instances_.size() - 1));
				double newObjective = newDistanceSum;
				
				if( newObjective >= bestObjective )
				{
					bestObjective = newObjective;
					bestNewDistanceSum = newDistanceSum;
					bestBag = j;
				}
			}
			
			// save the best bag
			selectedBags.push_back(bestBag);
			selectedBagSet.insert(bestBag);
			
			// the new objective is the objective achieved by the best bag
			currObjective = bestObjective;
			currDistanceSum = bestNewDistanceSum;
			
			// cover all of the instances from the best bag
			for(int j = 0; j < mimlExs[bestBag].instances_.size(); ++j)
				coveredInstances.push_back(mimlExs[bestBag].instances_[j]);
		}
		
		return selectedBags;
	}
	
	// selects bags by apply a farthest-first traversal to instances, and querying the bags the come from
	static vector<int> discoveryMethod_instance_farthest(vector<pMIMLExample>& mimlExs, int numToSelect)
	{
		vector< vector<float> > fvs; // flattened instance feature vectors
		map<int, int> instanceToBag; // maps the index of an instance to the index of the bag it comes from
		
		for(int i = 0; i < mimlExs.size(); ++i)
		{
			for(int j = 0; j < mimlExs[i].instances_.size(); ++j)
			{
				fvs.push_back(mimlExs[i].instances_[j]);
				instanceToBag[fvs.size()-1] = i;
			}
		}
		
		// note: it is necessary to select more instances than bags, because otherwise we may select the same bag twice.
		// without intertwining the instance-level farthest first code, it is not trivial to select just the right amount. 
		// as a hack, it will just select twice as many instances as bags, then later check if that was enough (and break if it isn't).
		int numInstancesToSelect = numToSelect * 4;
		vector<int> selectedInstances = pFarthestFirstTraversal::traverseMaxMin(fvs, numInstancesToSelect);
		
		// select the bags from which the farthest instances came, but skip any bag that has already been selected

		vector<int> selectedBags;
		set<int> bagSet;
		int currCandidateInstance = 0;
		while(selectedBags.size() < numToSelect && currCandidateInstance < selectedInstances.size())
		{
			int candidateBag = instanceToBag[selectedInstances[currCandidateInstance]];
			
			if( bagSet.count(candidateBag) == 0 )
			{
				bagSet.insert(candidateBag);
				selectedBags.push_back(candidateBag);
			}
			
			++currCandidateInstance;
		}
		
		if( selectedBags.size() != numToSelect )
		{
			cout << "error: not enough instances were selected to obtain the specified number of unique bags. increase numInstancesToSelect and try again" << endl;
			cout << "got " << selectedBags.size() << " bags" << endl;
			assert(false);
		}
		
		return selectedBags;
	}
	
	// the basic instance farthest first algorithm does not take into account the issue that by querying a particular instance,
	// all other instances in that bag will also be covered. this second version of instance farthest first repeatedly selects
	// the instance that is farthest from all instances covered by currently selected bags
	static vector<int> discoveryMethod_instance_farthest_v2(vector<pMIMLExample>& mimlExs, int numToSelect)
	{
		// start by selecting a random non-empty bag
		vector<int> selectedBags;
		hash_set<int> selectedBagSet;
		vector<int> nonEmptyBags;
		for(int i = 0; i < mimlExs.size(); ++i)
			if( mimlExs[i].instances_.size() != 0 ) 
				nonEmptyBags.push_back(i);
		int firstBag = nonEmptyBags[rand()%nonEmptyBags.size()];
		selectedBags.push_back(firstBag);
		selectedBagSet.insert(firstBag);
		
		for(int i = 1; i < numToSelect; ++i)
		{
			cout << "selecting bag " << i << endl;
			
			vector<float> candidateDistances(mimlExs.size(), -1); // -1 will never be the argmax over other choices
			
			for(int j = 0; j < mimlExs.size(); ++j)
			{
				if( selectedBagSet.count(j) != 0 ) continue; // check if bag j was already selected and skip it if so
				
				for(int k = 0; k < mimlExs[j].instances_.size(); ++k)
				{
					// after the following loop, this will store the distance from instance k in bag j to the closest instance covered so far
					float minDistToSet = FLT_MAX; 
					
					for(int p = 0; p < selectedBags.size(); ++p)
					for(int q = 0; q < mimlExs[selectedBags[p]].instances_.size(); ++q)
						minDistToSet = min(minDistToSet, distL2Squared(mimlExs[selectedBags[p]].instances_[q], mimlExs[j].instances_[k]));
					
					// for the purpose of selecting a bag, what we care about is the "best" instance in it, 
					// i.e. the one with the maximum distance to its nearest neighbor in the coevered set.
					// this update says that the distance associated with bag j is the largest among its instance distances
					candidateDistances[j] = max(candidateDistances[j], minDistToSet);
				}
			}
			
			int bagToSelect = pVectorUtils::argmax(candidateDistances);
			selectedBags.push_back(bagToSelect);
			selectedBagSet.insert(bagToSelect);
		}
		
		return selectedBags;
	}
	
	// generalization of v2... TODO: comment this better. 
	// NOTE: it is likely that this method will be discarded in favor of a more general/better versino of it
	static vector<int> discoveryMethod_instance_farthest_v3(vector<pMIMLExample>& mimlExs, int numToSelect, int k)
	{
		// start by selecting a non-empty bag
		vector<int> selectedBags;
		hash_set<int> selectedBagSet;
		vector<int> nonEmptyBags;
		for(int i = 0; i < mimlExs.size(); ++i)
			if( mimlExs[i].instances_.size() != 0 ) 
				nonEmptyBags.push_back(i);
		int firstBag = nonEmptyBags[rand()%nonEmptyBags.size()];
		selectedBags.push_back(firstBag);
		selectedBagSet.insert(firstBag);

		for(int i = 1; i < numToSelect; ++i)
		{
			cout << "selecting bag " << i << endl;
			
			vector<float> candidateDistances(mimlExs.size(), -1);
			
			for(int j = 0; j < mimlExs.size(); ++j) // consider selecting bag j
			{
				if( selectedBagSet.count(j) != 0 ) continue; // check if bag j was already selected and skip it if so
				
				// for every instance l in bag j, compute the distance to its nearest neighbor in the set of already covered instances
				vector<float> distFromInstToNearestCoveredNeighbor;
				for(int l = 0; l < mimlExs[j].instances_.size(); ++l)
				{
					float minDistToSet = FLT_MAX; 
					for(int p = 0; p < selectedBags.size(); ++p)
					for(int q = 0; q < mimlExs[selectedBags[p]].instances_.size(); ++q)
						minDistToSet = min(minDistToSet, distL2Squared(mimlExs[selectedBags[p]].instances_[q], mimlExs[j].instances_[l]));
					distFromInstToNearestCoveredNeighbor.push_back(minDistToSet);
				}
				
				// we want the k "best" instances from bag j to represent it, i.e. the k instances that are farthest from 
				// their nearest neighbor in the covered set. therefore, we sort the distances in descending order, then take the first k
				sort(distFromInstToNearestCoveredNeighbor.begin(), distFromInstToNearestCoveredNeighbor.end(), std::greater<float>());
				
				// compute the value of selecting bag j as the sum of the distances from the covered set to its best k instances
				float sumBestKDists = 0;
				for(int l = 0; l < min(k, (int)distFromInstToNearestCoveredNeighbor.size()); ++l) // if there are fewer than k instances in the bag, just take all of them
					sumBestKDists += distFromInstToNearestCoveredNeighbor[l];
					
				candidateDistances[j] = sumBestKDists;
			}
			
			int bagToSelect = pVectorUtils::argmax(candidateDistances);
			selectedBags.push_back(bagToSelect);
			selectedBagSet.insert(bagToSelect);
		}
		
		return selectedBags;
	}
	
	// get the distance from x to its nearest neighbor in fvs. this is a helper function for MIFF
	static float distanceToNearestNeighbor(vector<float>& x, vector< vector<float> >& fvs)
	{
		float minDist = FLT_MAX;
		for(int i = 0; i < fvs.size(); ++i)
			minDist = min(minDist, distL2Squared(fvs[i], x));
		return minDist;
	}
	
	// multi-instance farthest first (MIFF). the parameter k controls how many instances within each bag are considered.
	// when k=1, this method works as follows: first a random non-empty bag is selected, and all of its instances are covered.
	// next, it selects the bag which contains the instance who's nearest neighbor in the set of already covered instances is maximized.
	// when k>1, the criteria for which bag to select is a sum of distances instead of a single distance. the first term in that sum is the
	// same distance described above; the rest of the terms of distances from other instances in the same bag to the covered set, in a greedy 
	// fashion (instances within the bag are covered one at a time, up to k instances, and their distance to the covered set includes both
	// previously selected bags, and previously selected instances in the same bag).
	static vector<int> discoveryMethod_MIFF(vector<pMIMLExample>& mimlExs, int numToSelect, int k)
	{
		// start by selecting a random non-empty bag
		vector<int> selectedBags;
		hash_set<int> selectedBagSet;
		vector< vector<float> > coveredInstances;
		
		vector<int> nonEmptyBags;
		for(int i = 0; i < mimlExs.size(); ++i)
			if( mimlExs[i].instances_.size() != 0 ) 
				nonEmptyBags.push_back(i);
		int firstBag = nonEmptyBags[rand()%nonEmptyBags.size()];
		selectedBags.push_back(firstBag);
		selectedBagSet.insert(firstBag);
		
		for(int i = 0; i < mimlExs[firstBag].instances_.size(); ++i)
			coveredInstances.push_back(mimlExs[firstBag].instances_[i]);
		
		for(int i = 1; i < numToSelect; ++i)
		{
			cout << "selecting bag " << i << endl;
			
			// after considering all eligible bags, we will select the bag with the highest score
			vector<float> candidateScores(mimlExs.size(), -1);
			
			for(int j = 0; j < mimlExs.size(); ++j) // consider selecting bag j
			{
				if( selectedBagSet.count(j) != 0 ) continue; // check if bag j was already selected and skip it if so
				if( mimlExs[j].instances_.size() == 0 ) continue; // if the bag has 0 instances, skip it
	
				vector< vector<float> > coveredInstancesInThisBag;
				hash_set<int> coveredInstancesInThisBagSet;
				
				// find the k "best" instances in the bag. the bag's score for selection
				// is the sum of the distances from those best k instances to their nearest covered neighbor
				float bagScore = 0;
				for(int l = 0; l < k; ++l)
				{
					float maxDistToNearestCoveredNeighbor = 0;
					int bestInstanceToCover = 0;
					
					for(int m = 0; m < mimlExs[j].instances_.size(); ++m)
					{
						if( coveredInstancesInThisBagSet.count(m) != 0 ) continue; // don't consider instances that are already covered
						
						// considering instance m as one of the best in the bag, its distance to its nearest covered neighbor
						// is the min of the distance to its nearest covered neighbor in already selected bags, and other instances in the same bag
						float distFromInstanceMToNearestCoveredNeighbor = 
							min(distanceToNearestNeighbor(mimlExs[j].instances_[m], coveredInstances), 
								distanceToNearestNeighbor(mimlExs[j].instances_[m], coveredInstancesInThisBag));
						
						if( distFromInstanceMToNearestCoveredNeighbor >= maxDistToNearestCoveredNeighbor )
						{
							maxDistToNearestCoveredNeighbor = distFromInstanceMToNearestCoveredNeighbor;
							bestInstanceToCover = m;
						}
					}
					
					coveredInstancesInThisBag.push_back(mimlExs[j].instances_[bestInstanceToCover]);
					coveredInstancesInThisBagSet.insert(bestInstanceToCover);
					bagScore += maxDistToNearestCoveredNeighbor;
				}
				
				candidateScores[j] = bagScore;
			}
			
			// select the bag with the highest score, and update all necessary data structures (i.e. covered instances)
			int bagToSelect = pVectorUtils::argmax(candidateScores);
			selectedBags.push_back(bagToSelect);
			selectedBagSet.insert(bagToSelect);
			
			for(int i = 0; i < mimlExs[bagToSelect].instances_.size(); ++i)
				coveredInstances.push_back(mimlExs[bagToSelect].instances_[i]);
		}
		
		return selectedBags;
	}

	// first apply the cluster coverage method, then once all clusters are covered, apply MIFF
	static vector<int> discoveryMethod_clusterCoverageThenMIFF(vector<pMIMLExample>& mimlExs, int numToSelect, int numClusters, int miffK) 
	{
		// put the instances in a datastructure that can be used by k-means
		vector< vector<float> > fvs;
		for(int i = 0; i < mimlExs.size(); ++i)
			for(int j = 0; j < mimlExs[i].instances_.size(); ++j)
				fvs.push_back(mimlExs[i].instances_[j]);
				
		// run kmeans++ on the instances
		int featureDim = fvs[0].size();
		pKMeansCluster km(featureDim, numClusters, fvs.size());
		km.randomizeKMeansPlusPlus(fvs);
		km.runToConvergenceOrMaxIterations(fvs, 10000000); // the large # here means "no limit on # of iterations to converge"
		
		// for each bag, get a list of clusters it contains
		vector< vector<int> > clusterSets;
		for(int i = 0; i < mimlExs.size(); ++i)
		{
			vector<int> clustersInThisBag;
			for(int j = 0; j < mimlExs[i].instances_.size(); ++j)
			{
				int instanceCluster = km.getLabelForExample(mimlExs[i].instances_[j]);
				if( !pVectorUtils::vector_contains(clustersInThisBag, instanceCluster) )
					clustersInThisBag.push_back(instanceCluster);
			}
			clusterSets.push_back(clustersInThisBag);
		}
		
		// select bags by greedy max-coverage of clusters
		vector<int> selectedBags;
		set<int> selectedBagSet;
		set<int> knownClusters; // these are all of the clusters coverage so far by selected bags
		vector< vector<float> > coveredInstances; // these are the instances covered by selected bags. used by MIFF, not cluster coverage
		
		// keep going until either we have selected the requested number of bags, or all clusters are covered
		while(selectedBags.size() < numToSelect && knownClusters.size() != numClusters)
		{
			int bestCoverage = 0;
			vector<int> bagsTiedForBest;
			
			for(int j = 0; j < mimlExs.size(); ++j) // consider bag j as a candidate
			{
				if( selectedBagSet.count(j) != 0 ) continue; // bag j has already been selected. skip it
				
				set<int> newKnownClusters(knownClusters);
				for(int k = 0; k < clusterSets[j].size(); ++k)
					newKnownClusters.insert(clusterSets[j][k]);
					
				int newCoverage = newKnownClusters.size();
				if( newCoverage > bestCoverage )
				{
					bestCoverage = newCoverage;
					bagsTiedForBest.clear();
					bagsTiedForBest.push_back(j);
				}
				else if( newCoverage == bestCoverage )
				{
					bagsTiedForBest.push_back(j);
				}
			}
			
			// find the bag that is tied for best (according to cluster coverage), with the instance that is furthest 
			// from its nearest neighbor in the set of covered instances
			vector<float> bagScores(mimlExs.size(), -1);
			for(int j = 0; j < bagsTiedForBest.size(); ++j)
			{
				int currBag = bagsTiedForBest[j];
				for(int m = 0; m < mimlExs[currBag].instances_.size(); ++m)
					bagScores[currBag] = max(bagScores[currBag], distanceToNearestNeighbor(mimlExs[currBag].instances_[m], coveredInstances));
			} 
			
			int bestBag = pVectorUtils::argmax(bagScores);
			
			// update list and set of selected bags
			selectedBags.push_back(bestBag);
			selectedBagSet.insert(bestBag);
			
			// update the set of covered clusters with the newly selected bag
			for(int k = 0; k < clusterSets[bestBag].size(); ++k)
				knownClusters.insert(clusterSets[bestBag][k]);
				
			// update the instances covered by selected bags
			for(int k = 0; k < mimlExs[bestBag].instances_.size(); ++k)
				coveredInstances.push_back(mimlExs[bestBag].instances_[k]);
		}
		
		cout << "finished cluster coverage, got " << selectedBags.size() << " bags" << endl;
		
		// select the rest of the bags with MIFF
		while(selectedBags.size() < numToSelect)
		{
			cout << "selecting next bag by MIFF, # bags selected = " << selectedBags.size() << endl;
			
			// after considering all eligible bags, we will select the bag with the highest score
			vector<float> candidateScores(mimlExs.size(), -1);
			
			for(int j = 0; j < mimlExs.size(); ++j) // consider selecting bag j
			{
				if( selectedBagSet.count(j) != 0 ) continue; // check if bag j was already selected and skip it if so
				if( mimlExs[j].instances_.size() == 0 ) continue; // if the bag has 0 instances, skip it
				
				vector< vector<float> > coveredInstancesInThisBag;
				hash_set<int> coveredInstancesInThisBagSet;
			
				// find the k "best" instances in the bag. the bag's score for selection
				// is the sum of the distances from those best k instances to their nearest covered neighbor
				float bagScore = 0;
				for(int l = 0; l < miffK; ++l)
				{
					float maxDistToNearestCoveredNeighbor = 0;
					int bestInstanceToCover = 0;
					
					for(int m = 0; m < mimlExs[j].instances_.size(); ++m)
					{
						if( coveredInstancesInThisBagSet.count(m) != 0 ) continue; // don't consider instances that are already covered
						
						// considering instance m as one of the best in the bag, its distance to its nearest covered neighbor
						// is the min of the distance to its nearest covered neighbor in already selected bags, and other instances in the same bag
						float distFromInstanceMToNearestCoveredNeighbor = 
							min(distanceToNearestNeighbor(mimlExs[j].instances_[m], coveredInstances), 
								distanceToNearestNeighbor(mimlExs[j].instances_[m], coveredInstancesInThisBag));
						
						if( distFromInstanceMToNearestCoveredNeighbor >= maxDistToNearestCoveredNeighbor )
						{
							maxDistToNearestCoveredNeighbor = distFromInstanceMToNearestCoveredNeighbor;
							bestInstanceToCover = m;
						}
					}
					
					coveredInstancesInThisBag.push_back(mimlExs[j].instances_[bestInstanceToCover]);
					coveredInstancesInThisBagSet.insert(bestInstanceToCover);
					bagScore += maxDistToNearestCoveredNeighbor;
				}
				
				candidateScores[j] = bagScore;
			}
			
			// select the bag with the highest score, and update all necessary data structures (i.e. covered instances)
			int bagToSelect = pVectorUtils::argmax(candidateScores);
			selectedBags.push_back(bagToSelect);
			selectedBagSet.insert(bagToSelect);
			
			for(int i = 0; i < mimlExs[bagToSelect].instances_.size(); ++i)
				coveredInstances.push_back(mimlExs[bagToSelect].instances_[i]);
		}
		
		return selectedBags;
	}
	
	static void runExperiment(string discoveryMethod, int numToSelect, int numRepetitions, bool useBootstrap)
	{
		vector<string> speciesCodes;
		vector<string> speciesNames;
		vector<pAudioFile_MetaData> metaData;
		
		pMetaData_Parser::parseSpeciesList("species_list.txt", speciesCodes, speciesNames);
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
		pMetaData_Parser::parseBagLabelSets("bag_labels.txt", metaData);
		
		
		const int featureDim = pMetaData_Parser::getSegmentFeatureDim(metaData);
		int numClasses = speciesNames.size();
		int numBags = metaData.size();
		
		cout << "# casses = " << numClasses << endl;
		cout << "# bags = " << numBags << endl;
	
		cout << "rescaling features" << endl;
		pMetaData_Rescale::rescaleTo01(metaData, featureDim);
		
		// build a MIML dataset, so we can pass the data to class discovery algorithms in an abstract format
		vector<pMIMLExample> mimlExs;
		for(int i = 0; i < metaData.size(); ++i)
		{
			vector< vector<float> > instances;
			for(int j = 0; j < metaData[i].segments_.size(); ++j)
				instances.push_back(metaData[i].segments_[j].featureVector_);
			
			mimlExs.push_back(pMIMLExample(instances, metaData[i].speciesList_));
		}
		
		// run multiple randomized repetitions
		
		vector< vector<int> > classDiscoveryCurves;
		
		for(int rep = 0; rep < numRepetitions; ++rep)
		{
			cout << "*** repetition " << rep << endl;
			
			vector<pMIMLExample> mimlExsForThisRep;
			if( useBootstrap ) // this might not technically be bootstrapping since it is without replacement
			{
				cout << "bootstrapping" << endl;
				vector<int> allBagIndices;
				for(int i = 0; i < mimlExs.size(); ++i)
					allBagIndices.push_back(i);
				pRandom::shuffle(allBagIndices);
				
				const int numToBootstrap = float(mimlExs.size()) * .63;
				
				for(int i = 0; i < numToBootstrap; ++i)
					mimlExsForThisRep.push_back(mimlExs[allBagIndices[i]]); 
			}
			else
			{
				mimlExsForThisRep = mimlExs;
			}
			
			numToSelect = min(numToSelect, (int)mimlExsForThisRep.size());
			
			// apply the chosen class discovery method
			
			vector<int> selectedBags;
			
			if		( discoveryMethod == "random" )					selectedBags = discoveryMethod_random					(mimlExsForThisRep, numToSelect);
			else if ( discoveryMethod == "random_non_empty" )		selectedBags = discoveryMethod_randomNonEmpty			(mimlExsForThisRep, numToSelect);
			else if	( discoveryMethod == "farthest" )				selectedBags = discoveryMethod_instance_farthest		(mimlExsForThisRep, numToSelect);
			else if	( discoveryMethod == "farthest2" )				selectedBags = discoveryMethod_instance_farthest_v2		(mimlExsForThisRep, numToSelect);
			else if	( discoveryMethod == "farthest3" )				selectedBags = discoveryMethod_instance_farthest_v3		(mimlExsForThisRep, numToSelect, 3); // the last argument here is k, the number of instances to consider per bag
			else if	( discoveryMethod == "MIFFk1" )					selectedBags = discoveryMethod_MIFF						(mimlExsForThisRep, numToSelect, 1); // the last argument here is k, the number of instances to consider per bag
			else if	( discoveryMethod == "MIFFk2" )					selectedBags = discoveryMethod_MIFF						(mimlExsForThisRep, numToSelect, 2); // the last argument here is k, the number of instances to consider per bag
			else if	( discoveryMethod == "CCMIFF" )					selectedBags = discoveryMethod_clusterCoverageThenMIFF	(mimlExsForThisRep, numToSelect, 50, 1);
			else if	( discoveryMethod == "coverage_oracle" )		selectedBags = discoveryMethod_coverageOracle			(mimlExsForThisRep, numToSelect);
			else if	( discoveryMethod == "cluster_coverage" )		selectedBags = discoveryMethod_clusterCoverage			(mimlExsForThisRep, numToSelect, 50);
			else if	( discoveryMethod == "generalized_coverage" )	selectedBags = discoveryMethod_generalizedMaxCoverage	(mimlExsForThisRep, numToSelect);
			else if	( discoveryMethod == "instance_coverage" )		selectedBags = discoveryMethod_instance_coverage		(mimlExsForThisRep, numToSelect);
			else if	( discoveryMethod == "cluster_centers" )		selectedBags = discoveryMethod_clusterCenters			(mimlExsForThisRep, numToSelect);
			
			
			else assert(false && "unrecognized discovery method");
			
			assert(selectedBags.size() == numToSelect);	
			
			cout << "selected bags: " << endl;
			pVectorUtils::cout_vector(selectedBags); cout << endl;
			
			
			/// generate class discovery curve ///
			
			vector<int> numClassesFound; // numClassesFound[i] is the number of classes found after looking at i bags
			set<int> knownClasses;
			
			numClassesFound.push_back(0);
			
			for(int i = 0; i < selectedBags.size(); ++i)
			{
				for(int j = 0; j < mimlExsForThisRep[selectedBags[i]].labels_.size(); ++j)
					knownClasses.insert(mimlExsForThisRep[selectedBags[i]].labels_[j]);
				
				numClassesFound.push_back(knownClasses.size());
			}
			
			classDiscoveryCurves.push_back(numClassesFound);
		}
		
		vector<float> avgClassDiscoveryCurve;
		vector<float> stdevClassDiscoveryCurve;
		for(int i = 0; i < numToSelect; ++i)
		{
			float avg = 0;
			for(int rep = 0; rep < numRepetitions; ++rep)
				avg += classDiscoveryCurves[rep][i];
			avg /= float(numRepetitions);
			
			avgClassDiscoveryCurve.push_back(avg);
			
			float sse = 0;
			for(int rep = 0; rep < numRepetitions; ++rep)
				sse += (avg - classDiscoveryCurves[rep][i]) * (avg - classDiscoveryCurves[rep][i]);
			
			float stdev = sqrt(sse / float(numRepetitions));
			stdevClassDiscoveryCurve.push_back(stdev);
		}
		
		/// output results ///
		for(int i = 0; i < numToSelect; ++i)
			cout << i << "\t" << avgClassDiscoveryCurve[i] << "\t" << stdevClassDiscoveryCurve[i] << endl;
	}
	
	// this function does the work of one thread for runExperiment_multithread
	static void* runExperiment_pthread_run(void* ptr)
	{
		pMIMLClassDiscoveryExperiment_runExperiment_args& args = *((pMIMLClassDiscoveryExperiment_runExperiment_args*) ptr);
		
		cout << "runExperiment_pthread_run, method=" << args.discoveryMethod_ << " param1=" << args.param1_ << " param2=" << args.param2_ << endl;
		
		bool done = false;
		while(!done)
		{
			vector<pMIMLExample> mimlExsForThisRep;

			cout << "thread " << args.threadNum_ << " bootstrapping" << endl;
			vector<int> allBagIndices;
			for(int i = 0; i < (*args.mimlExs_).size(); ++i)
				allBagIndices.push_back(i);
			pRandom::shuffle(allBagIndices);
			
			const int numToBootstrap = float((*args.mimlExs_).size()) * .63;
			
			for(int i = 0; i < numToBootstrap; ++i)
				mimlExsForThisRep.push_back((*args.mimlExs_)[allBagIndices[i]]); 

			int numToSelect = min(args.numToSelect_, (int)mimlExsForThisRep.size());
			
			// apply the chosen class discovery method
			
			vector<int> selectedBags;
			
			if		( args.discoveryMethod_ == "random" )						selectedBags = discoveryMethod_random						(mimlExsForThisRep, args.numToSelect_);
			else if ( args.discoveryMethod_ == "random_non_empty" )				selectedBags = discoveryMethod_randomNonEmpty				(mimlExsForThisRep, args.numToSelect_);
			else if	( args.discoveryMethod_ == "farthest" )						selectedBags = discoveryMethod_instance_farthest			(mimlExsForThisRep, args.numToSelect_);
			else if	( args.discoveryMethod_ == "farthest2" )					selectedBags = discoveryMethod_instance_farthest_v2			(mimlExsForThisRep, args.numToSelect_);
			else if	( args.discoveryMethod_ == "farthest3" )					selectedBags = discoveryMethod_instance_farthest_v3			(mimlExsForThisRep, args.numToSelect_, args.param1_);
			else if	( args.discoveryMethod_ == "MIFF" )							selectedBags = discoveryMethod_MIFF							(mimlExsForThisRep, args.numToSelect_, args.param1_);
			else if	( args.discoveryMethod_ == "CCMIFF" )						selectedBags = discoveryMethod_clusterCoverageThenMIFF		(mimlExsForThisRep, args.numToSelect_, args.param1_, args.param2_);
			else if	( args.discoveryMethod_ == "coverage_oracle" )				selectedBags = discoveryMethod_coverageOracle				(mimlExsForThisRep, args.numToSelect_);
			else if	( args.discoveryMethod_ == "cluster_coverage" )				selectedBags = discoveryMethod_clusterCoverage				(mimlExsForThisRep, args.numToSelect_, args.param1_);
			else if	( args.discoveryMethod_ == "cluster_coverage_tiebreak" )	selectedBags = discoveryMethod_clusterCoverage_tiebreakMIFF	(mimlExsForThisRep, args.numToSelect_, args.param1_);
			else if	( args.discoveryMethod_ == "unified_CCMIFF" )				selectedBags = discoveryMethod_unified_clusterCoverage_MIFF	(mimlExsForThisRep, args.numToSelect_, args.param1_, args.param2_);
			else if	( args.discoveryMethod_ == "generalized_coverage" )			selectedBags = discoveryMethod_generalizedMaxCoverage		(mimlExsForThisRep, args.numToSelect_);
			else if	( args.discoveryMethod_ == "instance_coverage" )			selectedBags = discoveryMethod_instance_coverage			(mimlExsForThisRep, args.numToSelect_);
			else if	( args.discoveryMethod_ == "cluster_centers" )				selectedBags = discoveryMethod_clusterCenters				(mimlExsForThisRep, args.numToSelect_);
			
			else assert(false && "unrecognized discovery method");
			
			assert(selectedBags.size() == numToSelect);	
			
			cout << "thread " << args.threadNum_ << " selected bags: " << endl;
			pVectorUtils::cout_vector(selectedBags); cout << endl;
			
			/// generate class discovery curve ///
			
			vector<int> numClassesFound; // numClassesFound[i] is the number of classes found after looking at i bags
			set<int> knownClasses;
			
			numClassesFound.push_back(0);
			
			for(int i = 0; i < selectedBags.size(); ++i)
			{
				for(int j = 0; j < mimlExsForThisRep[selectedBags[i]].labels_.size(); ++j)
					knownClasses.insert(mimlExsForThisRep[selectedBags[i]].labels_[j]);
				
				numClassesFound.push_back(knownClasses.size());
			}
		
			// this thread finished generating a class discovery curve. 
			// if there are enough already, just discard it. otherwise, add it to the results datastructure
			pthread_mutex_lock(&pMIMLClassDiscoveryExperiment_mutex);
			if( args.classDiscoveryCurves_->size() >= args.numRepetitions_ )
				done = true;
			else												
				args.classDiscoveryCurves_->push_back(numClassesFound);
			pthread_mutex_unlock(&pMIMLClassDiscoveryExperiment_mutex);
		}
		
		return NULL;
	}
	
	// run many bootstrapped experiments in parallel
	static void runExperiment_multithread(string discoveryMethod, string resultsFilename, int numToSelect, int numRepetitions, int param1 = 0, int param2 = 0)
	{
		vector<string> speciesCodes;
		vector<string> speciesNames;
		vector<pAudioFile_MetaData> metaData;
		
		pMetaData_Parser::parseSpeciesList("species_list.txt", speciesCodes, speciesNames);
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
		pMetaData_Parser::parseBagLabelSets("bag_labels.txt", metaData);
		
		
		const int featureDim = pMetaData_Parser::getSegmentFeatureDim(metaData);
		int numClasses = speciesNames.size();
		int numBags = metaData.size();
		
		cout << "# casses = " << numClasses << endl;
		cout << "# bags = " << numBags << endl;
	
		cout << "rescaling features" << endl;
		pMetaData_Rescale::rescaleTo01(metaData, featureDim);
		
		// build a MIML dataset, so we can pass the data to class discovery algorithms in an abstract format
		vector<pMIMLExample> mimlExs;
		for(int i = 0; i < metaData.size(); ++i)
		{
			vector< vector<float> > instances;
			for(int j = 0; j < metaData[i].segments_.size(); ++j)
				instances.push_back(metaData[i].segments_[j].featureVector_);
			
			mimlExs.push_back(pMIMLExample(instances, metaData[i].speciesList_));
		}
		
		// results are stored here before final analysis
		vector< vector<int> > classDiscoveryCurves;
		
		// setup the threads
		//int numThreads = sysconf( _SC_NPROCESSORS_ONLN ); // set the number of threads to the number of processors
		int numThreads = 10;
		vector<pthread_t> threads(numThreads);
		vector<pMIMLClassDiscoveryExperiment_runExperiment_args>threadArgs (numThreads);
		
		for(int t = 0; t < numThreads; ++t)
		{
			pMIMLClassDiscoveryExperiment_runExperiment_args& args = threadArgs[t];
			
			args.threadNum_ = t;
			args.discoveryMethod_ = discoveryMethod;
			args.numToSelect_ = numToSelect;
			args.numRepetitions_ = numRepetitions;
			args.mimlExs_ = &mimlExs;
			args.param1_ = param1;
			args.param2_ = param2;
			args.classDiscoveryCurves_ = &classDiscoveryCurves;
			
			if( t != 0 )  // the 0'th thread will be this thread (i.e. the main thread). this reduces wasted time
				pthread_create(&threads[t], NULL, runExperiment_pthread_run, (void*)&args);
		}
		
		// before waiting for other threads to finish, make some trees in this thread as well
		runExperiment_pthread_run(&threadArgs[0]);
		
		for(int t = 0; t < numThreads; ++t)
			pthread_join(threads[t], NULL);
		
		cout << "(all threads finished)" << endl;
		
		vector<float> avgClassDiscoveryCurve;
		vector<float> stdevClassDiscoveryCurve;
		for(int i = 0; i < numToSelect; ++i)
		{
			float avg = 0;
			for(int rep = 0; rep < numRepetitions; ++rep)
				avg += classDiscoveryCurves[rep][i];
			avg /= float(numRepetitions);
			
			avgClassDiscoveryCurve.push_back(avg);
			
			float sse = 0;
			for(int rep = 0; rep < numRepetitions; ++rep)
				sse += (avg - classDiscoveryCurves[rep][i]) * (avg - classDiscoveryCurves[rep][i]);
			
			float stdev = sqrt(sse / float(numRepetitions));
			stdevClassDiscoveryCurve.push_back(stdev);
		}
		
		/// output results ///
		pTextFile resultsFile("species_discovery_results/" + resultsFilename, PFILE_WRITE);
		resultsFile << "num_selected" << "\t" << discoveryMethod << "_" << param1 << "_" << param2 << "\n";
		for(int i = 0; i < numToSelect; ++i)
		{
			cout		<< i << "\t" << avgClassDiscoveryCurve[i] << "\t" << stdevClassDiscoveryCurve[i] << endl;
			resultsFile << i << "\t" << avgClassDiscoveryCurve[i] << "\n";
		}	
		resultsFile.close();
	}
};