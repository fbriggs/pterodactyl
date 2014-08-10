// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// this is the "everything but the kitchen sink" version of my Random Forest.
// this implementation of random forest has a few options beyond Brieman's basic algorithm:
// (1) a maximum depth can be set for trees
// (2) you can specify whether trees store a histogram or the the majority label, in cases where max-depth or identical features force a leaf to be constructed from an impure set of examples
// additionally, this implementation has a basic single-threaded training method which constructs trees one at a time, and a parallel version which attempts to use 1 thread per available processor
// to construct trees.
// (3) out-of-bag estimation is not implemented inside this class, but the necessary data to do it is exposed so it can be implemented elsewhere
// (4) class-balanced-bootstrapping: bootstrapping is problematic for classification problems with extremely unbalanced classes and/or very few examples of some classes. class-balanced-bootstrapping attempts to help with this issue by randomly constructing an approximately class-balanced dataset

#pragma once

#include "pVectorUtils.h"
#include "pRFDecisionTree_pthread.h"
#include "pTextFile.h"
#include <vector>
#include <iostream>
#include <pthread.h>

using namespace std;

pthread_mutex_t pRandomForest_pthread_mutex = PTHREAD_MUTEX_INITIALIZER;

void* pRandomForest_pthread_run(void* ptr); // this function does the work of one thread

struct pRandomForest_pthread_args
{
	int threadNum_;
	int numTreesToMake_;
	int numExamples_;
	int numClasses_;
	int featureDim_;
	int numFeaturesToTry_;
	bool storeHistogramInLeaf_;
	int maxDepth_;
	bool classBalancedBootstrapping_;
	vector<pRFDecisionTree_pthread*>* trees_;
	vector<pExample>* trainingExamples_;
	vector< vector<bool> >* inBag_; 
};

class pRandomForest_pthread
{
public:
	vector<pRFDecisionTree_pthread*> trees_;		// the ensemble of trees
	int numClasses_;								// the number of classes in the classification problem
	int featureDim_;								// dimension of the feature vectors
	
	// one of the nice features of RF is out-of-bag (OOB) estimation. pRandomForest_pthread does not implement all of these capabilities,
	// but exposes the necessary infomration to implement OOB features externally. inBag_[i][j] = true if tree i in the forest got example j in its bagged sample of the dataset
	// TODO: inBag_ is currently not saved/loaded to file; it can only be used in the same process that trains the RF.
	vector< vector<bool> > inBag_;

	// before training on an input dataset, the RF will check the data for basic errors. the issues detected are:
	// NAN values
	// INF values
	// inconsistent feature vector dimensions
	// warning: no examples of one class
	static void verifyDataset(int numClasses, vector<pExample>& trainingExamples)
	{
		assert(trainingExamples.size() > 0);
		
		int dim = trainingExamples[0].featureVector_.size();
		
		vector<int> examplesPerClass(numClasses, 0);
		
		for(int i = 0; i < trainingExamples.size(); ++i)
		{
			assert(dim == trainingExamples[i].featureVector_.size());
				
			++examplesPerClass[trainingExamples[i].classLabel_];
				
			for(int j = 0; j < trainingExamples[i].featureVector_.size(); ++j)
			{
				assert( !isnan(trainingExamples[i].featureVector_[j]) );
				assert( !isinf(trainingExamples[i].featureVector_[j]) );
			}
		}
		
		for(int j = 0; j < numClasses; ++j)
		{
			if( examplesPerClass[j] == 0 )
				cout << "WARNING: 0 training examples for class " << j << endl;
		}
	}

	void train_multithread(int numTrees, int numClasses, vector<pExample>& trainingExamples, bool storeHistogramInLeaf = true, int maxDepth = 15, bool useClassBalancedBootstrapping = false)
	{
		verifyDataset(numClasses, trainingExamples);
		
		int numThreads = sysconf( _SC_NPROCESSORS_ONLN ); // set the number of threads to the number of processors
		
		numClasses_ = numClasses;
		featureDim_ = trainingExamples[0].featureVector_.size();
		int numFeaturesToTry = (int) log2(featureDim_) + 1; // the number of features to decide between in each node of the decision trees.
		int numExamples = trainingExamples.size();
		
		inBag_ = vector< vector<bool> >(numTrees, vector<bool>(trainingExamples.size(), false));
		
		vector<pthread_t> threads(numThreads);
		vector<pRandomForest_pthread_args>threadArgs (numThreads);
		for(int t = 0; t < numThreads; ++t)
		{
			pRandomForest_pthread_args& args = threadArgs[t];
			
			args.threadNum_ = t;
			args.numExamples_ = numExamples;
			args.numClasses_ = numClasses;
			args.featureDim_ = featureDim_;
			args.numFeaturesToTry_ = numFeaturesToTry;
			args.numTreesToMake_ = numTrees;
			args.storeHistogramInLeaf_ = storeHistogramInLeaf;
			args.maxDepth_ = maxDepth;
			args.trees_ = &trees_;
			args.trainingExamples_ = &trainingExamples;
			args.inBag_ = &inBag_;
			args.classBalancedBootstrapping_ = useClassBalancedBootstrapping;
			
			if( t != 0 )  // the 0'th thread will be this thread (i.e. the main thread). this reduces wasted time
				pthread_create(&threads[t], NULL, pRandomForest_pthread_run, (void*)&args);
		}
		
		// before waiting for other threads to finish, make some trees in this thread as well
		pRandomForest_pthread_run(&threadArgs[0]);
		
		for(int t = 0; t < numThreads; ++t)
			pthread_join(threads[t], NULL);
		
		cout << "(all threads finished)" << endl;
	}
	
	void train(int numTrees, int numClasses, vector<pExample>& trainingExamples, bool storeHistogramInLeaf = true, int maxDepth = 15,  bool useClassBalancedBootstrapping = false)
	{
		verifyDataset(numClasses, trainingExamples);
		
		numClasses_ = numClasses;
		featureDim_ = trainingExamples[0].featureVector_.size();
		int numFeaturesToTry = (int) log2(featureDim_) + 1; // the number of features to decide between in each node of the decision trees.
		int numExamples = trainingExamples.size();
		
		inBag_ = vector< vector<bool> >(numTrees, vector<bool>(trainingExamples.size(), false));
		
		// build the decision tree ensemble
		for(int t = 0; t < numTrees; ++t)
		{
			cout << "*"; cout.flush();
			
			// construct a different bootstrap sample from the training set for each tree
			vector<pExample*> bootstrappedExamples;
			
			if( useClassBalancedBootstrapping ) // "class-balanced bootstrapping" - choose a class uniformly at random, then pick an example from that class iid with replacement
			{
				// divide all of the examples up into separate classes (so its faster to sample from a particular class)
				vector< vector<pExample*> > examplesInEachClass(numClasses, vector<pExample*>());
				for(int i = 0; i < trainingExamples.size(); ++i)
					examplesInEachClass[trainingExamples[i].classLabel_].push_back(&trainingExamples[i]);
					
				for(int i = 0; i < numExamples; ++i) // bootstrapped sample- same size as original dataset
				{
					int randomClass = rand()%numClasses;
					int randomInstance = rand()%examplesInEachClass[randomClass].size();
					bootstrappedExamples.push_back(examplesInEachClass[randomClass][randomInstance]);
					inBag_[t][randomInstance] = true; // make a note that this instance is in the bag for tree i
				}
				
			}
			else // basic boostrapping - blindly sample iid with replacement
			{
				for(int i = 0; i < numExamples; ++i)
				{
					int randomInstance = rand() % numExamples;
					bootstrappedExamples.push_back(&trainingExamples[randomInstance]);
					inBag_[t][randomInstance] = true; // make a note that this instance is in the bag for tree i
				}
			}
			
			// make a tree from this bootstrapped sample and add it to the ensemble
			pRFDecisionTree_pthread* tree = new pRFDecisionTree_pthread(bootstrappedExamples, numClasses, featureDim_, numFeaturesToTry, storeHistogramInLeaf, 0, maxDepth);
			
			trees_.push_back(tree);
		}
	}
	
	// add up the output histograms for each tree
	vector<float> estimateClassProbabilities(vector<float>& feature)
	{
		assert(feature.size() == featureDim_ ); 
		
		vector<float> outClassProbabilities(numClasses_);
		
		for(int i = 0; i < trees_.size(); ++i)
		{
			vector<float> histo = trees_[i]->getClassHistogramForInput(feature);
			
			for(int c = 0; c < numClasses_; ++c)
				outClassProbabilities[c] += histo[c];
		}
		
		pVectorUtils::normalizeVectorToPDF(outClassProbabilities);
		
		return outClassProbabilities;
	}
	
	// predict the class label for a feature
	int classify(vector<float>& feature)
	{
		vector<float> classProbabilities = estimateClassProbabilities(feature);
		return pVectorUtils::argmax(classProbabilities);
	}
		
	~pRandomForest_pthread()
	{
		for(int i = 0; i < trees_.size(); ++i)
			delete trees_[i];
	}
	
	void save(string filename)
	{
		pTextFile f(filename, PFILE_WRITE);
		
		f << numClasses_ << "," << featureDim_ << "\n";
		
		for(int i = 0; i < trees_.size(); ++i)
		{
			f << trees_[i]->toString() << "\n";
		}
		
		f.close();
	}
	
	void load(string filename)
	{
		pTextFile f(filename, PFILE_READ);
		
		vector<string> parts = pStringUtils::split(f.readLine(),",");
		assert(parts.size() == 2 && "RF model file must have header: numClasses,featureDim");
		
		numClasses_ = pStringUtils::stringToInt(parts[0]);
		featureDim_ = pStringUtils::stringToInt(parts[1]);
		
		while(!f.eof())
		{
			string line = f.readLine();
			trees_.push_back(new pRFDecisionTree_pthread(line, numClasses_));
		}
		
		f.close();
	}
};

// this function is run in each thread
void* pRandomForest_pthread_run(void* ptr)
{
	pRandomForest_pthread_args& args = *((pRandomForest_pthread_args*) ptr);
	
	bool done = false;
	while(!done)
	{
		// construct a different bootstrap sample from the training set for each tree
		vector<pExample*> bootstrappedExamples;
		
		vector<bool> inBagForThisTree(args.numExamples_, false);
		
		if( args.classBalancedBootstrapping_ )
		{
			// divide all of the examples up into separate classes (so its faster to sample from a particular class)
			vector< vector<pExample*> > examplesInEachClass(args.numClasses_, vector<pExample*>());
			for(int i = 0; i < (*args.trainingExamples_).size(); ++i)
				examplesInEachClass[(*args.trainingExamples_)[i].classLabel_].push_back(&((*args.trainingExamples_)[i]));
				
			for(int i = 0; i < args.numExamples_; ++i) // bootstrapped sample- same size as original dataset
			{
				int randomClass = rand()%args.numClasses_;
				int randomInstance = rand()%examplesInEachClass[randomClass].size();
				bootstrappedExamples.push_back(examplesInEachClass[randomClass][randomInstance]);
				inBagForThisTree[randomInstance] = true; // make a note that this instance is in the bag for tree i
			}
		}
		else
		{
			for(int i = 0; i < args.numExamples_; ++i)
			{
				int randomInstance = rand() % args.numExamples_;
				bootstrappedExamples.push_back(&(*args.trainingExamples_)[randomInstance]);
				inBagForThisTree[randomInstance] = true;
			}
		}
		
		// make a tree from this bootstrapped sample and add it to the ensemble
		pRFDecisionTree_pthread* tree = new pRFDecisionTree_pthread(bootstrappedExamples, args.numClasses_, args.featureDim_, args.numFeaturesToTry_, args.storeHistogramInLeaf_, 0, args.maxDepth_);
		
		pthread_mutex_lock(&pRandomForest_pthread_mutex);
		
		if( args.trees_->size() >= args.numTreesToMake_ )	
		{
			done = true;
			// if we got here, a tree has been allocated, but it will never be used/saved to the RF's list of trees, which is where the deallocation is handled.
			// it should be deleted at this point, or there will be a memory leak
			delete tree;
		}
		else												
		{
			args.trees_->push_back(tree);
			(*args.inBag_)[args.trees_->size()-1] = inBagForThisTree;
		}
		
		pthread_mutex_unlock(&pRandomForest_pthread_mutex);
		cout << "*"; cout.flush();
	}
	
	return NULL;
}
