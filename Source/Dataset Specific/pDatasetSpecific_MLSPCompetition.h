// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// this file implements stuff I used to setup/test the dataset for the MLSP Kaggle competition.

#pragma once

#include "pMetaData.h"
#include "pTextFile.h"
#include "pTrainValTest.h"
#include "pRandom.h"
#include "pExample_MultiLabel.h"
#include "pBinaryRelevance.h"
#include "pDatasetSpecific_HJAndrews.h"

#include <set>
using namespace std;

class pDatasetSpecific_MLSPCompetition
{
public:
	// this is used to generate some figures in the MLSP many-author paper. The idea is to see which dates, times, and sites are covered by the 
	// recordings in the dataset
	static void visualizeRecordingCoverage()
	{
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		
		vector<int> recordingsByHour(24, 0);
		map<string, int> recordingsByDateStr;
		map<string, int> recordingsBySite;
		
		for(int i = 0; i < metaData.size(); ++i)
		{
			cout << metaData[i].filename_ << endl;
			
			string site = pDatasetSpecific_HJAndrews::getSiteFromFilename(metaData[i].filename_);
			int hour = pDatasetSpecific_HJAndrews::getHourFromFilename(metaData[i].filename_);
			pair<int,int> date = pDatasetSpecific_HJAndrews::getDateFromFilename(metaData[i].filename_);
			int month = date.first;
			int day = date.second;
			int year = pDatasetSpecific_HJAndrews::getYearFromFilename(metaData[i].filename_);
			
			cout << "\t" << "site: "	<< site << endl;
			cout << "\t" << "year: "	<< year << endl;
			cout << "\t" << "month: "	<< month << endl;
			cout << "\t" << "day: "		<< day << endl;
			cout << "\t" << "hour: "	<< hour << endl;
			
			string dateStr = pStringUtils::intToString(month) + "/" + pStringUtils::intToString(day) + "/" + pStringUtils::intToString(year-2000);
			
			++recordingsByHour[hour];
			
			if( recordingsByDateStr.count(dateStr) == 0 )
				recordingsByDateStr[dateStr] = 1;
			else
				recordingsByDateStr[dateStr] = recordingsByDateStr[dateStr] + 1;
				
			if( recordingsBySite.count(site) == 0 )
				recordingsBySite[site] = 1;
			else
				recordingsBySite[site] = recordingsBySite[site] + 1;
		}
		
		// output recordings per hour
		cout << "---recordings per hour----" << endl;
		for(int i = 0; i < 24; ++i)
			cout << i << "," << recordingsByHour[i] << endl;
		
		// output recordings by date string
		cout << "---recordings by date str---" << endl;
		for(map<string,int>::iterator it = recordingsByDateStr.begin(); it != recordingsByDateStr.end(); ++it)
			cout << (*it).first << "," << (*it).second << endl;
			
		// output recordings by site
		cout << "---recordings by site---" << endl;
		for(map<string,int>::iterator it = recordingsBySite.begin(); it != recordingsBySite.end(); ++it)
			cout << (*it).first << "," << (*it).second << endl;
	}
	
	// baseline-method binary relevance random forest on histogram of segments features
	static void baselineRFBRHisto()
	{
		// load the dataset
		vector<string> speciesCodes;
		vector<string> speciesNames;
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseSpeciesList("species_list.txt", speciesCodes, speciesNames);
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseBagLabelSets("bag_labels.txt", metaData);
		pMetaData_Parser::parseCrossValidationFolds("CVfolds_2.txt", metaData);
		pMetaData_Parser::parseBagLevelFeatures("histogram_of_segments.txt", metaData);
		
		int numClasses = speciesCodes.size();
		cout << "#classes = " << numClasses << endl;
		cout << "#bags = " << metaData.size() << endl;
		
		vector<pExample_MultiLabel> trainingExs;
		for(int i = 0; i < metaData.size(); ++i)
			if( metaData[i].crossValidationFold_ == kTrain )
				trainingExs.push_back(pExample_MultiLabel(metaData[i].bagLevelFeatureVector_, metaData[i].speciesList_));
				
		const int numTrees = 100;
		pBinaryRelevance br;
		br.train(trainingExs, numClasses, numTrees);
		
		pTextFile f("brrfhisto_predictions.txt", PFILE_WRITE);
		f << "rec_id,species,probability\n";
		
		for(int i = 0; i < metaData.size(); ++i)
		{
			if( metaData[i].crossValidationFold_ == kTest )
			{
				
				vector<float> scores = br.getClassScores(metaData[i].bagLevelFeatureVector_);
				
				for(int j = 0; j < numClasses; ++j)
					f << i << "," << j << "," << scores[j] << "\n";
	
			}
		}
		f.close();
	}
	
	// output the desired results in the format required to setup a competition on Kaggle.com
	static void makeSolutionFileForKaggle()
	{
		// load the dataset
		vector<string> speciesCodes;
		vector<string> speciesNames;
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseSpeciesList("species_list.txt", speciesCodes, speciesNames);
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseBagLabelSets("bag_labels.txt", metaData);
		pMetaData_Parser::parseCrossValidationFolds("CVfolds_2.txt", metaData);
		int numClasses = speciesCodes.size();
		cout << "#classes = " << numClasses << endl;
		cout << "#bags = " << metaData.size() << endl;
		
		pTextFile f("kaggle_solution.txt", PFILE_WRITE);
		f << "rec_id,species,probability,Usage\n";
		
		
		// put 33% of the test set in the "public test set"
		vector<int> testIndices;
		for(int i = 0; i < metaData.size(); ++i)
			if( metaData[i].crossValidationFold_ == kTest )
				testIndices.push_back(i);
			
		pRandom::shuffle(testIndices);
		
		set<int> publicTestSet;
		for(int i = 0; i < testIndices.size() / 3; ++i)
			publicTestSet.insert(testIndices[i]);
		
		for(int i = 0; i < metaData.size(); ++i)
		{
			if( metaData[i].crossValidationFold_ == kTest )
			{
				for(int j = 0; j < numClasses; ++j)
					f << i << "," << j << "," << (pVectorUtils::vector_contains(metaData[i].speciesList_, j) ? 1 : 0) << "," << (publicTestSet.count(i) == 1 ? "PublicTest" : "PrivateTest") << "\n";
			}
		}
		
		f.close();

	}

	// make an example file in the format that should be submitted
	static void makeSampleSubmissionFile()
	{
		// load the dataset
		vector<string> speciesCodes;
		vector<string> speciesNames;
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseSpeciesList("species_list.txt", speciesCodes, speciesNames);
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseBagLabelSets("bag_labels.txt", metaData);
		pMetaData_Parser::parseCrossValidationFolds("CVfolds_2.txt", metaData);
		int numClasses = speciesCodes.size();
		cout << "#classes = " << numClasses << endl;
		cout << "#bags = " << metaData.size() << endl;
		
		pTextFile f("sample_submission.txt", PFILE_WRITE);
		f << "bag_id,species,probability\n";
		
		for(int i = 0; i < metaData.size(); ++i)
		{
			if( metaData[i].crossValidationFold_ == kTest )
			{
				for(int j = 0; j < numClasses; ++j)
					f << i << "," << j << ",0\n";
			}
		}
		
		f.close();
	}
	
	// make a bag labels file with the test set labels hidden / replaced with a ?
	static void makeBagLabelsTestHidden()
	{
		// load the dataset
		vector<string> speciesCodes;
		vector<string> speciesNames;
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseSpeciesList("species_list.txt", speciesCodes, speciesNames);
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseBagLabelSets("bag_labels.txt", metaData);
		pMetaData_Parser::parseCrossValidationFolds("CVfolds_2.txt", metaData);
		int numClasses = speciesCodes.size();
		cout << "#classes = " << numClasses << endl;
		cout << "#bags = " << metaData.size() << endl;
		
		pTextFile f("bag_labels_test_hidden.txt", PFILE_WRITE);
		f << "bag_id,[labels]\n";
		
		for(int i = 0; i < metaData.size(); ++i)
		{
			f << i;
			if( metaData[i].crossValidationFold_ == kTrain )
			{
				for(int j = 0; j < metaData[i].speciesList_.size(); ++j)
					f << "," << metaData[i].speciesList_[j];
			}
			else
			{
				f << ",?"; // if its part of the test set, put a ? instead of a label set
			}
			
			f << "\n";
		}
		
		f.close();
	}

	static void makeTrainAndTestSets()
	{
		// load the dataset
		vector<string> speciesCodes;
		vector<string> speciesNames;
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseSpeciesList("species_list.txt", speciesCodes, speciesNames);
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseBagLabelSets("bag_labels.txt", metaData);
		int numClasses = speciesCodes.size();
		cout << "#classes = " << numClasses << endl;
		cout << "#bags = " << metaData.size() << endl;
		
		// split the dataset 50-50 into train / test, randomly
		vector<pTrainValTest> partition;
		for(int i = 0; i < metaData.size(); ++i)
			partition.push_back(i < metaData.size() / 2 ? kTrain : kTest);
		pRandom::templatedShuffle<pTrainValTest>(partition);
		
		// write a 2-fold CV file
		pTextFile cvFoldsFile("CVfolds_2.txt", PFILE_WRITE);
		cvFoldsFile << "bag_id,fold" << "\n";
		for(int i = 0 ; i < metaData.size(); ++i)
			cvFoldsFile << i << "," << int(partition[i]) << "\n";
		cvFoldsFile.close();
		
		// copy a few files from the training set to /segmentation_examples
		vector<int> trainIndices;
		for(int i = 0; i < metaData.size(); ++i)
			if( partition[i] == kTrain )
				trainIndices.push_back(i);
		pRandom::shuffle(trainIndices);
		
		for(int i = 0; i < 20; ++i) // select this many examples for segmentation labeling
		{
			pAudioFile_MetaData& md = metaData[trainIndices[i]];
			cout << md.filename_ << endl;
			
			string copyCommand = "cp ./spectrograms/" + md.filename_ + " ./segmentation_examples/";
			system(copyCommand.c_str());
		}
	}

	// mostly what this does is prune the original 40-class dataset down to an 18 class dataset
	// the point is to remove recordings which contain a species of which there is only one example
	// and to remove some labels that aren't birds (e.g. wind) which are in the original label set
	static void filterOriginalLabels()
	{
		// load the dataset
		vector<string> speciesCodes;
		vector<string> speciesNames;
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseSpeciesList("species_list.txt", speciesCodes, speciesNames);
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseBagLabelSets("bag_labels.txt", metaData);
		int numClasses = speciesCodes.size();
		cout << "#classes = " << numClasses << endl;
		cout << "#bags = " << metaData.size() << endl;
		
		// figure out how many examples there are of each class
		vector<int> examplesPerClass(numClasses);
		for(int i = 0; i < metaData.size(); ++i)
			for(int j = 0; j < metaData[i].speciesList_.size(); ++j)
				++examplesPerClass[metaData[i].speciesList_[j]];
		
		for(int j = 0; j < numClasses; ++j)
			cout << speciesNames[j] << "\t" << examplesPerClass[j] << endl;
			
		// some of the non-bird labels should be removed from the label set (e.g. wind, airplane)
		vector<int> classesToRemove;
		classesToRemove.push_back(21);		// wind
		classesToRemove.push_back(22);		// stream
		classesToRemove.push_back(23);		// airplane
		classesToRemove.push_back(24);		// rain
		classesToRemove.push_back(26);		// unknown bird
		
		// add any classes with <=2 example to the list of classes to remove
		for(int j = 0; j < numClasses; ++j)
			if( examplesPerClass[j] <=2 && !pVectorUtils::vector_contains(classesToRemove, j) )
				classesToRemove.push_back(j);
				
		// make a new list of classes with the gaps from removed classes compressed
		vector<string> filteredSpeciesCodes;
		vector<string> filteredSpeciesNames;
		map<int,int> classRemap;
		for(int j = 0; j < numClasses; ++j)
		{
			if( !pVectorUtils::vector_contains(classesToRemove, j) )
			{
				classRemap[j] = filteredSpeciesCodes.size();
				filteredSpeciesCodes.push_back(speciesCodes[j]);
				filteredSpeciesNames.push_back(speciesNames[j]);
			}
		}
		
		int filteredNumClasses = filteredSpeciesCodes.size();
		cout << "#filtered classes = " << filteredNumClasses << endl;
		
		// write an updated species index file
		pTextFile filteredClassFile("filtered_species_list.txt", PFILE_WRITE);
		
		filteredClassFile << "class_id,code,species" << "\n";
		
		for(int j = 0; j < filteredNumClasses; ++j)
			filteredClassFile << j << "," << filteredSpeciesCodes[j] << "," << filteredSpeciesNames[j] << "\n";
		filteredClassFile.close();
		
			
		// run through the dataset a second time and decide to accept or reject each bag based on its labels.
		// while doing this, write a new bag label set file, and a new bag_id2filename.txt
		pTextFile filteredBagLabelsFile("filtered_bag_labels.txt", PFILE_WRITE);
		pTextFile filteredBagIDFile("filtered_bag_id2filename.txt", PFILE_WRITE);
		filteredBagLabelsFile << "bag_id,[labels]\n";
		filteredBagIDFile << "bag_id,filename\n";
		
		int filtered_bag_id = 0;
		int numRejected = 0;
		for(int i = 0; i < metaData.size(); ++i)
		{
			// if the label set contains a species which occurs in 2 or fewer recordings, we will reject it
			bool containsARareSpecies = false;
			for(int j = 0; j < metaData[i].speciesList_.size(); ++j)
				if( examplesPerClass[metaData[i].speciesList_[j]] <= 2 )
					containsARareSpecies = true;
					
			if( containsARareSpecies ) 
			{
				++numRejected;
				continue;
			}
						
			vector<int> filteredLabelSet;
			for(int j = 0; j < metaData[i].speciesList_.size(); ++j)
			{		
				if( !pVectorUtils::vector_contains(classesToRemove, metaData[i].speciesList_[j]) )
					filteredLabelSet.push_back(classRemap[metaData[i].speciesList_[j]]);
			}
			
			filteredBagLabelsFile << filtered_bag_id << (filteredLabelSet.size() == 0 ? "," : "");
			for(int j = 0; j < filteredLabelSet.size(); ++j)
				filteredBagLabelsFile << "," << filteredLabelSet[j]; 
			filteredBagLabelsFile << "\n";
			
			filteredBagIDFile << filtered_bag_id << "," << metaData[i].filename_ << "\n";
			
			string copyCommand = "cp ./src_wavs/" + pStringUtils::firstPartOfSplit(metaData[i].filename_, ".") + ".wav ./filtered_wavs/";
			system(copyCommand.c_str()); 
			
			++filtered_bag_id;
		}
		
		filteredBagLabelsFile.close(); 
		filteredBagIDFile.close();
		
		cout << "num rejected = " << numRejected << endl;
	}
};