// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// this file contains functions spefically for processing data from H. J. Andrews (not general purpose bird sound analysis).
// note this file is intended to be "separable" from the main project, in so far as the rest of the code should all work without it.
// furthermore, it is not exposed through the command line. 

#pragma once

#include <vector>
#include <string>
#include <map>
#include <set>

#include "pStringUtils.h"
#include "pMIMLExample.h"
#include "pMIMLClassDiscoveryExperiment.h"


using namespace std;

class pDatasetSpecific_HJAndrews
{
public:
	
	// the MLSP dataset is a subset of a collection that was analyzed by Sarah Hadley, Jed Irvine and others,
	// which includes audio from 2009 and 2010. Some processing is needed to convert the data in the format it was in when they labeled it,
	// to the format used in the current version of Pterodactyl. also, for the MLSP dataset it was desireable to remove files with rare species
	// (for a cross-validated classification experiment). However, we don't want to do that for the species discovery problem!, so that 
	// data needs to be filtered differently for this problem.
	// NOTE: some simple modifications were done in Excel as well (i.e. subtract 1 from some 1-based indices in the original data). 
	static void prepareLabeled2009_2010data_for_species_discovery()
	{
		// these store the species in the original meta data (no filtering)
		vector<string> originalSpeciesCodes;
		vector<string> originalSpeciesNames;
		pMetaData_Parser::parseSpeciesList("original_meta_data/class_names_renum.txt", originalSpeciesCodes, originalSpeciesNames);
		
		// these store the species in the filtered version (the species file was edited manually to delete some)
		vector<string> speciesCodes;
		vector<string> speciesNames;
		pMetaData_Parser::parseSpeciesList("species_list.txt", speciesCodes, speciesNames);
		
		// construct a map from the original species indices to new species indices. -1 is used for species/classes that have been deleted
		map<int, int> classRemap;
		for(int i = 0; i < originalSpeciesNames.size(); ++i)
		{
			classRemap[i] = pVectorUtils::index_of(speciesCodes, originalSpeciesCodes[i]);
		
			cout << i << "\t" << classRemap[i] << endl;
		}
		
		// generate a new bag-labels file
		pTextFile originalBagLabels("original_meta_data/bag_labels.txt", PFILE_READ);
		originalBagLabels.readLine(); // skip header
		
		pTextFile newBagLabels("bag_labels.txt", PFILE_WRITE);
		newBagLabels << "bag_id,[labels]" << "\n";
		
		while(!originalBagLabels.eof())
		{
			vector<string> parts = pStringUtils::splitNonEmpty(originalBagLabels.readLine(), ",");
			
			int bag_id = pStringUtils::stringToInt(parts[0]) - 1;
			
			newBagLabels <<  bag_id;
			
			string restOfTheLine = "";
			for(int i = 1; i < parts.size(); ++i)
			{
				int newClassNum = classRemap[pStringUtils::stringToInt(parts[i])-1];
				
				if( newClassNum != -1 )
					restOfTheLine += "," + pStringUtils::intToString(newClassNum);
			}
			if( restOfTheLine == "" ) restOfTheLine = ",";
			
			newBagLabels << restOfTheLine << "\n";
			
		}
		
		originalBagLabels.close();
		newBagLabels.close();
	}

	// the point of this is to make a visualization of which days/hours of the day data was recorded / at which sites
	static void visualizeRecordingCoverage()
	{
		vector<string> months;
		months.push_back("June");
		months.push_back("July");
		
		map<string, int> daysInMonth;
		daysInMonth["June"] = 30;
		daysInMonth["July"] = 13; // its actually 31, but the data only goes up to july 13th
		
		map<string, int> startAtDay;
		startAtDay["June"] = 4; // the data starts on june 4th
		startAtDay["July"] = 1;
		
		map<int, string> monthNumToName;
		monthNumToName[6] = "June";
		monthNumToName[7] = "July";
		
		map<string, int> monthNameToNum;
		monthNameToNum["June"] = 6;
		monthNameToNum["July"] = 7;
		
		/// parse the original file list to determine when recordings were collected /// 
		
		set<string> stamps; // stores a string encoding of each recorings site, month, day and hour
		
		vector<string> sourceFilenames;
		pFileUtils::getFilesInDirectory("src_wavs_unsplit", sourceFilenames);
		for(int i = 0; i < sourceFilenames.size(); ++i)
		{
			string filename = sourceFilenames[i];
			
			cout << filename << endl;
		
			vector<string> parts = pStringUtils::split(filename, "_.");
			string site = parts[0];
			string date = parts[1];
			string time = parts[2];
			
			char monthStr[2];
			char dayStr[2];
			char hourStr[2];
		
			monthStr[0] = date[4];
			monthStr[1] = date[5];
			dayStr[0]	= date[6];
			dayStr[1]	= date[7];
			hourStr[0]	= time[0];
			hourStr[1]	= time[1];
			
			int month	= pStringUtils::stringToInt(monthStr);
			int day		= pStringUtils::stringToInt(dayStr); 
			int hour	= pStringUtils::stringToInt(hourStr);

			string stamp = site +  "-" + pStringUtils::intToString(month) + "-" + pStringUtils::intToString(day) + "-" + pStringUtils::intToString(hour);
			stamps.insert(stamp);
			cout << stamp << endl; 
		}
		
		/// generate output visualization ///
			
		vector<string> siteCodes;
		siteCodes.push_back("PC1");
		siteCodes.push_back("PC2");
		siteCodes.push_back("PC4");
		siteCodes.push_back("PC5");
		siteCodes.push_back("PC7");
		siteCodes.push_back("PC8");
		siteCodes.push_back("PC10");
		siteCodes.push_back("PC11");
		siteCodes.push_back("PC13");
		siteCodes.push_back("PC15");
		siteCodes.push_back("PC16");
		siteCodes.push_back("PC17");
		siteCodes.push_back("PC18");
		

		pTextFile f("coverage.html", PFILE_WRITE);
		
		f << "<html><body style='font-family: Arial;'><table style='padding: 0px' cellspacing=0 cellpadding=0>";

		// make the header
		f << "<tr>";
		f << "<td>Site</td>";
		for(int m = 0; m < months.size(); ++m)
		{
			for(int d = startAtDay[months[m]]; d <= daysInMonth[months[m]]; ++d)
			{
				if( d == startAtDay[months[m]] ) // to reduce clutter, only draw the month for the first day
					f << "<td align='center' colspan=24 style='border-bottom: solid #333; border-bottom-width:1px;' width=24>"<< monthNameToNum[months[m]] << "/" << d <<"</td>";
				else 
					f << "<td align='center' colspan=24 style='border-bottom: solid #333; border-bottom-width:1px;' width=24>"<< d <<"</td>";
			}
		}
		f << "</tr>";
		
		// draw each row
		for(int s = 0; s < siteCodes.size(); ++s)
		{
			f << "<tr>";
			f << "<td style='padding-right:4px;'>" << siteCodes[s] << "</td>";
			
			for(int m = 0; m < months.size(); ++m)
			{
				for(int d = startAtDay[months[m]]; d <= daysInMonth[months[m]]; ++d)
				{
				
					bool allEmpty = true;
					bool allFull = true;
					
					for(int hr = 0; hr < 24; ++hr)
					{
						string stamp = siteCodes[s] + "-" + pStringUtils::intToString(monthNameToNum[months[m]]) + "-" + pStringUtils::intToString(d) + "-" + pStringUtils::intToString(hr);
						if( stamps.count(stamp) > 0 )
							allEmpty = false;
						
						if( stamps.count(stamp) == 0 )
							allFull = false;
					}
				
					if( allEmpty )
					{
						string rightPad = "";
						if( m == months.size()-1 && d == daysInMonth[months[m]] ) rightPad = "border-right: solid 1px #333;";
						f << "<td colspan=24 style='width: 24px; border-bottom: solid #333; border-bottom-width:1px;border-left: solid #333; border-left-width:1px; "+rightPad+"' width=24><div style='width: 24px'></div></td>";
					}
					//else if (allFull)
					//{
					//	string rightPad = "";
					//	if( m == months.size()-1 && d == daysInMonth[months[m]] ) rightPad = "border-right: solid 1px #333;";
					//	f << "<td colspan=24 style='background-color: #444444; width: 24px; border-bottom: solid #333; border-bottom-width:1px;border-left: solid #333; border-left-width:1px; "+rightPad+"' width=24><div style='width: 24px'></div></td>";
					//}
					else
					{
						for(int hr = 0; hr < 24; ++hr)
						{
							int colorR = 255;
							int colorG = colorR;
							int colorB = colorG;
							
							// construct a stamp for this time
							string stamp = siteCodes[s] + "-" + pStringUtils::intToString(monthNameToNum[months[m]]) + "-" + pStringUtils::intToString(d) + "-" + pStringUtils::intToString(hr);
							cout << stamp << endl;
							
							if( stamps.count(stamp) > 0 )
							{
								colorR = 68;
								colorG = 68;
								colorB = 68;
								
							}
							string leftPad = "";
							string rightPad = ""; // add a border on the right for the very last column
							if( hr == 0 ) leftPad = "border-left: solid #333; border-left-width:1px";
							if( m == months.size()-1 && d == daysInMonth[months[m]] && hr == 23 ) rightPad = "border-right: solid #333; border-right-width:1px"; 
							
							f << "<td style='background-color:rgb("<<colorR << "," << colorG << "," << colorB << "); width: 1px; border-bottom: solid #333; border-bottom-width:1px;"+leftPad+rightPad+"' width=1><div></div></td>";
						}
					}
				}
			}
			
			f << "</tr>";
		}

		f << "</table></body></html>";
		
		f.close();
	}

	// the files from HJA have a specific format which includes the time of day the recording started.
	// this function returns the hour of day (0-23) that the recording started. the format is
	// <sitename>_<year><month><day>_<hour><minute><second>_<other info>.<extension>
	static int getHourFromFilename(string filename)
	{
		vector<string> parts = pStringUtils::split(filename, "_.");
		assert(parts.size() > 4);
		
		string timePart = parts[2];
		assert(timePart.size() == 6);
		
		string hourString = timePart.substr(0,2);
		
		cout << "hrString = " << hourString << endl;
		
		int hr = pStringUtils::stringToInt(hourString);
		assert( hr >= 0 && hr <= 23 );
		return hr;
	}
	
	static string getSiteFromFilename(string filename)
	{
		vector<string> parts = pStringUtils::split(filename, "_.");
		assert(parts.size() > 4);
		return parts[0];
	}
	
	static pair<int,int> getDateFromFilename(string filename)
	{
		vector<string> parts = pStringUtils::split(filename, "_.");
		assert(parts.size() > 4);
		
		string datePart = parts[1];
		assert(datePart.size() == 8);
		
		string monthString = datePart.substr(4,2);
		string dayString = datePart.substr(6,2);
	
		cout << "monthStr=" << monthString << endl;
		cout << "dayString=" << dayString << endl;
		
		int month = pStringUtils::stringToInt(monthString);
		int day = pStringUtils::stringToInt(dayString);
		assert( month >= 1 && month <= 12 );
		assert( day >= 1 && day <= 31 );
		return pair<int,int>(month, day);
	}
	
	static int getYearFromFilename(string filename)
	{
		vector<string> parts = pStringUtils::split(filename, "_.");
		assert(parts.size() > 4);
		
		string datePart = parts[1];
		assert(datePart.size() == 8);
		
		string yearString = datePart.substr(0,4);
		cout << "yearString=" << yearString << endl;
		return pStringUtils::stringToInt(yearString);
	}
	
	// this analysis counts the number of segments in each hour of the day, aggregated over the whole dataset.
	// the point is to answer the question: when are birds (or other things that are detected by the system) making sounds?
	static void analyzeSegmentsPerHour()
	{
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
		
		cout << "counting segments" << endl;
		vector<int> segmentsPerHour(24, 0);
		for(int i = 0; i < metaData.size(); ++i)
		{
			cout << "filename=" << metaData[i].filename_ << endl;
			int hourOfDay = getHourFromFilename(metaData[i].filename_);
			segmentsPerHour[hourOfDay] += metaData[i].segments_.size();
		}
		
		cout << "*****" << endl;
		
		for(int i = 0; i < 24; ++i)
			cout << i << "\t" << segmentsPerHour[i] << endl;
	}
	
	// counts the number of segments in each day. results are returned relative to the first day in the dataset (june 4).
	static void analyzeSegmentsPerDay(int numClusters)
	{	
		// count the number of total days in the dataset
		
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
		pMetaData_Parser::parseInstanceClusterAssignment("segment2cluster.txt", metaData);
		
		int startMonth = 6;
		int startDay = 4;
		int endMonth = 7;
		int endDay = 13;
		
		map<int, int> daysInMonth;
		daysInMonth[6] = 30;
		daysInMonth[7] = 31;

		int totalDays = 0;
		for(int month = startMonth; month <= endMonth; ++month)
		{
			if( month == startMonth ) totalDays += daysInMonth[month]-startDay+1;
			else if( month == endMonth ) totalDays += endDay;
			else totalDays += daysInMonth[month];
		}
		
		cout << "total days = " << totalDays << endl;
		
		vector<string> siteCodes;
		siteCodes.push_back("PC1");
		siteCodes.push_back("PC2");
		siteCodes.push_back("PC4");
		siteCodes.push_back("PC5");
		siteCodes.push_back("PC7");
		siteCodes.push_back("PC8");
		siteCodes.push_back("PC10");
		siteCodes.push_back("PC11");
		siteCodes.push_back("PC13");
		siteCodes.push_back("PC15");
		siteCodes.push_back("PC16");
		siteCodes.push_back("PC17");
		siteCodes.push_back("PC18");
		
		map<string, vector<int> > segmentsPerDay;
		for(int s = 0; s < siteCodes.size(); ++s)
			segmentsPerDay[siteCodes[s]] = vector<int>(totalDays+1, 0);
			
		map<string, vector< vector<int> > > segmentsPerDayByCluster;
		for(int s = 0; s < siteCodes.size(); ++s)
			segmentsPerDayByCluster[siteCodes[s]] = vector< vector<int> >(numClusters, vector<int>(totalDays+1, 0));
		
		for(int i = 0; i < metaData.size(); ++i)
		{
			cout << "filename=" << metaData[i].filename_ << endl;
			
			string site = getSiteFromFilename(metaData[i].filename_);
			
			pair<int,int> date = getDateFromFilename(metaData[i].filename_);
			int month = date.first;
			int day = date.second;
			
			int dateOffset;
			
			// TODO: this math is specific to the range of dates in 2009 (i.e. june-july).
			if( month == 6 ) dateOffset = day - startDay;
			else if ( month == 7 ) dateOffset = daysInMonth[6] - startDay + day; // is this right?
			else assert(false);
			
			cout << "dateOffset=" << dateOffset << endl;
			
			assert(dateOffset >= 0 && dateOffset <= totalDays);
			
			segmentsPerDay[site][dateOffset] += metaData[i].segments_.size();
			
			for(int j = 0; j < metaData[i].segments_.size(); ++j)
				segmentsPerDayByCluster[site][metaData[i].segments_[j].clusterAssignment_][dateOffset] += 1;
		}
		
		cout << "generating result files" << endl;
		
		for(int s = 0; s < siteCodes.size(); ++s)
		{
			cout << "making segments per day for site " << siteCodes[s] << endl;
		
			pTextFile f("segments_by_site/" + siteCodes[s] + ".txt", PFILE_WRITE);
			f << "day,num_segments" << "\n";
			for(int i = 0; i < totalDays; ++i)
				f << i << "\t" << segmentsPerDay[siteCodes[s]][i] << "\n";
			f.close();
			
			cout << "analyzing results at site " << siteCodes[s] << " by cluster" << endl;
			
			for(int cluster = 0; cluster < numClusters; ++cluster)
			{
				pTextFile f("segments_by_site_and_cluster/" + siteCodes[s] + "_cluster_" + pStringUtils::intToString(cluster) + ".txt", PFILE_WRITE);
				f << "day\t" << siteCodes[s] << "\n";
				for(int i = 0; i < totalDays; ++i)
					f << i << "\t" << segmentsPerDayByCluster[siteCodes[s]][cluster][i] << "\n";
				  
				f.close();
			}
		}

	}
	
	// counts the number of segments belonging to each cluster, within each hour of day
	// the point is to answer the question: when is a particular sound being made? (at a low human effort cost)
	static void analyzeSegmentsPerHourByCluster(int numClusters)
	{
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
		pMetaData_Parser::parseInstanceClusterAssignment("segment2cluster.txt", metaData);
		
		vector< vector<int> > segmentsPerHour(numClusters, vector<int>(24, 0));
		for(int i = 0; i < metaData.size(); ++i)
		{
			cout << "filename=" << metaData[i].filename_ << endl;
			int hourOfDay = getHourFromFilename(metaData[i].filename_);
			for(int j = 0; j < metaData[i].segments_.size(); ++j)
				segmentsPerHour[metaData[i].segments_[j].clusterAssignment_][hourOfDay] += 1;
		}
		
		for(int cluster = 0; cluster < numClusters; ++cluster)
		{
			pTextFile f("segments_per_hour_by_cluster/" + pStringUtils::intToString(cluster) + ".txt", PFILE_WRITE);
			for(int hr = 0; hr < 24; ++hr)
				f << hr << "\t" << segmentsPerHour[cluster][hr] << "\n"; 
			f.close();
		}
	}
	
	// In "Sampling environmental acoustic recordings to determine bird species richness", Wimmer et al. concluded
	// that the sampling method that produced the most species for a fixed amount of effort was to randomly select recordings from 3
	// hours after dawn. In our experiments, we consider this approach as a baseline. This function selects a random sample of such
	// recordings, and sets up a folder with the necessary data to have our experts label those recordings.
	// Note: this is not exactly "dawn", instead we just go from 5am to 8am.
	static void speciesDiscovery_randomDawn(int numRecordingsToSelect)
	{
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		
		vector<int> candidateRecordings;
		
		for(int i = 0; i < metaData.size(); ++i)
		{
			cout << "filename=" << metaData[i].filename_ << endl;
			int hourOfDay = getHourFromFilename(metaData[i].filename_);
			
			if( hourOfDay >= 5 && hourOfDay <= 8 )
			{	
				cout << "dawn" << endl;
				candidateRecordings.push_back(i);
			}
		}
		
		cout << "# dawn recordings = " << candidateRecordings.size() << endl;
		
		// select k from the dawn recordings at random
		pRandom::shuffle(candidateRecordings);
		vector<int> selectedRecordings;
		for(int i = 0; i < numRecordingsToSelect; ++i)
			selectedRecordings.push_back(candidateRecordings[i]);
			
		// copy the source files for the selected recordings to a sub-folder so they can be delivered to the experts,
		// and make a template .csv file to be filled out in Excel with labels
		pTextFile labelTemplate("species_discovery/dawn/labels_dawn.csv", PFILE_WRITE);
		labelTemplate << "filename" << "\n";
		for(int i = 0; i < selectedRecordings.size(); ++i)
		{
			pAudioFile_MetaData& rec = metaData[selectedRecordings[i]];
			
			cout << rec.filename_ << endl;
			string prefix = pStringUtils::firstPartOfSplit(rec.filename_, ".");
			
			labelTemplate  << prefix << "\n";
			
			system(string("cp spectrograms/" + rec.filename_ + " species_discovery/dawn/spectrograms/").c_str());
			system(string("cp src_wavs/" + prefix + ".wav species_discovery/dawn/src_wavs/").c_str());
		}
		
		labelTemplate.close();
	}
	
	// run our proposed method cluster-coverage with MIFF tie breaking
	static void speciesDiscovery_unifiedCCMIFF(int numRecordingsToSelect, int numClusters, int miffK, string destFolder)
	{
		cout << "running unified CCMIFF, " << destFolder << endl;
		
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
		
		const int featureDim = pMetaData_Parser::getSegmentFeatureDim(metaData);
		cout << "feature dim = " << featureDim << endl;
		cout << "rescaling features" << endl;
		pMetaData_Rescale::rescaleTo01(metaData, featureDim);
		
		vector<pMIMLExample> mimlExs;
		for(int i = 0; i < metaData.size(); ++i)
		{
			vector< vector<float> > segmentFeatures;
			for(int j = 0; j < metaData[i].segments_.size(); ++j)
				segmentFeatures.push_back(metaData[i].segments_[j].featureVector_);
			
			vector<int> emptyLabelSet;
			mimlExs.push_back(pMIMLExample(segmentFeatures, emptyLabelSet));
		}
		
		vector<int> selectedRecordings = pMIMLClassDiscoveryExperiment::discoveryMethod_unified_clusterCoverage_MIFF(mimlExs, numRecordingsToSelect, numClusters, miffK);
		
		// copy the source files for the selected recordings to a sub-folder so they can be delivered to the experts,
		// and make a template .csv file to be filled out in Excel with labels
		pTextFile labelTemplate("species_discovery/"+destFolder+"/labels_"+destFolder+".csv", PFILE_WRITE);
		labelTemplate << "order" << "," << "filename" << "\n";
		for(int i = 0; i < selectedRecordings.size(); ++i)
		{
			pAudioFile_MetaData& rec = metaData[selectedRecordings[i]];
			
			cout << rec.filename_ << endl;
			string prefix = pStringUtils::firstPartOfSplit(rec.filename_, ".");
			
			labelTemplate << i << "," << prefix << "\n";
			
			system(string("cp spectrograms/" + rec.filename_ + " species_discovery/"+destFolder+"/spectrograms/").c_str());
			system(string("cp filtered_spectrograms/" + rec.filename_ + " species_discovery/"+destFolder+"/filtered_spectrograms/").c_str());
			system(string("cp segmentation_outlines/" + rec.filename_ + " species_discovery/"+destFolder+"/segmentation_outlines/").c_str());
			system(string("cp src_wavs/" + prefix + ".wav species_discovery/"+destFolder+"/src_wavs/").c_str());
		}
		
		labelTemplate.close();
	}
	
	// run our proposed method cluster-coverage with MIFF tie breaking (with rain filter)
	static void speciesDiscovery_rainfilter_unifiedCCMIFF(int numRecordingsToSelect, int numClusters, int miffK, string destFolder, float thresholdRainProbability)
	{
		cout << "running rain-filtered unified CCMIFF, " << destFolder << endl;
		
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
		pMetaData_Parser::parseRainProbabilities(metaData);
		
		const int featureDim = pMetaData_Parser::getSegmentFeatureDim(metaData);
		cout << "feature dim = " << featureDim << endl;
		cout << "rescaling features" << endl;
		pMetaData_Rescale::rescaleTo01(metaData, featureDim);
		
		map<int, int> mapMimlExIndexToBagID; // since we are skipping some recordings, the indices of bags in the miml dataset and recordings will no longer match up. this map allows us to reverse this issue
		 
		vector<pMIMLExample> mimlExs;
		for(int i = 0; i < metaData.size(); ++i)
		{
			vector< vector<float> > segmentFeatures;
			for(int j = 0; j < metaData[i].segments_.size(); ++j)
				segmentFeatures.push_back(metaData[i].segments_[j].featureVector_);
			
			vector<int> emptyLabelSet;
			
			if( metaData[i].rainProbability_ < thresholdRainProbability )
			{
				mimlExs.push_back(pMIMLExample(segmentFeatures, emptyLabelSet));
				mapMimlExIndexToBagID[mimlExs.size()-1] = metaData[i].bag_id_;
			}
		}
		
		cout << "bags left after rain filter = " << mimlExs.size() << endl;
		
		vector<int> selectedBags = pMIMLClassDiscoveryExperiment::discoveryMethod_unified_clusterCoverage_MIFF(mimlExs, numRecordingsToSelect, numClusters, miffK);
		vector<int> selectedRecordings;
		for(int i = 0; i < selectedBags.size(); ++i)
			selectedRecordings.push_back(mapMimlExIndexToBagID[selectedBags[i]]);
		
		// copy the source files for the selected recordings to a sub-folder so they can be delivered to the experts,
		// and make a template .csv file to be filled out in Excel with labels
		pTextFile labelTemplate("species_discovery/"+destFolder+"/labels_"+destFolder+".csv", PFILE_WRITE);
		labelTemplate << "order" << "," << "filename" << "\n";
		for(int i = 0; i < selectedRecordings.size(); ++i)
		{
			pAudioFile_MetaData& rec = metaData[selectedRecordings[i]];
			
			cout << rec.filename_ << endl;
			string prefix = pStringUtils::firstPartOfSplit(rec.filename_, ".");
			
			labelTemplate << i << "," << prefix << "\n";
			
			system(string("cp spectrograms/" + rec.filename_ + " species_discovery/"+destFolder+"/spectrograms/").c_str());
			system(string("cp filtered_spectrograms/" + rec.filename_ + " species_discovery/"+destFolder+"/filtered_spectrograms/").c_str());
			system(string("cp segmentation_outlines/" + rec.filename_ + " species_discovery/"+destFolder+"/segmentation_outlines/").c_str());
			system(string("cp src_wavs/" + prefix + ".wav species_discovery/"+destFolder+"/src_wavs/").c_str());
		}
		
		labelTemplate.close();
	}
	
	
	// same as above, but just runs plain MIFF (no cluster coverage)
	static void speciesDiscovery_MIFF(int numRecordingsToSelect, int miffK, string destFolder)
	{
		cout << "running MIFF, " << destFolder << endl;
	
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
		
		const int featureDim = pMetaData_Parser::getSegmentFeatureDim(metaData);
		cout << "feature dim = " << featureDim << endl;
		cout << "rescaling features" << endl;
		pMetaData_Rescale::rescaleTo01(metaData, featureDim);
		
		vector<pMIMLExample> mimlExs;
		for(int i = 0; i < metaData.size(); ++i)
		{
			vector< vector<float> > segmentFeatures;
			for(int j = 0; j < metaData[i].segments_.size(); ++j)
				segmentFeatures.push_back(metaData[i].segments_[j].featureVector_);
			
			vector<int> emptyLabelSet;
			mimlExs.push_back(pMIMLExample(segmentFeatures, emptyLabelSet));
		}
		
		vector<int> selectedRecordings = pMIMLClassDiscoveryExperiment::discoveryMethod_MIFF(mimlExs, numRecordingsToSelect, miffK);
	
		// copy the source files for the selected recordings to a sub-folder so they can be delivered to the experts,
		// and make a template .csv file to be filled out in Excel with labels
		
		pTextFile labelTemplate("species_discovery/"+destFolder+"/labels.csv", PFILE_WRITE);
		labelTemplate << string("order,filename\n");
		
		for(int i = 0; i < selectedRecordings.size(); ++i)
		{
			pAudioFile_MetaData& rec = metaData[selectedRecordings[i]];
			
			cout << rec.filename_ << endl;
			string prefix = pStringUtils::firstPartOfSplit(rec.filename_, ".");
			
			labelTemplate << i << "," << prefix << "\n";
			
			system(string("cp spectrograms/" + rec.filename_ + " species_discovery/"+destFolder+"/spectrograms/").c_str());
			system(string("cp filtered_spectrograms/" + rec.filename_ + " species_discovery/"+destFolder+"/filtered_spectrograms/").c_str());
			system(string("cp segmentation_outlines/" + rec.filename_ + " species_discovery/"+destFolder+"/segmentation_outlines/").c_str());
			system(string("cp src_wavs/" + prefix + ".wav species_discovery/"+destFolder+"/src_wavs/").c_str());
		}
		
		labelTemplate.close();
	}
	
	
	// runs plain MIFF on rain-filtered subset
	static void speciesDiscovery_rainfilter_MIFF(int numRecordingsToSelect, int miffK, string destFolder, float thresholdRainProbability)
	{
		cout << "running MIFF, " << destFolder << endl;
	
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
		pMetaData_Parser::parseRainProbabilities(metaData);
		
		const int featureDim = pMetaData_Parser::getSegmentFeatureDim(metaData);
		cout << "feature dim = " << featureDim << endl;
		cout << "rescaling features" << endl;
		pMetaData_Rescale::rescaleTo01(metaData, featureDim);
		
		map<int, int> mapMimlExIndexToBagID; // since we are skipping some recordings, the indices of bags in the miml dataset and recordings will no longer match up. this map allows us to reverse this issue
		 
		vector<pMIMLExample> mimlExs;
		for(int i = 0; i < metaData.size(); ++i)
		{
			vector< vector<float> > segmentFeatures;
			for(int j = 0; j < metaData[i].segments_.size(); ++j)
				segmentFeatures.push_back(metaData[i].segments_[j].featureVector_);
			
			vector<int> emptyLabelSet;
			
			if( metaData[i].rainProbability_ < thresholdRainProbability )
			{
				mimlExs.push_back(pMIMLExample(segmentFeatures, emptyLabelSet));
				mapMimlExIndexToBagID[mimlExs.size()-1] = metaData[i].bag_id_;
			}
		}
		
		cout << "bags left after rain filter = " << mimlExs.size() << endl;
		
		vector<int> selectedBags = pMIMLClassDiscoveryExperiment::discoveryMethod_MIFF(mimlExs, numRecordingsToSelect, miffK);
		vector<int> selectedRecordings;
		for(int i = 0; i < selectedBags.size(); ++i)
			selectedRecordings.push_back(mapMimlExIndexToBagID[selectedBags[i]]);
			
		// copy the source files for the selected recordings to a sub-folder so they can be delivered to the experts,
		// and make a template .csv file to be filled out in Excel with labels
		
		pTextFile labelTemplate("species_discovery/"+destFolder+"/labels.csv", PFILE_WRITE);
		labelTemplate << string("order,filename\n");
		
		for(int i = 0; i < selectedRecordings.size(); ++i)
		{
			pAudioFile_MetaData& rec = metaData[selectedRecordings[i]];
			
			cout << rec.filename_ << endl;
			string prefix = pStringUtils::firstPartOfSplit(rec.filename_, ".");
			
			labelTemplate << i << "," << prefix << "\n";
			
			system(string("cp spectrograms/" + rec.filename_ + " species_discovery/"+destFolder+"/spectrograms/").c_str());
			system(string("cp filtered_spectrograms/" + rec.filename_ + " species_discovery/"+destFolder+"/filtered_spectrograms/").c_str());
			system(string("cp segmentation_outlines/" + rec.filename_ + " species_discovery/"+destFolder+"/segmentation_outlines/").c_str());
			system(string("cp src_wavs/" + prefix + ".wav species_discovery/"+destFolder+"/src_wavs/").c_str());
		}
		
		labelTemplate.close();
	}
	
	// select exapmles closest to cluster centers
	static void speciesDiscovery_clusterCenters(int numRecordingsToSelect, string destFolder)
	{
		cout << "running cluster center discovery method" << endl;
		
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
		
		const int featureDim = pMetaData_Parser::getSegmentFeatureDim(metaData);
		cout << "feature dim = " << featureDim << endl;
		cout << "rescaling features" << endl;
		pMetaData_Rescale::rescaleTo01(metaData, featureDim);
		
		vector<pMIMLExample> mimlExs;
		for(int i = 0; i < metaData.size(); ++i)
		{
			vector< vector<float> > segmentFeatures;
			for(int j = 0; j < metaData[i].segments_.size(); ++j)
				segmentFeatures.push_back(metaData[i].segments_[j].featureVector_);
			
			vector<int> emptyLabelSet;
			mimlExs.push_back(pMIMLExample(segmentFeatures, emptyLabelSet));
		}
		
		vector<int> selectedRecordings = pMIMLClassDiscoveryExperiment::discoveryMethod_clusterCenters(mimlExs, numRecordingsToSelect);
		
		// copy the source files for the selected recordings to a sub-folder so they can be delivered to the experts,
		// and make a template .csv file to be filled out in Excel with labels
		
		pTextFile labelTemplate("species_discovery/"+destFolder+"/labels.csv", PFILE_WRITE);
		labelTemplate << string("order,filename\n");
		
		for(int i = 0; i < selectedRecordings.size(); ++i)
		{
			pAudioFile_MetaData& rec = metaData[selectedRecordings[i]];
			
			cout << rec.filename_ << endl;
			string prefix = pStringUtils::firstPartOfSplit(rec.filename_, ".");
			
			labelTemplate << i << "," << prefix << "\n";
			
			system(string("cp spectrograms/" + rec.filename_ + " species_discovery/"+destFolder+"/spectrograms/").c_str());
			system(string("cp filtered_spectrograms/" + rec.filename_ + " species_discovery/"+destFolder+"/filtered_spectrograms/").c_str());
			system(string("cp segmentation_outlines/" + rec.filename_ + " species_discovery/"+destFolder+"/segmentation_outlines/").c_str());
			system(string("cp src_wavs/" + prefix + ".wav species_discovery/"+destFolder+"/src_wavs/").c_str());
		}
		
		labelTemplate.close();
	}
	
	// select exapmles closest to cluster centers (with rain filter)
	static void speciesDiscovery_rainfilter_clusterCenters(int numRecordingsToSelect, string destFolder, float thresholdRainProbability)
	{
		cout << "running cluster center discovery method" << endl;
		
		vector<pAudioFile_MetaData> metaData;
		pMetaData_Parser::parseBags("bag_id2filename.txt", metaData);
		pMetaData_Parser::parseSegmentFeatures("segment_features.txt", metaData);
		
		const int featureDim = pMetaData_Parser::getSegmentFeatureDim(metaData);
		cout << "feature dim = " << featureDim << endl;
		cout << "rescaling features" << endl;
		pMetaData_Rescale::rescaleTo01(metaData, featureDim);
		
		map<int, int> mapMimlExIndexToBagID; // since we are skipping some recordings, the indices of bags in the miml dataset and recordings will no longer match up. this map allows us to reverse this issue
		 
		vector<pMIMLExample> mimlExs;
		for(int i = 0; i < metaData.size(); ++i)
		{
			vector< vector<float> > segmentFeatures;
			for(int j = 0; j < metaData[i].segments_.size(); ++j)
				segmentFeatures.push_back(metaData[i].segments_[j].featureVector_);
			
			vector<int> emptyLabelSet;
			
			if( metaData[i].rainProbability_ < thresholdRainProbability )
			{
				mimlExs.push_back(pMIMLExample(segmentFeatures, emptyLabelSet));
				mapMimlExIndexToBagID[mimlExs.size()-1] = metaData[i].bag_id_;
			}
		}
		
		cout << "bags left after rain filter = " << mimlExs.size() << endl;
		
		vector<int> selectedBags = pMIMLClassDiscoveryExperiment::discoveryMethod_clusterCenters(mimlExs, numRecordingsToSelect);
		vector<int> selectedRecordings;
		for(int i = 0; i < selectedBags.size(); ++i)
			selectedRecordings.push_back(mapMimlExIndexToBagID[selectedBags[i]]);
			
		// copy the source files for the selected recordings to a sub-folder so they can be delivered to the experts,
		// and make a template .csv file to be filled out in Excel with labels
		
		pTextFile labelTemplate("species_discovery/"+destFolder+"/labels.csv", PFILE_WRITE);
		labelTemplate << string("order,filename\n");
		
		for(int i = 0; i < selectedRecordings.size(); ++i)
		{
			pAudioFile_MetaData& rec = metaData[selectedRecordings[i]];
			
			cout << rec.filename_ << endl;
			string prefix = pStringUtils::firstPartOfSplit(rec.filename_, ".");
			
			labelTemplate << i << "," << prefix << "\n";
			
			system(string("cp spectrograms/" + rec.filename_ + " species_discovery/"+destFolder+"/spectrograms/").c_str());
			system(string("cp filtered_spectrograms/" + rec.filename_ + " species_discovery/"+destFolder+"/filtered_spectrograms/").c_str());
			system(string("cp segmentation_outlines/" + rec.filename_ + " species_discovery/"+destFolder+"/segmentation_outlines/").c_str());
			system(string("cp src_wavs/" + prefix + ".wav species_discovery/"+destFolder+"/src_wavs/").c_str());
		}
		
		labelTemplate.close();
	}
	
	
	// I messed up and sorted the filenames for the random dawn sample before giving them to the expert to label. This function puts the
	// results we got back in their original (random) order
	static void discoveryAnalysis_reorderDawn()
	{
		map<string, string> labeledLines;
		pTextFile labeledDawn("dawn/labels_dawn.csv", PFILE_READ);
		labeledDawn.readLine(); // skip the header
		while(!labeledDawn.eof())
		{
			string line = labeledDawn.readLine();
			vector<string> parts = pStringUtils::split(line, ",");
			//cout << parts[0] << "--" << line << endl;
			labeledLines[parts[0]] = line;
		}
		labeledDawn.close();
		
		pTextFile orderFile("dawn/labels_dawn_original_order.csv", PFILE_READ);
		orderFile.readLine(); // skip header
		
		pTextFile fixedFile("dawn/labels_dawn_reorder.csv", PFILE_WRITE);
		fixedFile << "filename,species,notes" << "\n";
		while(!orderFile.eof())
		{
			string line = orderFile.readLine();
			string lookupLine = labeledLines[line];
			cout << line << "--" << lookupLine << endl;
			fixedFile << lookupLine << "\n";
		}
		orderFile.close();
		fixedFile.close();
	}
	
	static void discoveryAnalysis_countResults(string labelFilename)
	{
		pTextFile labelFile(labelFilename, PFILE_READ);
		labelFile.readLine(); // skip header
		int lineNum = 1;
		
		set<string> codesToSkip;
		
		codesToSkip.insert("UNWA"); // UNWA = uknown warbler		
		for(int i = 1; i <= 37; ++i)
			codesToSkip.insert("UNKN" + pStringUtils::intToString(i));

		codesToSkip.insert("UNKN"); // no # suffix
		codesToSkip.insert("UNTH");
		
		set<string> nonSpeciesCodes;
		nonSpeciesCodes.insert("BEEP");
		nonSpeciesCodes.insert("MICBUMP");
		nonSpeciesCodes.insert("SMGLITCH");
		nonSpeciesCodes.insert("MOTOR");
		nonSpeciesCodes.insert("RAIN");
		nonSpeciesCodes.insert("WIND");
		nonSpeciesCodes.insert("ARPL");
		nonSpeciesCodes.insert("HAMMER");
		nonSpeciesCodes.insert("BREAK");
		
		
		map<string, int> classCodeCount;
		set<string> knownSpecies;
		
		vector<int> classDiscoveryCurve;
		vector<int> speciesDiscoveryCurve;
		
		classDiscoveryCurve.push_back(0);
		speciesDiscoveryCurve.push_back(0);
		
		while(!labelFile.eof())
		{
			++lineNum;
			
			string line = labelFile.readLine();
			
			cout << lineNum << ":" << line << endl;
			
			vector<string> parts = pStringUtils::splitNonEmpty(line, ",\"");
			
			vector<string> speciesCodes;
			for(int i = 1; i < parts.size(); ++i) speciesCodes.push_back(parts[i]);
			
			for(int i = 0; i < speciesCodes.size(); ++i)
			{
				cout << speciesCodes[i] << "|";
				
				if( codesToSkip.count(speciesCodes[i]) != 0 ) continue; // skip this code
				
				if( classCodeCount.count(speciesCodes[i]) == 0 )
					classCodeCount[speciesCodes[i]] = 1;
				else
					classCodeCount[speciesCodes[i]] = classCodeCount[speciesCodes[i]] + 1;
					
				if( nonSpeciesCodes.count(speciesCodes[i]) == 0 ) // if this is a valid species (not some other noise)
					knownSpecies.insert(speciesCodes[i]);
			}
			
			classDiscoveryCurve.push_back(classCodeCount.size());
			speciesDiscoveryCurve.push_back(knownSpecies.size());
			
			cout << endl << endl;
		}
		labelFile.close();
		
		cout << "*** species code counts ***" << endl;
		for(map<string, int>::iterator it = classCodeCount.begin(); it != classCodeCount.end(); ++it)
		{
			cout << (*it).first << "\t" << (*it).second << endl;
		}
		
		cout << "*** class discovery curve ***" << endl;
		cout << "#labeled\t#classes" << endl;
		for(int i = 0; i < classDiscoveryCurve.size(); ++i)
			cout << i << "\t" << classDiscoveryCurve[i] << endl;
			
		cout << "*** species discovery curve ***" << endl;
		cout << "#labeled\t#species" << endl;
		for(int i = 0; i < speciesDiscoveryCurve.size(); ++i)
			cout << i << "\t" << speciesDiscoveryCurve[i] << endl;
	}
	
	// for each discovery method, we output a histogram of (class,#of recordings)
	// however, different methods get a different set of classes. this function joins 
	// all of these histograms together
	static void discoveryAnalysis_joinSelectionHistograms()
	{	
		vector<string> histoFilenames; 
		vector<string> methods;
		set<string> joinedClasses;
		map<string, int> joinedCounts;
		
		pFileUtils::getFilesInDirectory("selection_histograms", histoFilenames);
		for(int i = 0; i < histoFilenames.size(); ++i)
		{
			methods.push_back(pStringUtils::firstPartOfSplit(histoFilenames[i], "."));
			cout << "** " << methods[i] << endl;
			
			pTextFile histoFile("selection_histograms/" + histoFilenames[i], PFILE_READ);
			while(!histoFile.eof())
			{
				vector<string> parts = pStringUtils::split(histoFile.readLine(), "\t");
				string className = parts[0];
				int count = pStringUtils::stringToInt(parts[1]);
				joinedClasses.insert(className);
				joinedCounts[methods[i] + "-" + className] = count;
			}
			histoFile.close();
		}
		
		cout << "***** joined classes *****" << endl;
		cout << "class";
		for(int i = 0; i < methods.size(); ++i)
			cout << "\t" << methods[i];
		cout << endl;
		
		for(set<string>::iterator it = joinedClasses.begin(); it != joinedClasses.end(); ++it)
		{
			cout << *it;
			for(int i = 0; i < methods.size(); ++i)
				cout << "\t" << joinedCounts[methods[i] + "-" + (*it)];
			cout << endl;
		}
	}
	
	// Excel isn't really up to the task of making a good graph from the selection histograms.
	// instead, I will just generate some HTML code for a graph...
	static  void discoveryAnalysis_makeSelectionHistoGraph()
	{
		pTextFile f("joined_selection_histo.txt", PFILE_READ);
		
		string header = f.readLine();
		vector<string> headerParts = pStringUtils::splitNonEmpty(header, "\t");
		vector<string> classCodes;
		map<string, int> counts;
		
		while(!f.eof())
		{
			string line = f.readLine();
			vector<string> lineParts = pStringUtils::splitNonEmpty(line, "\t");
			string classCode = lineParts[0];
			classCodes.push_back(classCode);
			
			assert(lineParts.size() == headerParts.size());
			
			for(int i = 1; i < lineParts.size(); ++i)
			{
				string key = headerParts[i] + "-" + classCode;
				counts[key] = pStringUtils::stringToInt(lineParts[i]);
				
				//cout << "key=" << key << " count = " << counts[key] << endl;
			}
		}
	
		f.close();
		
		cout << "<!DOCTYPE html><html><head></head><body style='font-family: Arial;'>";
		
		/*
		
		for(int i = 1; i < headerParts.size(); ++i)
		{
			cout << "<h2>" << headerParts[i] << "</h2>" << endl;
			cout << "<table cellpadding=0 cellspacing=1>";
			for(int j = 0; j < classCodes.size(); ++j)
			{
				string key = headerParts[i] + "-" + classCodes[j];
				int count = counts[key];
				
				cout << "<tr>";
				cout << "<td>"+classCodes[j]+"</td>";
				cout << "<td>" << count << "</td>";
				cout << "<td><span style='display: inline-block; background-color: #0000aa; width: "<< 5 * count <<"px'>&nbsp;</span></td>";
				cout << "</tr>" << endl;
			}
			cout << "</table>";
		}
		
		*/
		
		// make one giant mega table
		
		cout << "<div style='height: 230px'></div>";
		cout << "<table cellpadding=0 cellspacing=0 border=0>";
		
		cout << "<tr>";
		for(int i = 0; i < headerParts.size(); ++i)
			cout << "<td style='width: 91px;'><div style='width: 91px;'><div style='-webkit-transform: translate(-32px,-40px) rotate(-90deg);'>"<<headerParts[i]<<"</div></div></td>";
		cout << "</tr>";
		
		cout << "<tr><td colspan=12 style='border-bottom: solid 1px black; border-top: solid 1px black;' align=center>Bird-Species</td></tr>";
			
		for(int i = 0; i < classCodes.size(); ++i)
		{
			if( classCodes[i] == "RAIN" )
				cout << "<tr><td colspan=12 style='border-bottom: solid 1px black; border-top: solid 1px black;' align=center>Other Sounds</td></tr>";
			
			cout << "<tr>";
			cout << "<td style='border-bottom: solid 1px #cccccc;'>"<< classCodes[i] <<"</td>";
			for(int j = 1; j < headerParts.size(); ++j)
			{
				string key = headerParts[j] + "-" + classCodes[i];
				int count = counts[key];
				
				cout << "<td style='border-left: solid 1px #cccccc; border-bottom: solid 1px #cccccc'>" << "<span style='display: inline-block; width: 24px; background-color: #eeeeee'>" << count << "</span>" << "<span style='display: inline-block; background-color: #444444; width: "<< count <<"px'>&nbsp;</span></td>";
			}
			cout << "</tr>";
		}
		
		cout << "</table>";
		
		cout << "</body></html>";
	}
};
