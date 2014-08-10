// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#pragma once

#include <vector>
#include <string>
#include <math.h>
#include "wav_in.h"
#include "wav_out.h"
#include <iostream>
#include <assert.h>

using namespace std;

class pWavData
	{
	public:
		int sampleRate_;
		unsigned int bitsPerSample_;
		unsigned int channels_;
		vector<float> samples_;
		
		pWavData() {}
		pWavData(int sampleRate, unsigned int bitsPerSample, unsigned int channels) :
		sampleRate_(sampleRate), bitsPerSample_(bitsPerSample), channels_(channels) {}
		
		pWavData(string filename)
		{
			WAV_IN infile(filename.c_str());
			sampleRate_ = infile.get_sample_rate_hz();
			bitsPerSample_ = infile.get_bits_per_sample();
			channels_ = infile.get_num_channels();
			
			if( channels_ == 2  ) // if 2 channels, we need to convert to mono
			{
				cout << "CONVERTING STEREO TO MONO!!!" << endl;
				
				vector<float> samplesLR;
				while(infile.more_data_available()) 
					samplesLR.push_back(infile.read_current_input() / 32768.0); // 32768 is the maximum magnitude for a 16bit signed integer
				
				assert(samplesLR.size() % 2 == 0); // should be an even number of samples here
				
				for(int i = 0; i < samplesLR.size(); i +=2)
				{
					double avg = .5 * (samplesLR[i] + samplesLR[i+1]);
					samples_.push_back(avg);
				}
				
				channels_ = 1;
			}
			else
			{
				// read data from input file into memory
				while(infile.more_data_available()) 
				{
					double sample = infile.read_current_input();
					samples_.push_back(sample / 32768.0); // this is the maximum magnitude for a 16bit signed integer
				}
				
			}
			
			// normalize 16-bit samples to -1, 1
			//for(int i = 0; i < samples_.size(); ++i)
			//	samples_[i] /= ;  
		}
		
		// split this wav into chunks of duration (seconds)
		void splitIntoChunks(float duration, vector<pWavData>& outChunks)
		{
			int samplesPerChunk = sampleRate_ * duration;
			int totalChunks = floor((float)samples_.size() / samplesPerChunk);
			
			int currSample = 0;
			for(int i = 0; i < totalChunks; ++i)
			{
				pWavData chunk(sampleRate_, bitsPerSample_, channels_);
				for(int j = 0; j < samplesPerChunk; ++j)
				{
					chunk.samples_.push_back(samples_[currSample]);
					++currSample;
				}
				outChunks.push_back(chunk);
			}
		}
		
		void writeWAV(string filename)
		{
			WAV_OUT outfile(sampleRate_, bitsPerSample_, channels_);
			for(int i = 0; i < samples_.size(); ++i)
				outfile.write_current_output(samples_[i] * 32768.0);
			outfile.save_wave_file(filename.c_str());
		}
	};