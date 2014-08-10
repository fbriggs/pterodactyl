// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// functions for extracting isolated segments from a spectrogram for further processing
#pragma once

#include "pImageBytes.h"
#include "pRectangle.h"

#include <stack>
using namespace std;

class p2DSegmentExtraction
{
public:
	// segments which don't meet these criteria will be discarded
	const static int minSegmentWidth = 7;
	const static int minSegmentHeight = 7;
	const static int minSegmentArea = 15 * 15;
	
	// take a spectrogram and a binary mask. output a vector of images, each corresponding
	// to one segment cutout. the mask is passed by value because it will be modified
	static void extractSegments(pImageBytes& spectrogram, pImageBytes binaryMask, vector<pImageBytes>& extracted, vector<pImageBytes>& croppedMasks, vector<pRectangle>& rects)
	{
		while(true)
		{			
			// find the first pixel in the mask that is white
			int x, y;
			for(y = 0; y < binaryMask.h_; ++y)
			{
				for(x = 0; x < binaryMask.w_; ++x)
				{
					if( binaryMask.getPixel(x, y, 0) > 0 )
						goto found_a_pixel;
				}
			}
			return;
			found_a_pixel: 
			int minX, maxX, minY, maxY;
			
			minX = x;
			minY = y;
			maxX = x;
			maxY = y;
			
			pImageBytes componentMask	= extractComponentMask(binaryMask, x, y, minX, maxX, minY, maxY);
			pImageBytes cropped			= extractSpectrogramBoxFromComponentMask(spectrogram, componentMask);
			pImageBytes croppedMask		= extractCroppedMaskFromComponentMask(componentMask);
			
			// get rid of segments that are too small
			if( croppedMask.w_ < 1 || croppedMask.h_ < 1 ) continue; 
			if( maxX - minX < minSegmentWidth || maxY - minY < minSegmentHeight || ( maxX - minX) * (maxY - minY) < minSegmentArea ) continue;
			
			extracted.push_back(cropped);
			croppedMasks.push_back(croppedMask);
			rects.push_back(pRectangle(minX, minY, maxX, maxY));
		}
	}
	
	// search for all connected pixels and generate the segment mask
	// also zero out any pixels in the binary mask from this components
	static pImageBytes extractComponentMask(pImageBytes& binaryMask, int startX, int startY, int& minX, int& maxX, int& minY, int& maxY)
	{
		pImageBytes componentMask(binaryMask.w_, binaryMask.h_, 1);
		stack< pair<int, int> > dfs;
		dfs.push(pair<int,int>(startX, startY));
		while(!dfs.empty())
		{
			pair<int,int> currLoc = dfs.top(); 
			dfs.pop();
			int currX = currLoc.first;
			int currY = currLoc.second;
			
			minX = min(currX, minX);
			minY = min(currY, minY);
			maxX = max(currX, maxX);
			maxY = max(currY, maxY);
			
			binaryMask.setPixel(currX, currY, 0, 0);
			componentMask.setPixel(currX, currY, 0, 255);
			
			if( currX + 1 < binaryMask.w_			&& binaryMask.getPixel_(currX + 1, currY, 0) > 0 ) dfs.push(pair<int,int>(currX + 1, currY));
			if( currX - 1 >= 0						&& binaryMask.getPixel_(currX - 1, currY, 0) > 0 ) dfs.push(pair<int,int>(currX - 1, currY));
			if( currY + 1 < binaryMask.h_			&& binaryMask.getPixel_(currX, currY + 1, 0) > 0 ) dfs.push(pair<int,int>(currX, currY + 1));
			if( currY + 1 >= 0						&& binaryMask.getPixel_(currX, currY - 1, 0) > 0 ) dfs.push(pair<int,int>(currX, currY - 1));
		}
		
		return componentMask;
	}
	
	
	// takes a spectrogram and mask, and returns a cropped subimage of the spectrogram 
	// that contains just the pixels from the mask
	static pImageBytes extractSpectrogramBoxFromComponentMask(pImageBytes& spectrogram, pImageBytes& componentMask)
	{
		// figure out the min-x and max-x values included in the component mask
		int minX = spectrogram.w_-1;
		int maxX = 0;
		for(int y = 0; y < componentMask.h_; ++y)
		for(int x = 0; x < componentMask.w_; ++x)
		{
			if( componentMask.getPixel(x, y, 0) > 0 )
			{
				minX = min(x, minX);
				maxX = max(x, maxX);
			}
		}
		
		int cropW = maxX - minX;
		pImageBytes cropped(cropW, spectrogram.h_, 1);
		for(int y = 0; y < cropped.h_; ++y)
		for(int x = 0; x < cropped.w_; ++x)
		{
			int srcX = minX + x;
			int srcY = y;
			if( componentMask.getPixel(srcX, srcY, 0) > 0 )
				cropped.setPixel(x, y, 0, spectrogram.getPixel(srcX, srcY, 0));
		}
		return cropped;
	}
	
	// component mask is a big image with a small segment in it. this function extracts that segment
	// from the full image by discarding parts of the spectrogram to the left and right of it (above and below are preserved).
	static pImageBytes extractCroppedMaskFromComponentMask(pImageBytes& componentMask)
	{
		// figure out the min-x and max-x values included in the component mask
		int minX = componentMask.w_-1;
		int maxX = 0;
		for(int y = 0; y < componentMask.h_; ++y)
			for(int x = 0; x < componentMask.w_; ++x)
			{
				if( componentMask.getPixel(x, y, 0) > 0 )
				{
					minX = min(x, minX);
					maxX = max(x, maxX);
				}
			}
		
		int cropW = maxX - minX;
		pImageBytes cropped(cropW, componentMask.h_, 1);
		for(int y = 0; y < cropped.h_; ++y)
			for(int x = 0; x < cropped.w_; ++x)
			{
				int srcX = minX + x;
				int srcY = y;
				if( componentMask.getPixel(srcX, srcY, 0) > 0 )
					cropped.setPixel(x, y, 0, 255);
			}
		
		return cropped;
	}
	
	// we started the process of extracting a segment from the rest of the spectrogram with extractCroppedMaskFromComponentMask, 
	// which discarded the parts to the left and right. this function discards the parts above and below (only used for some stuff like visualization).
	static pImageBytes packVertical(pImageBytes& src, pImageBytes& mask)
	{
		int minY = src.h_ - 1;
		int maxY = 0;
		
		for(int y = 0; y < src.h_; ++y)
		for(int x = 0; x < src.w_; ++x)
		{
			if( mask.getPixel(x, y, 0) > 0 )
			{
				minY = min(minY, y);
				maxY = max(maxY, y);
			}
		}
		
		pImageBytes packed(src.w_, maxY - minY, 1);
		
		for(int y = 0; y < packed.h_; ++y)
		for(int x = 0; x < packed.w_; ++x) 
			packed.setPixel(x,y, 0, src.getPixel(x, y + minY, 0));
		
		return packed;
	}
};