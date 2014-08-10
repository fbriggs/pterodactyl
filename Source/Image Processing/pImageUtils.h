// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once

#include "pTuple.h"
#include "pImageBytes.h"
#include "pMatrix.h"
#include <cmath>

struct pImageBytes_HeightSortComparator // used to sort images by height 
{ 
	bool operator ()(const pImageBytes* a, const pImageBytes* b) { return a->h_ < b->h_; } 
};

struct pImageBytes_WidthSortComparator // used to sort images by height 
{ 
	bool operator ()(const pImageBytes* a, const pImageBytes* b) { return a->w_ < b->w_; } 
};

// a collection of miscalaneous utils for working with images
class pImageUtils
{
public:
	// convert an image to a matrix
	static pMatrix image2matrix(pImageBytes& img)
	{
		pMatrix m(img.w_, img.h_);
		
		for(int i = 0; i < img.w_; ++i)
		for(int j = 0; j < img.h_; ++j)
			m.data_[i][j] = img.getPixelF(i,j,0);
		
		return m;
	}
	
	// convert a matrix to an image
	static pImageBytes matrix2image(pMatrix m)
	{
		m.normalizeTo01();
		
		pImageBytes img(m.rows_, m.cols_, 1);
		for(int i = 0; i < m.rows_; ++i)
		for(int j = 0; j < m.cols_; ++j)
			img.setPixelF(i,j,0, m.data_[i][j]);
		
		return img;
	}
	
	
	// clamps x to the range [0,255]
	static float clamp255(float x)									
	{ 
		return min(255.0f, max(0.0f, x)); 
	}
	
	// clamps x to the range [0,1]
	static float clamp01(float x)									
	{ 
		return min(1.0f, max(0.0f, x)); 
	}
	
	// return a greyscale, 1-channel version of src
	static pImageBytes greyscale(pImageBytes src)
	{
		pImageBytes g(src.w_, src.h_, 1);
		
		for(int y = 0; y < src.h_; ++y)
		for(int x = 0; x < src.w_; ++x)
		{
			float avg = 0;
			for(int i = 0; i < src.channels_; ++i)
				avg += src.getPixel(x, y, i);
			
			g.setPixel(x, y, 0, avg / (float)src.channels_);
		}
		
		return g;
	}
	
	// takes a greyscale image and returns a 3-channel version
	static pImageBytes colorize(pImageBytes src)
	{
		pImageBytes g(src.w_, src.h_, 3);
		for(int y = 0; y < src.h_; ++y)
		for(int x = 0; x < src.w_; ++x)
		for(int c = 0; c < 3; ++c)
			g.setPixel(x, y, c, src.getPixel(x, y, 0));
			
		return g;
	}
	
	// returns a difference image between the two input images
	static pImageBytes difference(pImageBytes& imgA, pImageBytes& imgB)
	{
		assert(imgA.w_ == imgB.w_ && imgA.h_ == imgB.h_ && imgA.channels_ == imgB.channels_);
		
		pImageBytes diffImg(imgA.w_, imgA.h_, imgA.channels_);
		
		for(int y = 0; y < imgA.h_; ++y)
		for(int x = 0; x < imgA.w_; ++x)
		for(int c = 0; c < imgA.channels_; ++c)
			diffImg.setPixel(x, y, c, fabs(imgA.getPixel(x, y, c) - imgB.getPixel(x, y, c)));
		
		return diffImg;
	}
	
	static pImageBytes threshold(pImageBytes& src, float t)
	{
		pImageBytes results(src.w_, src.h_, 1);
		
		for(int y = 0; y < src.h_; ++y)
		for(int x = 0; x < src.w_; ++x)
			results.setPixel(x, y, 0, src.getPixelF(x, y, 0) >= t ? 255 : 0);
		
		return results;
	}
	
	// find the max-value point
	static pTuple<int,int> findMaxima(pImageBytes& img)
	{
		assert(img.channels_ == 1);
		
		float maxVal = 0;
		pTuple<int,int> maxPoint = make_tuple(0,0);
		
		for(int x = 0; x < img.w_; ++x)
		for(int y = 0; y < img.h_; ++y)
		{
			if( img.getPixel(x, y, 0) > maxVal )
			{
				maxVal = img.getPixel(x, y, 0);
				maxPoint = make_tuple(x, y);
			}
		}
		
		return maxPoint;
	}
	
	static void drawMaskOutline(pImageBytes& img, pImageBytes& mask, pRGB color)
	{
		pImageBytes perimImg = perimeterImage(mask);
		
		for(int x = 0; x < img.w_; ++x)
			for(int y = 0; y < img.h_; ++y)
			{
				if( perimImg.getPixel(x, y, 0) > 0 )
				{
					pRGB rgb = img.getPixelRGB(x,y);
					
					img.setPixelRGB(x, y, pRGB(.2*rgb.r_ + .8*color.r_, .2*rgb.g_ + .8*color.g_, .2*rgb.b_ + .8*color.b_));
				}
			}
	}
	
	// generate a binary image of the perimiter of mask
	static pImageBytes perimeterImage(pImageBytes& mask)
	{
		pImageBytes perim(mask.w_, mask.h_, 1);
		for(int x = 0; x < mask.w_; ++x)
			for(int y = 0; y < mask.h_; ++y)
			{
				// check if all the surrounding pixels are in the same cluster
				int matchVal = mask.getPixel(x, y, 0);
				bool allSame = true;
				for(int u = x - 1; u <= x; ++u)
					for(int v = y - 1; v <= y; ++v)
						if( mask.getPixel_(u, v, 0) != matchVal )
							allSame = false;
				
				if( !allSame ) 
					perim.setPixel(x, y, 0, 255);
			}
		return perim;
	}
	
	// compute the area of a mask
	static float maskArea(pImageBytes& mask)
	{
		float sum = 0;
		for(int x = 0; x < mask.w_; ++x)
			for(int y = 0; y < mask.h_; ++y)
				if( mask.getPixel(x, y, 0) > 0 )
					sum += 1;
		return sum;
	}
	
	// compute the perimeter of a mask
	static float maskPerimeter(pImageBytes& mask)
	{
		float sum = 0;
		for(int x = 0; x < mask.w_; ++x)
		for(int y = 0; y < mask.h_; ++y)
		{
			// check if all the surrounding pixels are in the same cluster
			int matchVal = mask.getPixel(x, y, 0);
			bool allSame = true;
			for(int u = x - 1; u <= x; ++u)
			for(int v = y - 1; v <= y; ++v)
				if( mask.getPixel_(u, v, 0) != matchVal )
					allSame = false;
			
			if( !allSame ) sum += 1;
		}
		return sum;
	}
	
	static pImageBytes stackHorizontal(vector<pImageBytes>& images)
	{
		int newW = 0;
		int newH = 0;
		for(int i = 0; i < images.size(); ++i)
		{
			newW += images[i].w_;
			newH = max(newH, images[i].h_);
		}
		pImageBytes newImg(newW, newH, images[0].channels_);
		
		int startX = 0;
		for(int i = 0; i < images.size(); ++i)
		{
			for(int x = 0; x < images[i].w_; ++x)
			for(int y = 0; y < images[i].h_; ++y)
				newImg.setPixelRGB(x + startX, y, images[i].getPixelRGB(x,y));
			
			startX += images[i].w_;
		}
		
		return newImg;
	}
	
	static pImageBytes stackVertical(vector<pImageBytes>& images)
	{
		int newW = 0;
		int newH = 0;
		for(int i = 0; i < images.size(); ++i)
		{
			newW = max(newW, images[i].w_);
			newH += images[i].h_;
		}
		pImageBytes newImg(newW, newH, images[0].channels_);
		
		int startY = 0;
		for(int i = 0; i < images.size(); ++i)
		{
			for(int x = 0; x < images[i].w_; ++x)
			for(int y = 0; y < images[i].h_; ++y)
			for(int c = 0; c < images[0].channels_; ++c)
				newImg.setPixelF(x, y + startY, c, images[i].getPixelF(x,y,c));
				
			startY += images[i].h_;
		}
		return newImg;
	}

	// takes an array of images and arranges them in a space-efficient mosaic.
	// the image is formatted to a max width of maxW, and its height is determined by the content
	static pImageBytes irregularImageMosaic(vector<pImageBytes*>& images, int maxW)
	{	
		if( images.size() == 1 ) return *images[0];
		
		int padding = 6;
		
		// if max width isn't big enough, expand. this prevents the possibility of
		// maxW being too small to contain the images
		for(int i = 0; i < images.size(); ++i)
			maxW = max(maxW, images[i]->w_ + padding + 1);
		
		stable_sort(images.begin(), images.end(), pImageBytes_WidthSortComparator());
		stable_sort(images.begin(), images.end(), pImageBytes_HeightSortComparator());
		
		vector<pImageBytes> rows;
		
		pImageBytes padImg(maxW, padding, 1);
		//fillChannel(padImg, 0, 0.55);
		rows.push_back(padImg);
		
		int i = 0;
		while( i < images.size() )
		{		
			// assuming the current row starts with image i, find the maximum index of
			// image that can fit in the row
			int currRowWidth = 0;
			int maxIndexInRow = i;
			for(int j = i; j < images.size(); ++j)
			{
				currRowWidth += images[j]->w_+ padding;
				maxIndexInRow = j;
				if( currRowWidth >= maxW )
				{
					--maxIndexInRow;
					break;
				}
			}
		
			int maxH = 0;
			for(int j = i; j <= maxIndexInRow; ++j)
				maxH = max(maxH, images[j]->h_);
			
			pImageBytes rowImg(maxW, maxH, 1);
			//fillChannel(rowImg, 0, 0.55);
			
			// compute horizontal centering
			int trueRowWidth = 0;
			for(int j = i; j <= maxIndexInRow; ++j)
				trueRowWidth += images[j]->w_ + padding;
			
			int remainder = maxW - trueRowWidth;
				
			int xOffset = remainder / 2 + padding / 2;
			for(int j = i; j <= maxIndexInRow; ++j)
			{
				for(int x = 0; x < images[j]->w_; ++x)
				for(int y = 0; y < images[j]->h_; ++y)
					rowImg.setPixel(xOffset+x, y, 0, images[j]->getPixel(x,y,0));
				
				xOffset += images[j]->w_ + padding;
			}
			rows.push_back(rowImg);
			rows.push_back(padImg);
			
			i = maxIndexInRow + 1;
		}
		
		return stackVertical(rows);
	}
	
	// image gradients (sobel)
	static float dx(pImageBytes& img, int x, int y) 
	{ 
		return
		-img.getPixelF_(x - 1, y - 1, 0) - 2 * img.getPixelF_(x - 1, y, 0) - img.getPixelF_(x - 1, y + 1, 0)  +
		 img.getPixelF_(x + 1, y - 1, 0) + 2 * img.getPixelF_(x + 1, y, 0) + img.getPixelF_(x + 1, y + 1, 0);
	}
	
	static float dy(pImageBytes& img, int x, int y)
	{
		return
		-img.getPixelF_(x - 1, y - 1, 0) - 2 * img.getPixelF_(x, y - 1, 0) - img.getPixelF_(x + 1, y - 1, 0) +
		 img.getPixelF_(x - 1, y + 1, 0) + 2 * img.getPixelF_(x, y + 1, 0) + img.getPixelF_(x + 1, y + 1, 0);
	}
	
	static pVectorND grad(pImageBytes& img, int x, int y) { return pVectorND(dx(img, x, y), dy(img, x, y)); }
	
	static void fillChannel(pImageBytes& img, int channel, float val)
	{
		for(int x = 0; x < img.w_; ++x)
		for(int y = 0; y < img.h_; ++y)
			img.setPixelF(x, y, 0, val);
	}	
	
	static void fillRGB(pImageBytes& img, pRGB color)
	{
		fillChannel(img, 0, color.r_);
		fillChannel(img, 1, color.g_);
		fillChannel(img, 2, color.b_);
	}
};