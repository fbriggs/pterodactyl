// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// pImageBytes is the core class used for image processing operations. it stores the ray-pixel data of 1,3, or 4 channel image, and provides basic pixel setters and getters
#pragma once

#include "pColor.h"
#include <iostream>
#include <cstring>
#include <cassert>
#include <algorithm>

using namespace std;

// this class provides a low-level representation of an image
// using an array of bytes, and various setters and getters for manipulating pixels.
class pImageBytes
{
public:
	unsigned char* pixels_;	// the raw pixel data for the image
	int w_, h_, channels_;	// the width, height, and number of channels in the image
	int numBytes_;			// the total number of bytes in the image
	
	// deconstructor deallocates pixels_
	~pImageBytes() { clear(); } 
	
	// deallocate pixels and allow this image to be re-initialized
	void clear() { if( pixels_ != NULL ) delete[] pixels_; }
	
	// do-nothing constructor for deferred initialization
	pImageBytes() { pixels_ = NULL; }
	
	// construct an image with the specified width, height, and number of channels
	pImageBytes(int w, int h, int channels = 3) { init(w,h,channels); }
	
	// initialize the image by allocating pixels_ and setting all intensities to 0
	void init(int w, int h, int channels = 3)
	{
		w_ = w;
		h_ = h;
		channels_ = channels;
		numBytes_ = w * h * channels;
		pixels_ = new unsigned char[numBytes_];
		memset(pixels_, 0, w * h * channels);
	}
	
	// copy the contents of img into this image. both must have the 
	// same number of channels and dimensions
	void copy(const pImageBytes& img)
	{
		assert(w_ == img.w_ && h_ == img.h_ && channels_ == img.channels_);
		memcpy(pixels_, img.pixels_, numBytes_);
	}
	
	// copy constructor
	pImageBytes(const pImageBytes& img)
	{
		init(img.w_, img.h_, img.channels_);
		copy(img);
	}
	
	///----- misc utils  -----///
	
	// returns true if (x, y) is in the image and false otherwise
	bool testInside(int x, int y)							{ return (x >= 0 && y >= 0 && x < w_ && y < h_); }

	// get the byte index for pixel (x,y) channel c
	inline int byteIndex(int x, int y, int c = 0)			{ return (y * w_ + x) * channels_ + c; }
	
	///----- setters and getters for pixels  -----///
	
	// setters and getters for byte values in the range [0,255]
	inline unsigned char getPixel(int x, int y, int c = 0)			{ return pixels_[byteIndex(x,y,c)]; }
	inline void setPixel(int x, int y, int c, unsigned char v)		{ pixels_[byteIndex(x,y,c)] = v; }
	
	// setters and getters for float values in the range [0, 1].
	float getPixelF(int x, int y, int c = 0)				{ return getPixel(x,y,c) / 255.0f; }
	void setPixelF(int x, int y, int c, float v)			{ setPixel(x,y,c, v * 255.0f); }
	
	// accessors ending with an _ test if (x,y) is in the image fist before setting or getting.
	// setters do nothing if not, and getters return 0 if not.
	unsigned char getPixel_(int x, int y, int c = 0)		{ return testInside(x,y) ? getPixel(x, y, c) : 0; }
	float getPixelF_(int x, int y, int c = 0)				{ return getPixel_(x, y, c) / 255.0f; }
	void setPixel_ (int x, int y, int c, unsigned char v)	{ if( testInside(x, y) ) setPixel(x, y, c, v); }
	void setPixelF_(int x, int y, int c, float v)			{ if( testInside(x, y) ) setPixelF(x, y, c, v); }
	
	// setters and getters for RGB and RGBA data
	void setPixelRGB  (int x, int y, pRGB rgb)				{ setPixelF (x, y, 0, rgb.r_);	setPixelF (x, y, 1, rgb.g_);	setPixelF (x, y, 2, rgb.b_); }
	void setPixelRGBA (int x, int y, pRGBA rgb)				{ setPixelF (x, y, 0, rgb.r_);	setPixelF (x, y, 1, rgb.g_);	setPixelF (x, y, 2, rgb.b_);	setPixelF (x, y, 3, rgb.a_); }
	void setPixelRGB_ (int x, int y, pRGB rgb)				{ setPixelF_(x, y, 0, rgb.r_);	setPixelF_(x, y, 1, rgb.g_);	setPixelF_(x, y, 2, rgb.b_); }
	void setPixelRGBA_(int x, int y, pRGBA rgb)				{ setPixelF_(x, y, 0, rgb.r_);	setPixelF_(x, y, 1, rgb.g_);	setPixelF_(x, y, 2, rgb.b_);	setPixelF_(x, y, 3, rgb.a_); }
	
	pRGB  getPixelRGB  (int x, int y)						{ return pRGB (getPixelF (x, y, 0),	getPixelF (x, y, 1),	getPixelF (x, y, 2)); }
	pRGBA getPixelRGBA (int x, int y)						{ return pRGBA(getPixelF (x, y, 0),	getPixelF (x, y, 1),	getPixelF (x, y, 2), getPixelF(x, y, 3)); }
	pRGB  getPixelRGB_ (int x, int y)						{ return pRGB (getPixelF_(x, y, 0),	getPixelF_(x, y, 1),	getPixelF_(x, y, 2)); }
	pRGBA getPixelRGBA_(int x, int y)						{ return pRGBA(getPixelF_(x, y, 0),	getPixelF_(x, y, 1),	getPixelF_(x, y, 2), getPixelF(x, y, 3)); }
	
	// get a pixel and use edge-value extension to handle pixels out of bounds
	unsigned char getPixel_extend(int x, int y, int c = 0)
	{
		// if the pixel is inside, get its normal value
		if( testInside(x,y) ) return getPixel(x, y, c);
		
		// corners
		if( x < 0 && y < 0 )		return getPixel(0, 0, c);
		if( x < 0 && y >= h_ )		return getPixel(0, h_ - 1, c);
		if( x >= w_ && y < 0 )		return getPixel(w_ - 1, 0, c);
		if( x >= w_ && y >= h_ )	return getPixel(w_ - 1, h_ - 1, c);
		
		// sides
		if( x < 0 )					return getPixel(0, y, c);
		if( y < 0 )					return getPixel(x, 0, c);
		if( x >= w_ )				return getPixel(w_ - 1, y, c);
		if( y >= h_ )				return getPixel(x, h_ - 1, c);
		
		assert(false); // the above cases are exhaustive...
	}
	
	float getPixelF_extend(int x, int y, int c = 0)
	{
		return getPixel_extend(x, y, c) / 255.0f;
	}
};