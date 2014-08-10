// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// this file defines data structures and methods related to colors
#pragma once

#include <math.h>
#include <algorithm>
using namespace std;

// a 3-channel RGB color
class pRGB
{
public:
	float r_, g_, b_;
	
	pRGB() { r_ = 0; g_ = 0; b_ = 0; }
	pRGB(float r, float g, float b) { r_ = r; g_ = g; b_ = b; }
};

// a 4-channel RGBA color
class pRGBA
{
public: 
	float r_, g_, b_, a_;
	
	pRGBA() { r_ = 0; g_ = 0; b_ = 0; a_ = 0; }
	pRGBA(float r, float g, float b, float a) { r_ = r; g_ = g; b_ = b; a_ = a; }
};

// a 3-channel HSV color
class pHSV
{
public:
	float h_, s_, v_;
	
	pHSV() { h_ = 0; s_ = 0; v_ = 0; }
	pHSV(float h, float s, float v) { h_ = h; s_ = s; v_ = v; }
};


class pColorTransform
{
public:
	// transform hsv in [0,1] to rgb in [0,1]
	static pRGB hsv2rgb(float h, float s, float v)
	{
		float r, g, b;
		h = fmod(h * 360.0f, 360.0f);
		
		int i;
		float f, p, q, t;
		
		if( s == 0 )
		{
			r = g = b = v;
			return pRGB(r,g,b); 
		}
		
		h /= 60;			
		i = floor( h );
		f = h - i;						
		p = v * ( 1 - s );
		q = v * ( 1 - s * f );
		t = v * ( 1 - s * ( 1 - f ) );
		
		switch( i )
		{
			case 0:		r = v;	g = t;	b = p;	break;
			case 1:		r = q;	g = v;	b = p;	break;
			case 2:		r = p;	g = v;	b = t;	break;
			case 3:		r = p;	g = q;	b = v;	break;
			case 4:		r = t;	g = p;	b = v;	break;
			default:	r = v;	g = p;	b = q;	break;
		}
		return pRGB(r,g,b);
	}
	
	static float max3(float a, float b, float c) {return max(a, max(b,c)); }
	static float min3(float a, float b, float c) {return min(a, min(b,c)); }
	
	static pHSV rgb2hsv(pRGB color)
	{
		pHSV hsv;
		
		float minC, maxC, delta;
		minC = min3( color.r_, color.g_, color.b_ );
		maxC = max3( color.r_, color.g_, color.b_ );
		
		if( maxC == 0) 
			return pHSV(0,0,0);
		
		delta = maxC - minC;
		
		if( delta == 0)
			return pHSV(0, 0, maxC);
		
		hsv.v_ = maxC;					// v
		hsv.s_ = delta / maxC;			// s
		
		if( color.r_ == maxC )
			hsv.h_ = ( color.g_ - color.b_ ) / delta;		// between yellow & magenta
		else if( color.g_ == maxC )
			hsv.h_ = 2 + ( color.b_ - color.r_ ) / delta;	// between cyan & yellow
		else
			hsv.h_ = 4 + ( color.r_ - color.g_ ) / delta;	// between magenta & cyan
		hsv.h_ *= 60;										// degrees
		if( hsv.h_ < 0 )
			hsv.h_ += 360;
		
		
		return hsv;
	}
	
};