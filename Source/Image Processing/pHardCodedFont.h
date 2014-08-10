// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// this is a small-font text rendering system for pImageBytes. it doesn't require any resources or dependencies; bitmaps for each character
// are hard-coded into pHardCodedFont.cpp
#pragma once

#include "pImageBytes.h"

#include <iostream>
using namespace std;

class pHardCodedFont
{
public:
	// draw a string s in the image img at (x,y) using color
	static void drawString(pImageBytes& img, pRGB color, string s, int x, int y, bool whiteBox = false)
	{
		for(int i = 0; i < s.size(); ++i)
		{
			if( whiteBox )
			{
				for(int a = 0; a < 8; ++a)
				for(int b = 0; b < 10; ++b)
					img.setPixelRGB_(x + a  + i * 8, y + b, pRGB(1,1,1));
			}
			
			
			for(int a = 0; a < 8; ++a)
			for(int b = 0; b < 8; ++b)
				if( s[i] != ' ' && getCharPixel(s[i], a, b))
					img.setPixelRGB_(x + a  + i * 8, y + b, color);
		}
	}
	
	// c is a character in the limited character set for the font
	// x, y is the index of a pixel in the character (all characters are 8x8)
	static int getCharPixel(char c, int x, int y);
	
	const static int numCharsInFont_ = 73; // max char index in pirate encoding is this -1
	
	// convert a char in the ascii format to the special pirate encoding
	// used just by this class. the point of doing this is that the font
	// doesn't contain every ascii character, so these indices skip gaps
	static int asciiToPirate(char ascii)
	{
		if( ascii >= 'A' && ascii <= 'Z' )	return ascii - 'A';
		if( ascii >= 'a' && ascii <= 'z' )	return ascii - 'a' + 26;
		if( ascii >= '0' && ascii <= '9' )	return ascii - '0' + 52;
		
		switch(ascii)
		{
			case '.':	return 61 + 1;
			case ',':	return 61 + 2;
			case '\'':	return 61 + 3;
			case '"':	return 61 + 4;
			case '!':	return 61 + 5;
			case '?':	return 61 + 6;
			case ':':	return 61 + 7;
			case '*':	return 61 + 8;
			case '-':	return 61 + 9;
			case '+':	return 61 + 10;
			case '=':	return 61 + 11;
			default: cout << "fail on " << ascii << endl; assert(false && "asciiToPirate character not found");
		}
	}
	
	// convert a pirate code to a character
	static char pirateToAscii(int pirate)
	{
		if( pirate >= 0 && pirate < 26 )	return 'A' + pirate;
		if( pirate >= 26 && pirate < 52 )	return 'a' + pirate - 26;
		if( pirate >= 52 && pirate < 62 )	return '0' + pirate - 52;
		
		switch(pirate)
		{
		case 61 + 1:		return '.';
		case 61 + 2:		return ',';
		case 61 + 3:		return '\'';
		case 61 + 4:		return '"';
		case 61 + 5:		return '!';
		case 61 + 6:		return '?';
		case 61 + 7:		return ':';
		case 61 + 8:		return '*';
		case 61 + 9:		return '-';
		case 61 + 10:		return '+';
		case 61 + 11:		return '=';
		default: assert(false && "pirateToAscii code note found");
		}
	}
};
