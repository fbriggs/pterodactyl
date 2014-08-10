// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#include "pFileIO_BMP.h"
#include "pColor.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <unistd.h>
#include <stack>
#include <utility>
#include <fstream>

///----- defines for Win32 data types -----///

#if !defined(WIN32) || defined(_CONSOLE)

typedef unsigned char BYTE;	
typedef unsigned short int WORD;
typedef unsigned int DWORD;		
typedef int LONG;

typedef struct tagBITMAPFILEHEADER {
    WORD bfType;
    DWORD bfSize;
    WORD bfReserved1;
    WORD bfReserved2;
    DWORD bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
    DWORD biSize;
    LONG biWidth;
    LONG biHeight;
    WORD biPlanes;
    WORD biBitCount;
    DWORD biCompression;
    DWORD biSizeImage;
    LONG biXPelsPerMeter;
    LONG biYPelsPerMeter;
    DWORD biClrUsed;
    DWORD biClrImportant;
} BITMAPINFOHEADER;

// constants for the biCompression field
#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

typedef struct tagRGBTRIPLE {
    BYTE rgbtBlue;
    BYTE rgbtGreen;
    BYTE rgbtRed;
} RGBTRIPLE;

typedef struct { // tagRGBQUAD
    BYTE rgbBlue;
    BYTE rgbGreen;
    BYTE rgbRed;
    BYTE rgbReserved;
} RGBQUAD;

#endif // !defined(WIN32) || defined(_CONSOLE)

// magic numbers
#define BMP_BF_TYPE 0x4D42	// word BM
#define BMP_BF_OFF_BITS 54	// 14 for file header + 40 for info header (not sizeof(), but packed size)
#define BMP_BI_SIZE 40		// packed size of info header

///----- helper functions -----///

// reads a WORD from a file in little endian format
static WORD WordReadLE(ifstream& fp)
{
    WORD lsb, msb;
	
    lsb = fp.get();
    msb = fp.get();
    return (msb << 8) | lsb;
}

// writes a WORD to a file in little endian format
static void WordWriteLE(WORD x, ofstream& fp)
{
    BYTE lsb, msb;
	
    lsb = (BYTE) (x & 0x00FF);
    msb = (BYTE) (x >> 8);
    fp.put(lsb);
    fp.put(msb);
}

// reads a DWORD word from a file in little endian format
static DWORD DWordReadLE(ifstream& fp)
{
    DWORD b1, b2, b3, b4;
	
    b1 = fp.get();
    b2 = fp.get();
    b3 = fp.get();
    b4 = fp.get();
    return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}

// writes a DWORD to a file in little endian format
static void DWordWriteLE(DWORD x, ofstream& fp)
{
    unsigned char b1, b2, b3, b4;
	
    b1 = (x & 0x000000FF);
    b2 = ((x >> 8) & 0x000000FF);
    b3 = ((x >> 16) & 0x000000FF);
    b4 = ((x >> 24) & 0x000000FF);
    fp.put(b1);
    fp.put(b2);
    fp.put(b3);
    fp.put(b4);
}

// reads a LONG word from a file in little endian format
static LONG LongReadLE(ifstream& fp)
{
    LONG b1, b2, b3, b4;
	
    b1 = fp.get();
    b2 = fp.get();
    b3 = fp.get();
    b4 = fp.get();
    return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}

// writes a LONG to a file in little endian format
static void LongWriteLE(LONG x, ofstream& fp)
{
    char b1, b2, b3, b4;
	
    b1 = (x & 0x000000FF);
    b2 = ((x >> 8) & 0x000000FF);
    b3 = ((x >> 16) & 0x000000FF);
    b4 = ((x >> 24) & 0x000000FF);
    fp.put(b1);
    fp.put(b2);
    fp.put(b3);
    fp.put(b4);
}

int bitcount(DWORD w)
{
	w = (0x55555555 & w) + (0x55555555 & (w>> 1));
	w = (0x33333333 & w) + (0x33333333 & (w>> 2));
	w = (0x0f0f0f0f & w) + (0x0f0f0f0f & (w>> 4));
	w = (0x00ff00ff & w) + (0x00ff00ff & (w>> 8));
	w = (0x0000ffff & w) + (0x0000ffff & (w>>16));
	return w;
}

///----- read and write BMPs -----///


pImageBytes pFileIO_BMP::flipVertical(pImageBytes& src)
{
	pImageBytes dest(src.w_, src.h_, src.channels_);
	
	for(int x = 0; x < src.w_; ++x)
	for(int y = 0; y < src.h_; ++y)
	for(int c = 0; c < src.channels_; ++c)
		dest.setPixel(x, src.h_ - y - 1, c, src.getPixel(x, y, c));
	
	return dest;
}

pImageBytes pFileIO_BMP::read(string filename)
{
	bool needsFlipVertical = false;
	
	
	ifstream fp(filename.c_str(), ios::in | ios::binary);
	assert(fp.good() && "readBMP error: file not found");
	
	BITMAPFILEHEADER bmfh;
    BITMAPINFOHEADER bmih;
	
	bmfh.bfType = WordReadLE(fp);
    bmfh.bfSize = DWordReadLE(fp);
    bmfh.bfReserved1 = WordReadLE(fp);
    bmfh.bfReserved2 = WordReadLE(fp);
    bmfh.bfOffBits = DWordReadLE(fp);
	
	assert(bmfh.bfType == BMP_BF_TYPE && "readBMP error: not a BMP?");
	
	bmih.biSize = DWordReadLE(fp);
    bmih.biWidth = LongReadLE(fp);
	
    bmih.biHeight = LongReadLE(fp);
	if( bmih.biHeight < 0 )
	{
		needsFlipVertical = true;
		bmih.biHeight = abs(bmih.biHeight);
	}

	
	
    bmih.biPlanes = WordReadLE(fp);
    bmih.biBitCount = WordReadLE(fp);
    bmih.biCompression = DWordReadLE(fp);
    bmih.biSizeImage = DWordReadLE(fp);
    bmih.biSizeImage = bmih.biHeight*bmih.biWidth;  
    bmih.biXPelsPerMeter = LongReadLE(fp);
    bmih.biYPelsPerMeter = LongReadLE(fp);
    bmih.biClrUsed = DWordReadLE(fp);
    bmih.biClrImportant = DWordReadLE(fp);
	
	assert(bmih.biSize == BMP_BI_SIZE && "readBMP error: header info inconsistency");
	assert((bmih.biWidth > 0)  && (bmih.biHeight > 0) && (bmih.biPlanes == 1));
	
	// init the image
	pImageBytes img(bmih.biWidth, bmih.biHeight, 3);
	
	int scanlinelength;
    if ((img.w_ * bmih.biBitCount) % 32 == 0)	scanlinelength = img.w_ * bmih.biBitCount / 8;
    else										scanlinelength = (img.w_ * bmih.biBitCount / 32 + 1) * 4;
	
    RGBQUAD    *palette;
    BYTE       *scanlineByte;
    WORD       *scanlineWord;
    DWORD      *scanlineDword;
    RGBTRIPLE  *scanline24;
    RGBQUAD    *scanline32;
    int         index;
    pRGB		pixel;
    WORD		red, green, blue;
    DWORD       bluemask, greenmask, redmask;
    int         bluewidth, greenwidth;
    bool        done;
	
	switch (bmih.biBitCount)
    {
	case 1:	// 1 bit - monochrome, index mode
		// read the palette
		palette = new RGBQUAD[2];
		memset(palette, 0, 2 * sizeof(RGBQUAD));
		fp.read((char*) palette, 2 * sizeof(RGBQUAD));
		
		// read the data line by line
		scanlineByte = new BYTE[scanlinelength];
		fp.seekg((long) bmfh.bfOffBits, ifstream::beg);
		
		for(int y = 0; y < img.h_; ++y) 
		{
			fp.read((char*) scanlineByte, scanlinelength);
			assert(!fp.fail() && "unknown stream error");
			
			for(int x = 0; x < img.w_; ++x)
			{
				index = (scanlineByte[x/8] >> (7 - (x % 8))) & 0x01;
				
				pixel.r_ = palette[index].rgbRed   / 255.0;
				pixel.g_ = palette[index].rgbGreen / 255.0;
				pixel.b_ = palette[index].rgbBlue  / 255.0;
				img.setPixelRGB(x, y, pixel);
			}
		}
		
		delete[] scanlineByte;
		delete[] palette;
		break;
		
	case 4:	assert(false && "readBMP error: 4-bit format detected, not implemented... seriously, who uses a 4-bit BMP?"); break;
	case 8:	// 8 bit - 256 color, index mode
		// read the palette 
		palette = new RGBQUAD[256];
		memset(palette, 0, 256 * sizeof(RGBQUAD));
		fp.read((char*) palette, 256 * sizeof(RGBQUAD));
		
		// detect a grayscale palette, and modify image accordingly
		done = false;
		for(int i = 0; !done && i < 256; ++i)
			done = done || (palette[i].rgbRed  != palette[i].rgbGreen) || (palette[i].rgbBlue != palette[i].rgbGreen);
		
		if (!done)
		{
			img.clear();
			img.init(bmih.biWidth, bmih.biHeight, 3);
		}
		
		// read the data line by line
		if(bmih.biCompression == BI_RGB)		// uncompressed data
		{
			scanlineByte = new BYTE[scanlinelength];
			fp.seekg((long) bmfh.bfOffBits, ifstream::beg);
			
			for(int y = 0; y < img.h_; ++y) 
			{
				fp.read((char*) scanlineByte, scanlinelength);
				
				assert(!fp.fail() && "unknown stream error");
				
				for(int x = 0; x < img.w_; ++x)
				{
					index = scanlineByte[x];
					pixel.r_ = palette[index].rgbRed   / 255.0;
					pixel.g_ = palette[index].rgbGreen / 255.0;
					pixel.b_ = palette[index].rgbBlue  / 255.0;
					img.setPixelRGB(x, y, pixel);
				}
			}
			delete[] scanlineByte;
		}
		else if(bmih.biCompression == BI_RLE8)	// 8-bit RLE compression
		{
			unsigned char rleCode[2];
			int curx = 0;
			int cury = 0;
			done = false;
			
			fp.seekg((long) bmfh.bfOffBits, ifstream::beg);
			while(!done && fp)
			{
				fp.read((char*)rleCode, 2);
				
				if(rleCode[0] == 0 && rleCode[1] < 3)	// escape
				{
					if(rleCode[1] == 0)
					{
						curx = 0;
						++cury;
						if (cury >= img.h_)
							done = true;
					}
					else if(rleCode[1] == 1)
					{
						done = true;
					}
					else
					{
						curx += (int) fp.get();
						cury += (int) fp.get();
					}
				}
				else if(rleCode[0] == 0)			// absolute mode
				{
					for(int i = 0; i < rleCode[1]; ++i)
					{
						index = fp.get();
						
						pixel.r_ = palette[index].rgbRed   / 255.0;
						pixel.g_ = palette[index].rgbGreen / 255.0;
						pixel.b_ = palette[index].rgbBlue  / 255.0;
						img.setPixelRGB(curx, cury, pixel);
						
						++curx;
					}
					
					if (rleCode[1] % 2 != 0)
						fp.get();
				}
				else						// encoded mode
				{
					pixel.r_ = palette[rleCode[1]].rgbRed   / 255.0;
					pixel.g_ = palette[rleCode[1]].rgbGreen / 255.0;
					pixel.b_ = palette[rleCode[1]].rgbBlue  / 255.0;
					
					for(int i = 0; i < rleCode[0]; ++i)
					{
						img.setPixelRGB(curx, cury, pixel);
						++curx;
					}
				}
			}
			
			assert(!fp.fail() && "unknown stream error");
		}
		
		delete[] palette;
		break;
			
	case 16:	// 16 bit - 2^16 color, rgb mode 
		// set the color masks and shifts
		if(bmih.biCompression == BI_RGB)		// using default values
		{
			bluemask   = 0x001F;
			bluewidth  = 5;
			greenmask  = 0x03E0;
			greenwidth = 5;
			redmask    = 0x7C00;
		}
		else if(bmih.biCompression == BI_BITFIELDS)	// user specified
		{
			redmask    = DWordReadLE(fp);
			greenmask  = DWordReadLE(fp);
			bluemask   = DWordReadLE(fp);
			bluewidth  = bitcount(bluemask);
			greenwidth = bitcount(greenmask);
		}
		
		// get data line by line
		scanlineWord = new WORD[scanlinelength];
		fp.seekg((long) bmfh.bfOffBits, ifstream::beg);
		
		for(int y = 0; y < img.h_; ++y) 
		{
			fp.read((char*) scanlineWord, scanlinelength);
			assert(!fp.fail() && "unknown stream error");
			
			for(int x = 0; x < img.w_; ++x)
			{
				red   = (scanlineWord[x] & redmask) >> (bluewidth + greenwidth);
				green = (scanlineWord[x] & greenmask) >> bluewidth;
				blue  = (scanlineWord[x] & bluemask);
				
				pixel.r_ = red   / 255.0;
				pixel.g_ = green / 255.0;
				pixel.b_ = blue  / 255.0;
				img.setPixelRGB(x, y, pixel);
			}
		}
		
		delete[] scanlineWord;
		break;
		
	case 24:	// 24 bit - 2^24 color, rgb mode 
		// read the data line by line
		scanline24 = new RGBTRIPLE[scanlinelength];
		fp.seekg((long) bmfh.bfOffBits, ifstream::beg);
		
		for(int y = 0; y < img.h_; ++y) 
		{
			fp.read((char*) scanline24, scanlinelength);
			
			assert(!fp.fail() && "unknown stream error");
			
			for(int x = 0; x < img.w_; ++x)
			{
				pixel.r_ = scanline24[x].rgbtRed   / 255.0;
				pixel.g_ = scanline24[x].rgbtGreen / 255.0;
				pixel.b_ = scanline24[x].rgbtBlue  / 255.0;
				img.setPixelRGB(x, y, pixel);
			}
		}
		
		delete[] scanline24;
		break;

	case 32:	// 32 bit - 2^32 color, rgb mode
		// read data line by line
		if (bmih.biCompression == BI_RGB)		// default encoding
		{
			scanline32 = new RGBQUAD[scanlinelength];
			fp.seekg((long) bmfh.bfOffBits, ifstream::beg);
			
			for(int y = 0; y < img.h_; ++y) 
			{
				fp.read((char*) scanline32, scanlinelength);
				assert(!fp.fail() && "unknown stream error");
				
				for(int x = 0; x < img.w_; ++x)
				{
					pixel.r_ = scanline32[x].rgbRed   / 255.0;
					pixel.g_ = scanline32[x].rgbGreen / 255.0;
					pixel.b_ = scanline32[x].rgbBlue  / 255.0;
					img.setPixelRGB(x, y, pixel);
				}
			}
			
			delete[] scanline32;
		}
		else if(bmih.biCompression == BI_BITFIELDS)	// user specified
		{
			// get masks and shifts
			redmask    = DWordReadLE(fp);
			greenmask  = DWordReadLE(fp);
			bluemask   = DWordReadLE(fp);
			bluewidth  = bitcount(bluemask);
			greenwidth = bitcount(greenmask);
			
			scanlineDword = new DWORD[scanlinelength];
			fp.seekg((long) bmfh.bfOffBits, ifstream::beg);
			
			for(int y = 0; y < img.h_; ++y) 
			{
				fp.read((char*) scanlineDword, scanlinelength);
				assert(!fp.fail() && "unknown stream error");
				
				for(int x = 0; x < img.w_; ++x)
				{
					red   = (scanlineDword[x] & redmask) >> (bluewidth + greenwidth);
					green = (scanlineDword[x] & greenmask) >> bluewidth;
					blue  = (scanlineDword[x] & bluemask);
					
					pixel.r_ = red   / 255.0;
					pixel.g_ = green / 255.0;
					pixel.b_ = blue  / 255.0;
					img.setPixelRGB(x, y, pixel);
				}
			}
			
			delete[] scanlineDword;
		}
		break;
    }
	
	
	if( needsFlipVertical )
		return flipVertical(img);
	
	return img;
}

void pFileIO_BMP::write(pImageBytes& img, string filename)
{
	ofstream fp(filename.c_str(), ios::out | ios::binary);
	assert(fp.good() && "writeBMP error: could not create file");
	
	// write header
	
	BITMAPFILEHEADER bmfh;
    BITMAPINFOHEADER bmih;
    int lineLength;

	if (img.channels_ == 1) lineLength = img.w_; else lineLength = img.w_ * 3;
	if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
	
	bmfh.bfType = BMP_BF_TYPE;
    bmfh.bfSize = BMP_BF_OFF_BITS + lineLength * img.h_;
    bmfh.bfReserved1 = 0;
    bmfh.bfReserved2 = 0;
    bmfh.bfOffBits = BMP_BF_OFF_BITS;
	
    if (img.channels_ == 1)
		bmfh.bfOffBits += 256 * 4;
	
    WordWriteLE(bmfh.bfType, fp);
    DWordWriteLE(bmfh.bfSize, fp);
    WordWriteLE(bmfh.bfReserved1, fp);
    WordWriteLE(bmfh.bfReserved2, fp);
    DWordWriteLE(bmfh.bfOffBits, fp);
	
    bmih.biSize = BMP_BI_SIZE;
    bmih.biWidth = img.w_;
    bmih.biHeight = img.h_;
    bmih.biPlanes = 1;
    bmih.biBitCount = (img.channels_ == 1) ? 8 : 24;
    bmih.biCompression = BI_RGB;
    bmih.biSizeImage = lineLength * (DWORD) bmih.biHeight;
    bmih.biXPelsPerMeter = 2925;
    bmih.biYPelsPerMeter = 2925;
    bmih.biClrUsed = (img.channels_ == 1) ? 256 : 0;
    bmih.biClrImportant = 0;
	
    DWordWriteLE(bmih.biSize, fp);
    LongWriteLE(bmih.biWidth, fp);
    LongWriteLE(bmih.biHeight, fp);
    WordWriteLE(bmih.biPlanes, fp);
    WordWriteLE(bmih.biBitCount, fp);
    DWordWriteLE(bmih.biCompression, fp);
    DWordWriteLE(bmih.biSizeImage, fp);
    LongWriteLE(bmih.biXPelsPerMeter, fp);
    LongWriteLE(bmih.biYPelsPerMeter, fp);
    DWordWriteLE(bmih.biClrUsed, fp);
    DWordWriteLE(bmih.biClrImportant, fp);
	
    // write pixels
	
    if (img.channels_ == 1)
    {
		// write 8-bit grayscale palette
		unsigned char palettecolor[4];
		for(int i = 0; i < 256; ++i)
		{
			memset(palettecolor, (unsigned char) i, 4);
			fp.write((char*)palettecolor, 4);
		}
		
		// write image data
		for(int y = 0; y < img.h_; ++y) 
		{
			int nbytes = 0;
			for(int x = 0; x < img.w_; ++x) 
			{
				fp.put(img.getPixel(x, y, 0)); nbytes++;
			}
			
			while((nbytes % 4) != 0) 
			{
				fp.put((unsigned char) 0);
				nbytes++;
			}
		}
    }
    else
    {
		for(int y = 0; y < img.h_; ++y) 
		{
			int nbytes = 0;
			for(int x = 0; x < img.w_; ++x) 
			{
				fp.put(img.getPixel(x, y, 2)); nbytes++;
				fp.put(img.getPixel(x, y, 1)); nbytes++;
				fp.put(img.getPixel(x, y, 0)); nbytes++;
			}
			
			while((nbytes % 4) != 0) 
			{
				fp.put((unsigned char) 0);
				nbytes++;
			}
		}
    }
	
    assert(!fp.fail() && "unknown stream error");
}