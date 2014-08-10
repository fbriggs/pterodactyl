// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// this code handles reading and writing BMP format images. warning: the cpp is pretty nasty (so is the BMP file format)!
#pragma once

#include "pImageBytes.h"
#include <string>
using namespace std;

class pFileIO_BMP
{
public:
	// load the BMP file at filename and return a pointer to a newly allocated pImageBytes.
	// the caller of this function is reponsible for deleting the returned image pointer!!!
	static pImageBytes read(string filename);

	// write a BMP file to the specified destination
	static void write(pImageBytes& img, string filename);
	
	// some BMP images are upside-down.... damn you BMP format
	static pImageBytes flipVertical(pImageBytes& src);
};