// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// a collection of utilities for working with C++ strings
#pragma once

#include <vector>
#include <string>
#include <ctype.h>

using namespace std;

class pStringUtils
{
public:
	// like split, but removes empty results
	static vector<string> splitNonEmpty(const string s, const string delims)
	{
		vector<string> partsNE;
		vector<string> parts = split(s, delims);
		for(int i = 0; i < parts.size(); ++i)
			if( parts[i] != "" ) partsNE.push_back(parts[i]);
		return partsNE;
	}
	
	// split s into a vector of strings by delim and writes the results to outParts (which
	// should be empty before calling this
	static vector<string> split(const string s, const char delim)
	{
		vector<string> parts;
		
		string part = "";
		unsigned int sLen = s.length();
		
		for(unsigned int i = 0; i < sLen; ++i)
		{
			if( s[i] == delim )
			{
				parts.push_back(part);
				part = "";
			}
			else
				part += s[i];
		}
		
		parts.push_back(part);
		return parts;
	}
	
	// split with a string specifying a list of delimeters
	static vector<string> split(const string s, const string delims)
	{
		return split(replaceChar(s, delims, delims[0]), delims[0]);
	}
	
	// split, then only take the first part
	static string firstPartOfSplit(const string s, const char delim)
	{
		vector<string> parts = split(s, delim);
		return parts[0];
	}
	
	static string firstPartOfSplit(const string s, const string delims)
	{
		vector<string> parts = split(s, delims);
		return parts[0];
	}
	
	static string secondPartOfSplit(const string s, const string delims)
	{
		vector<string> parts = split(s, delims);
		return parts[1];
	}
	
	static string removeChar(string s, char c)
	{
		string r = "";
		for(int i = 0; i < s.size(); ++i)
			if( s[i] != c ) r += s[i];
		return r;
	}
	
	
	static string replaceChar(string s, string charsToReplace, char c)
	{
		for(int i = 0; i < s.size(); ++i)
		{
			for(int j = 0; j < charsToReplace.size(); ++j)
				if( s[i] == charsToReplace[j] )
					s[i] = c;
		}
		return s;
	}

	// Conversion between strings and integers
	static string intToString(int x)
	{
		char buffer[1024];
		sprintf(buffer, "%d", x);
		return string(buffer);
	}
	
	static int stringToInt(string s)
	{
		return atoi( s.c_str() );
	}

	// Conversion between strings and floats.  Precision is the number of
	// digits after the decimal point (-1 is unlimited digits).
	static string floatToString(float x, int precision = -1)
	{
		char buffer[1024];
		sprintf(buffer, "%f", x);
		
		string result(buffer);
		size_t decimalPos = result.find_first_of('.');
		
		if(precision < 0)
			return result;
		else
		{
			if(decimalPos == string::npos)
				return result;
			else
				return result.substr(0, decimalPos + 1 + precision);
		}
	}
	
	static float stringToFloat(string s)
	{
		return atof(s.c_str());
	}

	// Returns true if c is a space, tab character or newline
	static bool isWhitespace(char c)
	{
		return c == ' ' || c == '\n' || c == '\r' || c == '\t';
	}

	// Returns the argument string stripped of leading and trailing whitespace
	// and newline characters. If the string is entirely whitespace or newline
	// characters, the empty string is returned.
	static string pack(const string &val)
	{
		string result = val;
		
		size_t first_char = result.find_first_not_of(" \t\n");
		size_t last_char = result.find_last_not_of(" \t\n");
		
		if(first_char == string::npos || last_char == string::npos) {
			result = "";
		} else {
			result.erase(0, first_char);
			result.erase(last_char - first_char + 1, string::npos);
		}
		return result;
	}
};
