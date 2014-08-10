// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// this file defines an enumeration of distance metrics, and functions for computing
// distances in those metrics

#pragma once

#include <math.h>
#include <vector>
using namespace std;

enum pDistanceMetric
{
	DIST_ERROR, // this is if the parsing falied
	DIST_L1,
	DIST_L2,
	DIST_HELLINGER,
	DIST_KULLBACK_LIEBLER
};

pDistanceMetric distanceMetricStringToEnum(string distName)
{
	if( distName == "L1" )					return DIST_L1;
	if( distName == "L2" )					return DIST_L2;
	if( distName == "Hellinger" )			return DIST_HELLINGER;
	if( distName == "Kullback-Liebler" )	return DIST_KULLBACK_LIEBLER;
	return DIST_ERROR; // no matching distance found
}

string distanceMetricToString(pDistanceMetric d)
{
	switch( d )
	{
	case DIST_ERROR:						return "ERROR PARSING DISTANCE";
	case DIST_L1:							return "L1";
	case DIST_L2:							return "L2";
	case DIST_HELLINGER:					return "Hellinger";
	case DIST_KULLBACK_LIEBLER:				return "Kullback-Liebler";
	default:								assert(false);
	}
}

float distL1(vector<float>& x, vector<float>& y)
{
	float dist = 0;
	for(int i = 0; i < x.size(); ++i)
		dist += fabs(x[i] - y[i]);
	return dist;
}

float distL2Squared(vector<float>& x, vector<float>& y)
{
	float dist = 0;
	for(int i = 0; i < x.size(); ++i)
++		dist += (x[i] - y[i]) * (x[i] - y[i]);
	return dist;
}

float distHellinger(vector<float>& x, vector<float>& y)
{
	float dist = 0;
	for(int i = 0; i < x.size(); ++i)
	{
		float diff = sqrt(x[i]) - sqrt(y[i]);
		dist += diff * diff;
	}
	return dist;
}

float distKullbackLiebler(vector<float>& x, vector<float>& y)
{	
	float dist = 0;
	for(int i = 0; i < x.size(); ++i)
		dist += x[i] * log(x[i] / y[i]);
	return dist;
}


float distWithMetric(vector<float>& x, vector<float>& y, pDistanceMetric metric)
{
	switch( metric )
	{
	case DIST_L1:							return distL1(x, y);
	case DIST_L2:							return distL2Squared(x, y);
	case DIST_HELLINGER:					return distHellinger(x, y);
	case DIST_KULLBACK_LIEBLER:				return distKullbackLiebler(x, y);
	default:								assert(false);
	}
}