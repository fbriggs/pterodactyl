// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once

#include <iostream>
#include <vector>
#include <math.h>
#include "pStringUtils.h"
using namespace std;

// pVectorND defines an n-dimensional vector with overloaded operators. this class is used in conjunction with pMatrix,
// e.g. to multiply a matrix and a vector, etc. examples:
//
// special case constructors for 2 or 3d
// -------------------------------------
// pVectorND(1,2)
//
// vector of zeros:
// ----------------
// pVectorND().zeros(5)
//
// concatenation:
// --------------
// pVectorND(1,2) ^ 3
// 
// combining lots of operators:
// ----------------------------
//
// cout << (pVectorND(1,2) ^ (pVectorND().zeros(5)) ^ 3 ^ 4) << endl;
//
class pVectorND
{
public:
	vector<float> data_;

	/// constructors ///
	
	// empty constructor
	pVectorND() {}
	
	// this constructor takes an int specifying the vector dimensional and makes it the 0 vector
	pVectorND(int dim) { zeros(dim); }
	
	// special case constructors for 2 and 3d vectors
	pVectorND(float x, float y)				{ data_.push_back(x); data_.push_back(y); }
	pVectorND(float x, float y, float z)	{ data_.push_back(x); data_.push_back(y); data_.push_back(z); }
	
	// copy constructor
	pVectorND(const pVectorND &copy)
	{
		for(int i = 0; i < copy.data_.size(); ++i) data_.push_back(copy.data_[i]);
	}
	
	// copy from a vector of floats
	pVectorND(const vector<float>& v)
	{
		for(int i = 0; i < v.size(); ++i) data_.push_back(v[i]);
	}
	
	/// vector operations ///
	
	// return the number of elements in the vector
	int size() const { return data_.size(); }
	
	// fill the vector with n zeros (erases whatever was there before)
	const pVectorND&  zeros(int n) { data_.clear(); for(int i = 0; i < n; ++i) data_.push_back(0); return *this;}
	
	// append elements to vector
	void push_back(float f) { data_.push_back(f); }
	
	// remove all elements from the vector
	void clear() { data_.clear(); }
	
	// length (norm)
	float length() const { return sqrt(lengthSquared()); }
	
	// length squared
	float lengthSquared() const
	{
		float sum = 0;
		for(int i = 0; i < data_.size(); ++i) sum += data_[i] * data_[i];
		return sum;
	}
	
	// make this vector's length 1
	const pVectorND& normalize()
	{
		float l = length();
		assert(l != 0);
		for(int i = 0; i < data_.size(); ++i) data_[i] /= l;
		return *this;
	}
	
	// dot product
	float dot(const pVectorND &rhs) const
	{
		float sum = 0; 
		for(int i = 0; i < data_.size(); ++i) sum += data_[i] * rhs.data_[i];
		return sum;
	}
	
	// projection
	pVectorND projectOnto(const pVectorND &v) const
	{
		return v * (dot(v) / v.lengthSquared());
	}
	
	// convert to string
	string toString(int precision = -1) const
	{
		string s  = "<";
		for(int i = 0 ; i < data_.size(); ++i)
		{
			s += pStringUtils::floatToString(data_[i], precision);
			if( i != data_.size() - 1 ) s += ", ";
		}
		s += ">";
		return s;
	}
	
	/// operators - all basic arithmetic, plus streams, append and element access ///
	
	// stream operator (ex: cout << pVectorND(1,2,3) << endl; )
	friend std::ostream& operator<<(std::ostream& out, const pVectorND &v)
	{
		out << v.toString();
		return out;
	}
	
	// ^ operator - appends the element to this vector
	pVectorND operator^(const float &e) const
	{
		pVectorND r(data_);
		r.data_.push_back(e);
		return r;
	}
	
	// ^ operator - concatenate two vectors into one
	pVectorND operator^(const pVectorND& v) const
	{
		pVectorND r(data_);
		for(int i = 0; i < v.data_.size(); ++i) r.data_.push_back(v.data_[i]);
		return r;
	}
	
	// [] operators access data_ as shorthand
	float& operator[](unsigned int i) { return data_[i]; }
	const float operator[](unsigned int i) const { return data_[i]; }
	
	// assignment operator
	const pVectorND& operator=(const pVectorND &rhs)
	{
		data_.clear();
		for(int i = 0; i < rhs.data_.size(); ++i) data_.push_back(rhs.data_[i]);
		return *this;
	}

	// += operator
	const pVectorND& operator+=(const pVectorND &rhs)
	{
		for(int i = 0; i < data_.size(); ++i) data_[i] += rhs.data_[i];
		return *this;
	}
	
	// -= operator
	const pVectorND& operator-=(const pVectorND &rhs)
	{
		for(int i = 0; i < data_.size(); ++i) data_[i] -= rhs.data_[i];
		return *this;
	}
	
	// *= operator
	const pVectorND& operator*=(float scale)
	{
		for(int i = 0; i < data_.size(); ++i) data_[i] *= scale;
		return *this;
	}
	
	// /= operator
	const pVectorND& operator/=(float scale)
	{
		assert( scale != 0 );
		for(int i = 0; i < data_.size(); ++i) data_[i] /= scale;
		return *this;
	}

	// == equality operator
	bool operator==(const pVectorND &rhs) const 
	{
		for(int i = 0; i < data_.size(); ++i) 
			if( data_[i] != rhs.data_[i] ) return false;
		return true;
	}
	
	// != inequality operator
	bool operator!=(const pVectorND &rhs) const
	{
		return !(*this == rhs);
	}
	
	// + operator
	pVectorND operator+(const pVectorND &rhs) const
	{
		pVectorND r(data_);
		for(int i = 0; i < data_.size(); ++i) r.data_[i] += rhs.data_[i];
		return r;
	}
	
	// - operator
	pVectorND operator-(const pVectorND &rhs) const
	{
		pVectorND r(data_);
		for(int i = 0; i < data_.size(); ++i) r.data_[i] -= rhs.data_[i];
		return r;
	}
	
	// negation
	pVectorND operator-() const
	{
		pVectorND r(data_);
		for(int i = 0; i < data_.size(); ++i) r.data_[i] *= -1;
		return r;
	}
	
	// scalar multiplication
	friend pVectorND operator*(const pVectorND &v, float scale)
	{
		pVectorND r(v.data_);
		for(int i = 0; i < v.data_.size(); ++i) r.data_[i] *=  scale;
		return r;
	}
	
	// more scalar multiplication
	friend pVectorND operator*(float scale, const pVectorND &v)
	{
		return v * scale; // use the operator defined above
	}
	
	// scalar division
	friend pVectorND operator/(const pVectorND &v, float scale)
	{
		assert(scale != 0);
		pVectorND r(v.data_);
		for(int i = 0; i < v.data_.size(); ++i) r.data_[i] /= scale;
		return r;
	}
};