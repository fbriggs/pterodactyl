// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#pragma once

#include <math.h>
#include <iostream>

class pComplex
{
public:	
	double real_, imag_;

	// zero constructor
	pComplex()
	{
		real_ = 0;
		imag_ = 0;
	}
	
	// initialization constructor
	pComplex(double real, double imag)
	{
		real_ = real;
		imag_ = imag;
	}
	
	 // copy constructor
	pComplex(const pComplex &copy)
	{
		real_ = copy.real_;
		imag_ = copy.imag_;
	}

	const pComplex& operator=(const pComplex &rhs) //assignment
	{
		real_ = rhs.real_;
		imag_ = rhs.imag_;
		
		return *this;
	}

	const pComplex& operator+=(const pComplex &rhs)
	{
		real_ += rhs.real_;
		imag_ += rhs.imag_;
		
		return *this;
	}
	
	const pComplex& operator-=(const pComplex &rhs)
	{
		real_ -= rhs.real_;
		imag_ -= rhs.imag_;
		
		return *this;
	}
	
	const pComplex& operator*=(const pComplex &rhs)
	{
		double a = real_;
		double b = imag_;
		double c = rhs.real_;
		double d = rhs.imag_;
		
		double newR = a * c - b * d;
		double newI = b * c + a * d;
		
		real_= newR;
		imag_ = newI;
		
		return *this;
	}
		
	const pComplex& operator*=(double scale)
	{
		real_ *= scale;
		imag_ *= scale;
		
		return *this;
	}
	
	
	const pComplex& operator/=(double scale)
	{
		real_ /= scale;
		imag_ /= scale;
		
		return *this;
	}
	
	bool operator==(const pComplex &rhs) const
	{
		return real_ == rhs.real_ && imag_ == rhs.imag_;
	}
	
	bool operator!=(const pComplex &rhs) const
	{
		return real_ != rhs.real_ || imag_ != rhs.imag_;
	}
	
	pComplex operator+(const pComplex &rhs) const
	{
		return pComplex(real_ + rhs.real_, imag_ + rhs.imag_);
	}

	pComplex operator-(const pComplex &rhs) const
	{
		return pComplex(real_ - rhs.real_, imag_ - rhs.imag_);
	}
	
	pComplex operator+(const double r) const
	{
		return pComplex(real_ + r, imag_);
	}
	
	pComplex operator-(const double r) const
	{
		return pComplex(real_ - r, imag_);
	}
	
	pComplex operator-() const
	{
		return pComplex(-real_, -imag_);
	}
	
	friend pComplex operator*(const pComplex &v, double scale)
	{
		return pComplex(v.real_*scale, v.imag_*scale);
	}
	
	friend pComplex operator*(const pComplex &u, const pComplex &v)
	{
		double a = u.real_;
		double b = u.imag_;
		double c = v.real_;
		double d = v.imag_;
		
		return pComplex(a * c - b * d, b * c + a * d);
	}
	
	
	friend pComplex operator/(const pComplex &u, const pComplex &v)
	{
		double a = u.real_;
		double b = u.imag_;
		double c = v.real_;
		double d = v.imag_;
		
		return pComplex(a * c - b * d, b * c + a * d);
	}
	
	friend pComplex operator*(double scale, const pComplex &v)
	{
		return pComplex(v.real_*scale, v.imag_*scale);
	}

	friend pComplex operator/(const pComplex &v, double scale)
	{
		return pComplex(v.real_/scale, v.imag_/scale);
	}

	double mag() const
	{
		return sqrt(real_*real_ + imag_*imag_);
	}
	
	double magSquared() const
	{
		return real_*real_ + imag_*imag_;
	}
	
	double phase() const
	{
		double theta = atan2(imag_, real_);
		if( theta < 0 ) theta += 2 * M_PI;
		return theta;
	}
	
	const pComplex& normalize()
	{
		double l = mag();
		real_ /= l;
		imag_ /= l;
		
		return *this;
	}
	
	friend std::ostream& operator<<(std::ostream& out, const pComplex &v)
	{
		out << "[" << v.real_ << ", " << v.imag_ << "]";
		return out;
	}
};