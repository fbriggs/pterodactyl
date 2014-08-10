// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#pragma once

#include "pVectorND.h"
#include <float.h>
#include <iostream>
#include <vector>

using namespace std;

// pMatrix defines a m x n matrix with overloaded operators and basic operations such as transpose,
// determinant and inverse. examples:
//
// stream output and transpose
// ---------------------------
//
// cout << pMatrix(1,2,3,4).transpose() << endl;
//
// special case constructor for 3x3
// ---------------------------------
//
// pMatrix A(	1, 2, 0,
//				4, 5, 6,
//				7, 8, 9);
//
// solve Ax = b
// ------------
//
// pVectorND b(1,2,3);
// pVectorND xHat = (((A.transpose() * A).inverse()) * A) * b;
// cout << xHat << endl;
//
class pMatrix
{
public:
	// stores the data in the matrix. it is ok to modify the contents, but don't mess with the dimensions
	// note that it is also possible to acces data with like so: pMatrix A; cout <<  A(1,2);
	vector< vector<float> > data_;
	
	// these are public to avoid the annoyance of getters, but should be treated as read-only
	int rows_, cols_;
	
	// empty constructor, makes a null matrix.
	pMatrix() { rows_ = 0; cols_ = 0; }
	
	// zero constructor, takes matrix size
	pMatrix(int rows, int cols)
	{
		zeros(rows, cols);
	}
	
	// special case constructor for 2x2
	pMatrix(	float a, float b, 
				float c, float d)
	{
		zeros(2,2);
		data_[0][0] = a;	data_[0][1] = b;
		data_[1][0] = c;	data_[1][1] = d;
	}
	
	// special case constructor for 3x3
	pMatrix(	float m00, float m01, float m02,
				float m10, float m11, float m12,
				float m20, float m21, float m22 )
	{
		zeros(3,3);
		data_[0][0] = m00;	data_[0][1] = m01;	data_[0][2] = m02;
		data_[1][0] = m10;	data_[1][1] = m11;	data_[1][2] = m12;
		data_[2][0] = m20;	data_[2][1] = m21;	data_[2][2] = m22;
	}

	// copy constructor
	pMatrix(const pMatrix &c) { copy(c); }

	// overwrite the contents of the current matrix and copy in all data from
	// matrix c, including its dimensions if they differ from the current dimensions
	void copy(const pMatrix &c)
	{
		data_.clear();
		
		rows_ = c.rows_;
		cols_ = c.cols_;
		
		for(int i = 0; i < rows_; ++i)
		{
			vector<float> row;
			for(int j = 0; j < cols_; ++j)
				row.push_back(c.data_[i][j]);
			data_.push_back(row);
		}
	}
	
	// fill the matrix with zeros
	const pMatrix& zeros(int rows, int cols)
	{
		data_.clear();
		
		for(int i = 0; i < rows; ++i)
		{
			vector<float> row;
			for(int j = 0; j < cols; ++j)
				row.push_back(0);
			data_.push_back(row);
		}
		
		rows_ = rows;
		cols_ = cols;
		
		return *this;
	}
	
	// makes an n x n identity matrix
	const pMatrix& identity(int n)
	{
		zeros(n, n);
		for(int i = 0; i < n; ++i) data_[i][i] = 1;
		
		return *this;
	}
	
	// assignment operator
	const pMatrix& operator=(const pMatrix &rhs) { copy(rhs); return *this; }

	// operator +=
	const pMatrix& operator+=(const pMatrix &rhs)
	{
		assert(rows_ == rhs.rows_ && cols_ == rhs.cols_ && "matrix += operator failed, dimension mismatch");
		
		for(int i = 0; i < rows_; ++i)
			for(int j = 0; j < cols_; ++j)
				data_[i][j] += rhs.data_[i][j];
		
		return *this;
	}
	
	// operator -=
	const pMatrix& operator-=(const pMatrix &rhs)
	{
		assert(rows_ == rhs.rows_ && cols_ == rhs.cols_ && "matrix += operator failed, dimension mismatch");
		
		for(int i = 0; i < rows_; ++i)
			for(int j = 0; j < cols_; ++j)
				data_[i][j] -= rhs.data_[i][j];
		
		return *this;
	}
	
	// operator *= (scalar)
	const pMatrix& operator*=(float scale)
	{
		for(int i = 0; i < rows_; ++i)
			for(int j = 0; j < cols_; ++j)
				data_[i][j] *= scale;
		
		return *this;
	}
	
	// operator /= (scalar)
	const pMatrix& operator/=(float scale)
	{
		assert(scale != 0 && "matrix /= operator failed, divide by zero");
		for(int i = 0; i < rows_; ++i)
			for(int j = 0; j < cols_; ++j)
				data_[i][j] *= scale;
		
		return *this;
	}

	// 2d subscript operators for shorter access to data_
	float& operator()(unsigned int i, unsigned int j)		{ return data_[i][j]; }
	float operator()(unsigned int i, unsigned int j) const	{ return data_[i][j]; }

	// operator ==
	bool operator==(const pMatrix &rhs)
	{
		if( rows_ != rhs.rows_ || cols_ != rhs.cols_ ) return false;
		
		for(int i = 0; i < rows_; ++i)
			for(int j = 0; j < cols_; ++j)
				if( data_[i][j] != rhs.data_[i][j] ) 
					return false;
		
		return true;
	}
	
	// operator !=
	bool operator!=(const pMatrix &rhs) { return !(*this == rhs); }
	
	// addition, reusing code from +=
	pMatrix operator+(const pMatrix &rhs)	{ pMatrix A(*this); A += rhs; return A; }
	
	// substraction
	pMatrix operator-(const pMatrix &rhs)	{ pMatrix A(*this); A -= rhs; return A; }
	
	// negatation
	pMatrix operator-()						{ pMatrix A(*this); A *= -1; return A; }
	
	// standard matrix multiplication
	pMatrix operator*(const pMatrix &rhs)
	{
		assert( cols_ == rhs.rows_ && "matrix multiplication dimension mismatch" );
		
		pMatrix A(rows_, rhs.cols_);
		
		for(int i = 0; i < A.rows_; ++i)
			for(int j = 0; j < A.cols_; ++j)
				for(int k = 0; k < cols_; ++k)
					A.data_[i][j] += data_[i][k] * rhs.data_[k][j];
		
		return A;
	}
	
	// friend operators for scalar multiplication and division
	friend pMatrix operator*(const pMatrix &m, float scale)		{ pMatrix A(m); A *= scale; return A; }
	friend pMatrix operator*(float scale, const pMatrix &m)		{ pMatrix A(m); A *= scale; return A; }
	friend pMatrix operator/(const pMatrix &m, float scale)		{ pMatrix A(m); A /= scale; return A; }
	
	// multiply a matrix by a vector
	friend pVectorND operator*(const pMatrix &m, const pVectorND &v)
	{
		assert(m.cols_ == v.size() && "dim mismatch multiplying vector by matrix");
		
		pVectorND r(m.rows_);
		
		for(int i = 0; i < m.rows_; ++i)
			for(int j = 0; j < m.cols_; ++j)
				r[i] += v[j] * m.data_[i][j];
		
		return r;
	}
	
	// compute the determinant of this matrix
	float determinant() const
	{
		// special cases for 2x2 and 3x3 matrices
		if( rows_ == 2 && cols_ == 2 ) return	data_[0][0]*data_[1][1] - data_[0][1]*data_[1][0];
		if( rows_ == 3 && cols_ == 3 ) return	data_[0][0]*data_[1][1]*data_[2][2] +
												data_[0][1]*data_[1][2]*data_[2][0] +
												data_[0][2]*data_[1][0]*data_[2][1] -
												data_[0][0]*data_[1][2]*data_[2][1] -
												data_[0][1]*data_[1][0]*data_[2][2] -
												data_[0][2]*data_[1][1]*data_[2][0];
	
		assert(false && "determinant of n x n not implemented");
	}
	
	float trace() const
	{
		assert(rows_ == cols_);
		
		float sum = 0;
		for(int i = 0; i < rows_; ++i)
			sum += data_[i][i];
		return sum;
	}
	
	// find the inverse of the matrix
	pMatrix inverse() const
	{
		assert( rows_ == cols_ && "attempt to take inverse of non-square matrix" );
	
		// special case for 2x2
		if( rows_ == 2 && cols_ == 2)
		{
			float det = determinant();
			return pMatrix(	data_[1][1] / det, -data_[0][1] / det, -data_[1][0] / det, data_[0][0] / det);
		}
		
		// for larger matrices, we will use the method of putting the original matrix next to the identity matrix,
		// then performing row operations until the original matrix is the identity, and what started as the identity is the inverse
		
		pMatrix I; I.identity(cols_); 
		pMatrix A(*this); // copy this matrix
		
		int i, j, k;
		
        // start from first column to the next 
        for(i = 0; i < rows_; i++) 
		{
            float div = A.data_[i][i]; 
            if( fabs(div) < 1e-10 ) // pivotal element too small
			{  
				for(j = i + 1; j < rows_; j++)
					if( fabs(div = A.data_[j][i]) > 1e-10 )
					{
                        A.swapRows(i, j);
                        I.swapRows(i, j);
                        break ;
					}
				if( j == rows_ ) 
				{ 
					cout << *this << endl;
					assert(false && "inverted singular matrix");
				}
            }
			
            for(j = 0; j < cols_; j++) // divide entire row by pivotal element 
			{	
				A.data_[i][j] /= div;
				I.data_[i][j] /= div;
            }
			
            for(k = 0; k < rows_; k++)
			{ 
                if( k == i ) continue;
                float mult = A.data_[k][i];
                if( fabs(mult) < 1e-10 ) continue;
                for(j = 0; j < cols_; j++) 
				{
					A.data_[k][j] -= mult * A.data_[i][j];
					I.data_[k][j] -= mult * I.data_[i][j];                
				}
            }
        }
		
        return I;
	}
	
	// returns the transpose of this matrix
	pMatrix transpose() const
	{
		pMatrix A(cols_, rows_);
		
		for(int i = 0; i < A.rows_; ++i)
			for(int j = 0; j < A.cols_; ++j)
				A.data_[i][j] = data_[j][i];
		
		return A;
	}
	
	// swamp rows i and j
	const pMatrix& swapRows(int i, int j)
	{
		float temp;
		for(int p = 0; p < cols_; ++p)
		{
			temp = data_[i][p];
			data_[i][p] = data_[j][p];
			data_[j][p] = temp;
		}
		
		return *this;
	}
	
	// stream outputput operator. prints in the format matlab reads
	friend std::ostream& operator<<(std::ostream& stream, const pMatrix &m)
	{
		stream << "[";
		for(int i = 0; i < m.rows_; ++i)
		{
			for(int j = 0; j < m.cols_; ++j)
				stream << m.data_[i][j] << (j == m.cols_ - 1 ? ";" : ",");
			if( i != m.rows_ - 1 ) stream << endl;
		}
		stream << "]" << endl;
		return stream;
	}
	
	// normalize the values in the matrix to sum to 1
	pMatrix& normalizeSumTo1()
	{
		// normalize g to sum to 1
		float sum = 0;
		for(int i = 0; i < rows_; ++i)
			for(int j = 0; j < cols_; ++j)
				sum += data_[i][j];
		
		for(int i = 0; i < rows_; ++i)
			for(int j = 0; j < cols_; ++j)
				data_[i][j] /= sum;
		
		return *this;
	}
	
	// normalize all values to the range 0-1
	void normalizeTo01()
	{
		float minVal = FLT_MAX;
		for(int i = 0; i < rows_; ++i)
			for(int j = 0; j < cols_; ++j)
				minVal = min(minVal, data_[i][j]);
		for(int i = 0; i < rows_; ++i)
			for(int j = 0; j < cols_; ++j)
				data_[i][j] -= minVal;
		
		float maxVal = -FLT_MAX;
		for(int i = 0; i < rows_; ++i)
			for(int j = 0; j < cols_; ++j)
				maxVal = max(maxVal, data_[i][j]);
		
		assert(maxVal != 0);
		
		for(int i = 0; i < rows_; ++i)
			for(int j = 0; j < cols_; ++j)
				data_[i][j] /= maxVal;
	}	
};