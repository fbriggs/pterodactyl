// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#pragma once

/*
Example use:
 
pTuple<int, bool, string, char> tup = make_tuple(3, true, string("hello"), 'x');
int a; bool b; string c; char d;
extract(tup, a,b,c,d);
cout << a << " " << tup.third_ << endl; 
*/

struct nil_type {};

template<class T1, class T2 = nil_type, class T3 = nil_type, class T4 = nil_type> 
class pTuple
{
public:
	T1 first_;
	T2 second_;
	T3 third_;
	T4 fourth_;

	pTuple(T1& e1)							: first_(e1) {}
	pTuple(T1& e1, T2& e2)					: first_(e1), second_(e2) {}
	pTuple(T1& e1, T2& e2, T3& e3)			: first_(e1), second_(e2), third_(e3) {}
	pTuple(T1& e1, T2& e2, T3& e3, T4& e4)	: first_(e1), second_(e2), third_(e3), fourth_(e4) {}
};

template<class T1>									pTuple<T1>					make_tuple(T1 e1)						{ return pTuple<T1>(e1); }
template<class T1, class T2>						pTuple<T1, T2>				make_tuple(T1 e1, T2 e2)				{ return pTuple<T1, T2>(e1, e2); }
template<class T1, class T2, class T3>				pTuple<T1, T2, T3>			make_tuple(T1 e1, T2 e2, T3 e3)			{ return pTuple<T1, T2, T3>(e1, e2, e3); }
template<class T1, class T2, class T3, class T4>	pTuple<T1, T2, T3, T4>		make_tuple(T1 e1, T2 e2, T3 e3, T4 e4)	{ return pTuple<T1, T2, T3, T4>(e1, e2, e3, e4); }

template<class T1>
void extract(pTuple<T1> tup, T1& t1)
{
	t1 = tup.first_;
}

template<class T1, class T2>
void extract(pTuple<T1, T2> tup, T1& t1, T2& t2)
{
	t1 = tup.first_;
	t2 = tup.second_;
}

template<class T1, class T2, class T3>
void extract(pTuple<T1, T2, T3> tup, T1& t1, T2& t2, T3& t3)
{
	t1 = tup.first_;
	t2 = tup.second_;
	t3 = tup.third_;
}

template<class T1, class T2, class T3, class T4>
void extract(pTuple<T1, T2, T3, T4> tup, T1& t1, T2& t2, T3& t3, T4& t4)
{
	t1 = tup.first_;
	t2 = tup.second_;
	t3 = tup.third_;
	t4 = tup.fourth_;
}
