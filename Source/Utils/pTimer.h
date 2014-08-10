// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#pragma once

#include "pTimer.h"

#include <sys/time.h>
#include <vector>
#include <map>

typedef unsigned int PTIMER_MARK;

static const PTIMER_MARK PTIMER_INVALID_MARK = 1 << 31;

class pTimer
{
public:

	static pTimer* instance();
	
	void start();
	
	double getElapsed();
	
	PTIMER_MARK setMark();
	
	void resetMark(PTIMER_MARK mark);
	
	double timeSinceMark(PTIMER_MARK mark);
	
	void wait(double milliseconds);

	void suspendMark(PTIMER_MARK mark);
	
	void resumeMark(PTIMER_MARK mark);
	
	void advanceMark(PTIMER_MARK mark, double milliseconds);
	
	void printTimeSinceMark(PTIMER_MARK mark);
private:

	static pTimer *instance_;
	
	double initialTime_;
	
	std::vector<double> userMarks_;
	std::map<PTIMER_MARK, double> suspendedMarks_;
	
	pTimer();
	
	pTimer(const pTimer& copy);
	pTimer& operator=(const pTimer& rhs);
	
	double currentTime();
};
