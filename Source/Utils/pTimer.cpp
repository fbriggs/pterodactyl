// Copyright (C) 2013 Forrest Briggs
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#include "pTimer.h"
#include <iostream>
using namespace std;

pTimer* pTimer::instance_ = NULL;

pTimer* pTimer::instance()
{
	if(!instance_)
		instance_ = new pTimer();
	return instance_;
}

pTimer::pTimer()
{
	start();
}

void pTimer::start()
{
	initialTime_ = currentTime();
}

double pTimer::getElapsed()
{
	return currentTime() - initialTime_;
}

PTIMER_MARK pTimer::setMark()
{
	userMarks_.push_back(currentTime());
	return (PTIMER_MARK)userMarks_.size() - 1;
}

void pTimer::resetMark(PTIMER_MARK mark)
{
	std::map<PTIMER_MARK, double>::iterator markItr = suspendedMarks_.find(mark);
	
	if(markItr != suspendedMarks_.end())
		resumeMark(mark);
	userMarks_[mark] = currentTime();
}

double pTimer::timeSinceMark(PTIMER_MARK mark)
{
	std::map<PTIMER_MARK, double>::iterator markItr = suspendedMarks_.find(mark);
	
	if(markItr != suspendedMarks_.end())
		return markItr->second - userMarks_[mark];
	else
		return currentTime() - userMarks_[mark];
}

void pTimer::wait(double milliseconds)
{
	double start = currentTime();
	while(true) 
		if(currentTime() - start > milliseconds)
			break;
}

void pTimer::suspendMark(PTIMER_MARK mark)
{
	suspendedMarks_.insert(std::pair<PTIMER_MARK, double> (mark, currentTime()));
}

void pTimer::resumeMark(PTIMER_MARK mark)
{
	std::map<PTIMER_MARK, double>::iterator markItr = suspendedMarks_.find(mark);
	
	if(markItr != suspendedMarks_.end())
	{
		userMarks_[mark] += (currentTime() - markItr->second);
		suspendedMarks_.erase(markItr);
	}
}

void pTimer::advanceMark(PTIMER_MARK mark, double milliseconds)
{
	userMarks_[mark] += milliseconds;
}

double pTimer::currentTime()
{
	timeval time;
	
	gettimeofday(&time, NULL);
	
	return time.tv_sec * 1000.0 + time.tv_usec / 1000.0;
}

void pTimer::printTimeSinceMark(PTIMER_MARK mark)
{
	float seconds =  timeSinceMark(mark) / 1000;
	float minutes = seconds / 60;
	cout << seconds << " seconds = " << minutes << " minutes" << endl;
}