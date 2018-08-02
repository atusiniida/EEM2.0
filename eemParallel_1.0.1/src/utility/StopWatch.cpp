//
// StopWatch.cpp
//
#include "StopWatch.h"

#include <cstdio>

using namespace utility;

#ifdef VC
#define snprintf _snprintf
#endif

StopWatch::StopWatch()
	: miliSecond(0), second(0), minute(0), hour(0), elapsedTime(0.0)
{
	startDate.getTime();
	// stopDate
}

void StopWatch::stop()
{
	stopDate.getTime();
	elapsedTime = stopDate - startDate;
	second = (int) elapsedTime;
	miliSecond = (int)(elapsedTime*1.0e+3) - second;
	if(second < 60){
		minute = 0;
		hour = 0;
		return;
	}
	minute = second/60;
	second = second%60;
	if(minute < 60){
		hour = 0;
		return;
	}
	hour = minute/60;
	minute = minute%60;
	return;
}

#define BUFSIZE 256

std::string StopWatch::toString() const
{
	static char buffer[BUFSIZE];
	if(! stopDate.isTimeSet()){
		buffer[0] = '\0';
		return std::string(buffer);
	}
	if(hour != 0){
		snprintf(buffer, sizeof(buffer), "%dhr %dmin %dsec", hour, minute, second);
		return std::string(buffer);
	}
	if(minute !=0 ){
		snprintf(buffer, sizeof(buffer), "%dmin %dsec", minute, second);
		return std::string(buffer);
	}
	snprintf(buffer, sizeof(buffer), "%d.%dsec", second, miliSecond);
	return std::string(buffer);
}
