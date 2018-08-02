//
// StopWatch.h
//
#ifndef StopWatch_H
#define StopWatch_H

#include <string>
#include "Date.h"

namespace utility {

class StopWatch{
private:
	Date startDate;
	Date stopDate;
	int miliSecond;
	int second;
	int minute;
	int hour;
	double elapsedTime;
public:
	StopWatch();
	virtual ~StopWatch(){}
	void stop();
	std::string toString() const;
	const Date& getStartDate() const { return startDate;}
	const Date& getStopDate() const { return stopDate;}
	double getElapsedTime() const { return elapsedTime; }
};

} // namespace utility

#endif // StopWatch_H
