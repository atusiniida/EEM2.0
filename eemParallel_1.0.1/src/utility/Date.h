//
// Date.h
//
#ifndef Date_H
#define Date_H

#ifdef VC
#include <sys/timeb.h>
#else
#include <sys/time.h>
#endif

#include <string>

namespace utility {

class Date {
#ifdef VC
	_timeb tv;
#else
	struct timeval tv;
#endif

	bool is_time_set;

public:
	Date();
	virtual ~Date();
	int getTime();
	bool isTimeSet() const;
	double operator-(const Date &start) const;
	std::string toString() const;
};

} // namespace utility

#endif // Date_H
