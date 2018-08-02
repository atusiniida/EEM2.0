//
// Date.cpp
///
#include "Date.h"

#ifdef VC
#include <sys/types.h>
#include <iostream>
#include <sstream>
#endif

#ifdef VC
static double timeval_to_sec(const _timeb &tv) {
	return (double)tv.time + (double)tv.millitm * 1.0e-3;//milli sec
}
#else
static double timeval_to_sec(const struct timeval &tv) {
	return (double)tv.tv_sec + (double)tv.tv_usec * 1.0e-6;//micro sec
}
#endif

namespace utility {

Date::Date()
	: is_time_set(false)
{
}

Date::~Date() {
}

int Date::getTime() {
#ifdef VC
	_ftime(&tv);// TODO change to safe
	is_time_set = true;

	return 0;//success=0, failure=-1
#else
	int iret = gettimeofday(&tv, 0);
	if(iret==0) is_time_set = true;

	return iret;
#endif
}

bool Date::isTimeSet() const {
	return is_time_set;
}

double Date::operator-(const Date &start) const {
	double t_start = timeval_to_sec(start.tv);
	double t_stop = timeval_to_sec(tv);
	return t_stop - t_start;
}

std::string Date::toString() const {
#ifdef VC
	std::stringstream ss;
	ss << tv.time;
	return ss.str();
#else
	return std::string(ctime(&(tv.tv_sec)));
#endif
}

} // namespace utility
