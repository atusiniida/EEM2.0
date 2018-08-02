//
// MyException.h
//
#ifndef MyException_H
#define MyException_H

#include <string>
#include <exception>

using std::string;
using std::exception;

namespace utility{

class MyException : public exception {
public:
	MyException();
	MyException(const string& msg);
	virtual ~MyException() throw();

private:
	static long long serialVersionUID;// = 1L;
	string errMsg;

public:
	const string& getMessage() const;
	string what_message(const string&);
};
} // namespace utility
#endif // MyException_H
