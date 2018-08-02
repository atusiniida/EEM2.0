//
// MyException.cpp
//
#include "MyException.h"

#include <iostream>

using namespace std;
using namespace utility;

MyException::MyException()
{
}
MyException::MyException(const string& msg)
{
	errMsg = msg;
}

MyException::~MyException() throw ()
{
}

const string& MyException::getMessage() const
{
	return errMsg;
}

string MyException::what_message(const string& msg)
{
	string str;

	str  = exception::what();
	str += " " + msg;

	return str;
}
