//
// MyMap.h
//
#ifndef MyMap_H
#define MyMap_H

//namespace utility {

#ifdef HAVE_UNORDERED_MAP

#include <unordered_map>
#define MyMap unordered_map

#else

#include <map>
#define MyMap map

#endif

//} // namespace utility

#endif // MyMap_H
