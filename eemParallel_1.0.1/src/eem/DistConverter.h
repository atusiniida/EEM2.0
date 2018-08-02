//
// DistConverter.h
//
#ifndef DistConverter_H
#define DistConverter_H

#include "../utility/MyMap.h"
#include <string>

using std::MyMap;
using std::string;

namespace eem{
class CoherenceBasedEEMsearch;
class CoherenceBasedEEM;
class DistConverter{
public:
	DistConverter();
	DistConverter(CoherenceBasedEEMsearch* parent);
	virtual ~DistConverter();

private:
	int maxItrForBisection;
	int maxItrForBisection2;
	double deltaForBisection;
	double upperAbsoluteDist;
	double lowerAbsoluteDist;
	CoherenceBasedEEMsearch *parent;
	CoherenceBasedEEM       *e;
	MyMap<string, double> absolute2relative;

public:
	double convertAbsolute2relativeDist(double absoluteDist);// generate CoherenceBasedEEM
	double convertRelative2absoluteDist(double relativeDist);

};
} // namespace eem
#endif // DistConverter_H
