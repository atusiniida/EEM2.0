//
// DistConverter.cpp
//
#include "DistConverter.h"

#include <iostream>
#include <sstream>
#include <cmath>

#include "../utility/MyFunc.h"
#include "../utility/MyException.h"

#include "CoherenceBasedEEM.h"
#include "CoherenceBasedEEMsearch.h"

using namespace std;
using namespace utility;
using namespace eem;

DistConverter::DistConverter()
{
	parent = NULL;
	e = NULL;
}
DistConverter::~DistConverter()
{
	if(e) delete e;
}

DistConverter::DistConverter(CoherenceBasedEEMsearch *parent){

	this->parent = parent;

	maxItrForBisection = 30;
	maxItrForBisection2 = 5;
	deltaForBisection  = 0.05;
	upperAbsoluteDist  = 0.0;
	lowerAbsoluteDist  = 0.0;

	e = NULL;
}
//
// generate CoherenceBasedEEM
//
double DistConverter::convertAbsolute2relativeDist(double absoluteDist)
{
	string s = MyFunc::toString(absoluteDist);

	if(absolute2relative.find(s) != absolute2relative.end()){
		cerr << absoluteDist << " in absolute distance is equal to " << absolute2relative[s] << " in relative distance." << endl;
		return absolute2relative[s];
	}
	//
	//  generate   CoherenceBasedEEM
	//
	if(e) delete e;
	e = new CoherenceBasedEEM(parent, parent->allGenes);
	e->setAbsoluteRadius(absoluteDist);
	try{
		e->findModuleGenes();
	}
	catch (MyException& err) {
		cerr << err.what_message("Unable to convert distance!") << endl;
	}
	double relativeDist = ((double) e->getModuleGenes().size())/parent->allGenes.size();

	cerr << absoluteDist << " in absolute distance is equal to " << relativeDist << " in relative distance." << endl;

	absolute2relative[s] = relativeDist;

	return relativeDist;
}

double  DistConverter::convertRelative2absoluteDist(double  relativeDist)
{
	double upperAbsoluteDist;
	double lowerAbsoluteDist;

	if( this->upperAbsoluteDist != 0.0){
		upperAbsoluteDist  = this->upperAbsoluteDist;
	}else{
		upperAbsoluteDist = 2;
	}
	if(this->lowerAbsoluteDist != 0.0){
		lowerAbsoluteDist  = this->lowerAbsoluteDist;
	}else{
		lowerAbsoluteDist = 0;
	}
	cerr << "converting " << relativeDist <<  " in relative distance to absolute distance..." << endl;


	double upperRelativeDist = convertAbsolute2relativeDist(upperAbsoluteDist);
	if(upperRelativeDist < relativeDist){
		stringstream ss;
		ss << "Unable to convert distance! increase upperAbsoluteDist(upperRelativeDist=" << upperRelativeDist << " relativeDist=" << relativeDist << ")." << endl;
		cerr << ss.str() << endl;
	}
	double lowerRelativeDist = convertAbsolute2relativeDist(lowerAbsoluteDist);
	if(lowerRelativeDist > relativeDist){
		stringstream ss;
		ss << "Unable to convert distance! increase lowerAbsoluteDist(lowerRelativeDist=" << lowerRelativeDist << " relativeDist=" << relativeDist << ")." << endl;
		cerr << ss.str() << endl;
	}


	int i,j=0;
	double midAbsoluteDist =  0.0;
	for(i = 0; i < maxItrForBisection; i++){
		midAbsoluteDist = (upperAbsoluteDist + lowerAbsoluteDist) / 2;
		double midRelativeDist =  convertAbsolute2relativeDist(midAbsoluteDist);
		if(midRelativeDist > upperRelativeDist || midRelativeDist < lowerRelativeDist){
			j++;
			if(j == maxItrForBisection2){
				cerr << "Unable to convert distance! convertAbsolute2relativeDist might not be a monolonically increasing function." << endl;
			}
		}
		if( abs(midRelativeDist - relativeDist) / relativeDist < deltaForBisection){
			return midAbsoluteDist;
		}
		if(midRelativeDist < relativeDist){
			lowerAbsoluteDist = midAbsoluteDist;
			lowerRelativeDist = midRelativeDist;
		}else{
			upperAbsoluteDist = midAbsoluteDist;
			upperRelativeDist = midRelativeDist;
		}
	}
	cerr << "WARN: Unable to convert distance! Increase maxItrForBisection or deltaForBisection." << endl;

	return midAbsoluteDist;
}
