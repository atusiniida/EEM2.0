//
// AbstractGeneSetAnalysis.cpp
//
#include "AbstractGeneSetAnalysis.h"

#include <iostream>
#include "../utility/MyFunc.h"

using namespace std;
using namespace utility;
using namespace eem;

AbstractGeneSetAnalysis::AbstractGeneSetAnalysis()
{
	maxSeedGeneSize = 2000;
	minSeedGeneSize = 10 ;
}
AbstractGeneSetAnalysis::~AbstractGeneSetAnalysis()
{
}
//--
//
//--
void AbstractGeneSetAnalysis::setGeneSets(const MyMap<string, vector<string> >& geneSets)
{
	MyMap<string, vector<string> >::const_iterator it;

	for(it=geneSets.begin(); it != geneSets.end(); it++){
		vector<string> tmp;
		MyFunc::isect( it->second, allGenes, tmp );

		if(tmp.size() < minSeedGeneSize || tmp.size() > maxSeedGeneSize){
			std::cerr << it->first << ": seed geneset size (" << tmp.size() << ") is out of range!" << endl;
			continue;
		}
		//
		// map => multimap ?
		//
		seedGeneSets[it->first] = tmp;// Java: seedGeneSets.put(it->first, tmp);
	}
}
void AbstractGeneSetAnalysis::setMaxGeneSetSize(int i)
{
	maxSeedGeneSize = i;
}
void AbstractGeneSetAnalysis::setMinGeneSetSize(int i)
{
	minSeedGeneSize = i;
}
void AbstractGeneSetAnalysis::perform(MPI_Comm comm)
{
	cerr << "Unsupported Operation  Exception : AbstractGeneSetAnalysis::perform()" << endl;
}
