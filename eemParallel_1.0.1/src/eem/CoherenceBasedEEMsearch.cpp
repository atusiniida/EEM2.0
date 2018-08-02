//
// CoherenceBasedEEMsearch.cpp
//
#include "CoherenceBasedEEMsearch.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>

#include "../utility/MyException.h"

#include "DistConverter.h"
#include "CoherenceBasedEEM.h"

using namespace std;
using namespace utility;
using namespace eem;

CoherenceBasedEEMsearch::CoherenceBasedEEMsearch()
{
	distfunc = NULL;
	distConverter = NULL;
}
CoherenceBasedEEMsearch::~CoherenceBasedEEMsearch()
{
	if(distfunc) delete distfunc;
	if(distConverter) delete distConverter;
}


CoherenceBasedEEMsearch::CoherenceBasedEEMsearch(const MyMat& E)
{
	//super();

	itrForPvalue2Calculation = 300;
	Pvalue1Cutoff = -log10(0.05);
	absoluteRadius = 0.0;
	relativeRadius = 0.0;
	coreGeneSize = 10;
	originalExp = E;
	Exp = E;
	Exp.normalizeRows();

	//Cor = null;
	Cor.clear();//----- 2013.11.20: Cor clear

	allGenes = Exp.getRowNames();

	distfunc = new Correlation(&Cor);//---------- 2013.11.21: argment "&Cor"
	distConverter =  new DistConverter(this);
};

//
// generate CoherenceBasedEEM
//
void CoherenceBasedEEMsearch::setEEM(MPI_Comm comm)
{
	if( seedGeneSets.empty() ){
		cerr << "seedGeneSets must be set!" << endl;
	}
	if(!Cor.use_argumentConstruct()){
		calculateCor(comm);
	}
	if(absoluteRadius == 0.0){
		setRelativeRadius(0.05, comm);
	}

	MyMap<string, vector<string> >::iterator it;
	for(it=seedGeneSets.begin(); it != seedGeneSets.end(); it++){
		//
		// generate EEM
		//
		EEM* e =  new CoherenceBasedEEM(this, it->second);
		eems[it->first] = e;
		seeds.push_back(it->first);
	}
}

void CoherenceBasedEEMsearch::setCoreGeneSize(int i)
{
	if( i > 3 ){
		coreGeneSize = i;
	}
}

void CoherenceBasedEEMsearch::setCor(const Dist& D)
{
	vector<string>::const_iterator it = Exp.getRowNames().begin();
	for(; it != Exp.getRowNames().end(); it++){
		if( find(D.getNames().begin(), D.getNames().end(), *it) == D.getNames().end() ){
			cerr << "Correlation data must include data for all genes in expression data" << endl;
		}
	}
	Cor = D;
}

// parallelized
void CoherenceBasedEEMsearch::calculateCor(MPI_Comm comm)
{
	cerr << "Calculating correlations..." << endl;
	Cor.call_argumentConstruct(Exp, 'C', comm);
}

void CoherenceBasedEEMsearch::setAbsoluteRadius(double r, MPI_Comm comm)
{
	if(!Cor.use_argumentConstruct()){
		calculateCor(comm);
	}
	absoluteRadius = r;
	relativeRadius = distConverter->convertAbsolute2relativeDist(r);
}

void CoherenceBasedEEMsearch::setRelativeRadius(double r, MPI_Comm comm)
{
	if(!Cor.use_argumentConstruct()){
		calculateCor(comm);
	}
	relativeRadius = r;
	absoluteRadius = distConverter->convertRelative2absoluteDist(r);
}

string CoherenceBasedEEMsearch::toString()
{
	vector<string> tmp;

	tmp.push_back("CoherenceBasedEEM");

	stringstream ss;
	ss << "coreGeneSize=" << coreGeneSize;
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "absoluteRadius= " << absoluteRadius;
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "relativeRadius= " << relativeRadius;
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "maxGeneSetSize= " << maxSeedGeneSize;
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "minGeneSetSize= " << minSeedGeneSize;
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "itrForPvalue2Calculation= " << itrForPvalue2Calculation;
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "Pvalue1Cutoff= " << Pvalue1Cutoff;
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "seeds=" << seeds.size();
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "candidates=" << candidates.size();
	tmp.push_back(ss.str());

	return MyFunc::join(" ", tmp);
}
