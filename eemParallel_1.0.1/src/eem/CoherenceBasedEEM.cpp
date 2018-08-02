//
// CoherenceBasedEEM.cpp
//
#include "CoherenceBasedEEM.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <map>

#include "../utility/MyException.h"
#include "../utility/MyFunc.h"

#include "CoherenceBasedEEMsearch.h"

using namespace std;
using namespace utility;
using namespace eem;

CoherenceBasedEEM::CoherenceBasedEEM()
{
}
CoherenceBasedEEM::~CoherenceBasedEEM()
{
}

CoherenceBasedEEM::CoherenceBasedEEM(CoherenceBasedEEMsearch *parent, const vector<string>& gene)
{


	this->parent = parent;
	absoluteRadius = parent->absoluteRadius;
	relativeRadius = parent->relativeRadius;

	itrForPvalue2Calculation =  parent->itrForPvalue2Calculation;
	seedGenes = gene;

	centerGene = "";
	Pvalue1 = -1.0;
	Pvalue2 = -1.0;
	Pvalue = -1.0;
	coreGeneSize = parent->coreGeneSize;

	if(!parent->nullDistrbutionForRecycle.empty()){
		if(parent->nullDistrbutionForRecycle.find(seedGenes.size()) != parent->nullDistrbutionForRecycle.end() ){
			nullDist = parent->nullDistrbutionForRecycle[seedGenes.size()];
		}else{
			nullDist.clear();// = new ArrayList<Integer>();
			parent->nullDistrbutionForRecycle[seedGenes.size()] = nullDist;
		}
	}else{
		nullDist.clear();// = new ArrayList<Integer>();
	}
}

string CoherenceBasedEEM::toString() const
{
	stringstream ss;
	ss << setprecision(17);

	vector<string> tmp;
	tmp.push_back("CoherenceBasedEEM");

	ss << "absoluteRadius=" << absoluteRadius;
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "relativeRadius=" << relativeRadius;
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "itrForPvalue2Calculation=" << itrForPvalue2Calculation;
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "coreGeneSize=" << coreGeneSize;
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "seedGenes=[";
	for(size_t i=0; i < seedGenes.size(); i++) {
		if (i > 0) ss << ", ";
		ss << seedGenes[i];
	}
	ss << "]";
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "moduleGenes=[";
	for(size_t i=0; i < moduleGenes.size(); i++) {
		if (i > 0) ss << ", ";
		ss << moduleGenes[i];
	}
	ss << "]";
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "center=[";
	for(size_t i=0; i < center.size(); i++) {
		if (i > 0) ss << ", ";
		ss << center[i];
	}
	ss << "]";
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "Pvalue=" << Pvalue;
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "Pvalue1=" << Pvalue1;
	tmp.push_back(ss.str());

	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "Pvalue2=" << Pvalue2;
	tmp.push_back(ss.str());

	return MyFunc::join(" ", tmp);
}

void CoherenceBasedEEM::findCenter()
{
	int i,j;
	moduleGenes.clear();// = new ArrayList<String>();
	centerGene = "";
	double distSum = 0;

	for(i =0; i < seedGenes.size(); i++){
		vector<string> tmp;
		double distSumTmp = 0.0;
		for(j =0; j < seedGenes.size(); j++){
			if(i == j){
				tmp.push_back(seedGenes[j]);
				continue;
			}
			double tmp2 = parent->distfunc->get(seedGenes[i], seedGenes[j]);
			if(tmp2 < absoluteRadius){
				tmp.push_back(seedGenes[j]);
				distSumTmp += tmp2;
			}
		}
		if(moduleGenes.empty() || tmp.size() > moduleGenes.size() || (tmp.size() ==  moduleGenes.size() && distSumTmp < distSum) ){
			moduleGenes  = tmp;
			center = parent->Exp.getRow(seedGenes[i]);
			centerGene = seedGenes[i];
			distSum = distSumTmp;
		}
	}
	MyMap<string, double> dist2center;
	for(i =0; i < moduleGenes.size(); i++){
		dist2center[moduleGenes[i]] = parent->distfunc->get(moduleGenes[i], centerGene);
	}
	MyFunc::sortKeysByAscendingOrderOfValues(dist2center, moduleGenes);
}

void CoherenceBasedEEM::refineCenter() throw ()
{
	int i;

	vector<string> coreGenes ;
	map<double, int> centerSum2count;// = new HashMap<Double, Integer>();

	//initialize center and core gene
	if(coreGeneSize < moduleGenes.size()){
		//  coreGenes = moduleGenes.subList(0, coreGeneSize-1);
		coreGenes = vector<string>(moduleGenes.begin(), moduleGenes.begin()+coreGeneSize-1);
	}else{
		coreGenes = moduleGenes;
	}
	center =  parent->Exp.getRowMeans(coreGenes);
	centerSum2count[MyFunc::sum(center)]= 1;

	//update center and core gene;
	vector<double> bestCenter;     // = new ArrayList<Double>();
	vector<string> bestModuleGenes;// = new ArrayList<String>();
	for(int k = 0; k < 10; k++){

		MyMap<string, double> dist2center;//  = new HashMap<String, Double>();
		for(i =0; i < seedGenes.size(); i++){
			double tmp = parent->distfunc->get(center, parent->Exp.getRow(seedGenes[i]));
			if(tmp < absoluteRadius){
				dist2center[seedGenes[i]] = tmp;
			}
		}
		MyFunc::sortKeysByAscendingOrderOfValues(dist2center, moduleGenes);
		if(moduleGenes .size() > coreGeneSize){
			//  coreGenes =  moduleGenes.subList(0, coreGeneSize-1);
			coreGenes = vector<string>(moduleGenes.begin(), moduleGenes.begin()+coreGeneSize-1);
		}else{
			coreGenes = moduleGenes;
		}
		center =  parent->Exp.getRowMeans(coreGenes);
		double sum = MyFunc::sum(center);
		if(centerSum2count.find(sum) != centerSum2count.end()){
			return;
		}else{
			centerSum2count[sum] = 1;
		}
		if(bestModuleGenes.size() < moduleGenes.size()){
			bestModuleGenes = moduleGenes;
			bestCenter = center;
		}
	}
	moduleGenes = bestModuleGenes;
	center = bestCenter;
}


void CoherenceBasedEEM::refineCenter0( ) throw ()
{
	int i,j;
	if(moduleGenes.size() < 3){
		cerr << "unnable to find an enough number of module genes!" << endl;
	}

	vector<string> coreGenes;
	if(coreGeneSize < moduleGenes.size()){
		MyMap<string, double> dist2center;

		for(i =0; i < moduleGenes.size(); i++){
			dist2center[moduleGenes[i]] = parent->distfunc->get(moduleGenes[i], centerGene);
		}
		MyFunc::sortKeysByAscendingOrderOfValues(dist2center, coreGenes);
		//   coreGenes = coreGenes.subList(0, coreGeneSize-1);
		vector<string>::iterator it_start = coreGenes.begin() + coreGeneSize-1;
		coreGenes.erase(it_start, coreGenes.end());
	}else{
		coreGenes = moduleGenes;
	}

	vector<vector<string> > CoreTriplet = MyFunc::getAllCombination(coreGenes, 3);
	double distSum = 0;

	for(i=0; i<CoreTriplet.size();i++){
		double distSumTmp = 0;
		vector<string> tmp;
		vector<double> tripletCenter = parent->Exp.getRowMeans(CoreTriplet[i]);
		for(j =0; j < seedGenes.size(); j++){
			double tmp2 = parent->distfunc->get( tripletCenter, parent->Exp.getRow(seedGenes[j]) );
			if(tmp2 < absoluteRadius){
				tmp.push_back(seedGenes[j]);
				distSumTmp += tmp2;
			}
		}
		if(tmp.size() > moduleGenes.size() || (tmp.size() ==  moduleGenes.size() && distSumTmp < distSum)){
			moduleGenes  = tmp;
			center = tripletCenter;
			distSum = distSumTmp;
		}
	}
}

void CoherenceBasedEEM::findModuleGenes() throw ()
{
	findCenter();
	try {
		refineCenter();
	}
	catch (MyException& err) {
		cerr << err.what_message("CoherenceBasedEEM::findModuleGenes") << endl;
	}
}

vector<string> CoherenceBasedEEM::getGlobalCoherentGenes()
{
	vector<string> tmp;// = new ArrayList<String>();

	for(int j =0; j < parent->allGenes.size(); j++){
		double tmp2 = parent->distfunc->get( center, parent->Exp.getRow(parent->allGenes[j]) );
		if(tmp2 < absoluteRadius){
			tmp.push_back(parent->allGenes[j]);
		}
	}
	return tmp;
}

void CoherenceBasedEEM::generateNullDist(int n)
{
	int j= 0;
	while(nullDist.size()<n){
		vector<string> sample;
		MyFunc::sample(parent->allGenes, seedGenes.size(), sample);
		CoherenceBasedEEM e( parent, sample );
		try{
			e.findModuleGenes();
			nullDist.push_back(e.getGlobalCoherentGenes().size());
		}
		catch (MyException& err) {
			//err.printStackTrace();//Java
			j++;
			if(j >= 10){
				cerr << err.what_message("Cannot generate null distribution!") << endl;
				throw MyException();
			}
		}
	}
}

void CoherenceBasedEEM::calculatePvalue1()
{
	if(moduleGenes.empty()){
		Pvalue1 = 0.0;
	}else{
		double p = MyFunc::calculatePvalueForSetOverlap(parent->allGenes.size(), seedGenes.size(), (int)(parent->allGenes.size()*relativeRadius), moduleGenes.size() );
		Pvalue1 = -log10( p );
	}
	Pvalue = Pvalue1;
}

vector<int> CoherenceBasedEEM::get_densityDist()
{
	return densityDist;
}

MyMap<int, double> CoherenceBasedEEM::getDensityDist()
{
	MyMap<int, double> DistMap;// = new HashMap<Integer, Double>();

	vector<int>::iterator it;
	//for(Integer I : nullDist){
	for(it=nullDist.begin(); it != nullDist.end(); it++){
		int I = *it;
		if(DistMap.find(I) != DistMap.end()){
			DistMap[I] = DistMap[I] + 1.0;
		}else{
			DistMap[I] = 1.0;
		}
	}
	return  DistMap;
}

void CoherenceBasedEEM::calculatePvalue2()
{
	generateNullDist(itrForPvalue2Calculation);

	MyMap<int, double> mDist =  getDensityDist();
	double P = 0.0;

	MyMap<int, double>::iterator it;
	for(it=mDist.begin(); it != mDist.end(); it++){
		int I = it->first;
		if(moduleGenes.size() > I){
			continue;
		}
		double p = MyFunc::calculatePvalueForSetOverlap( parent->allGenes.size(), seedGenes.size(), I, moduleGenes.size() );
		P += p * mDist[I];
	}
	Pvalue2 = -log10( P/(double)itrForPvalue2Calculation );
	Pvalue = Pvalue2;
}
//
// processing the same as constructor. CoherenceBssedEEM(parent, gene). 
//
void CoherenceBasedEEM::setParentNullDist(const vector<int>& null_dist)
{
	parent->nullDistrbutionForRecycle[(int)seedGenes.size()] = null_dist;
}


CoherenceBasedEEM::CoherenceBasedEEM(int irank, MPI_Comm comm)
{
	double tmp[7];
	MPI_Status status;
	MPI_Recv(tmp, 7, MPI_DOUBLE, irank, 1101, comm, &status);
	absoluteRadius = tmp[0]; // double
	relativeRadius = tmp[1]; // double
	itrForPvalue2Calculation = (int) tmp[2]; // int
	coreGeneSize = (int) tmp[3]; // int
	Pvalue = tmp[4]; // double
	Pvalue1 = tmp[5]; // double
	Pvalue2 = tmp[6]; // double

	MyFunc::vstrRecv(seedGenes, irank, comm); // vector<string>
	MyFunc::vstrRecv(moduleGenes, irank, comm); // vector<string>
	MyFunc::vdblRecv(center, irank, comm); // vector<double>
}
void CoherenceBasedEEM::send_to(int irank, MPI_Comm comm)
{
	double tmp[7];
	tmp[0] = absoluteRadius; // double
	tmp[1] = relativeRadius; // double
	tmp[2] = itrForPvalue2Calculation; // int
	tmp[3] = coreGeneSize; // int
	tmp[4] = Pvalue; // double
	tmp[5] = Pvalue1; // double
	tmp[6] = Pvalue2; // double
	MPI_Send(tmp, 7, MPI_DOUBLE, irank, 1101, comm);

	MyFunc::vstrSend(seedGenes, irank, comm); // vector<string>
	MyFunc::vstrSend(moduleGenes, irank, comm); // vector<string>
	MyFunc::vdblSend(center, irank, comm); // vector<double>
}
