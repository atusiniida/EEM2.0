//
// ExpressionModule.cpp
//
#include "ExpressionModule.h"

#include <ostream>
#include <cfloat>
#include <cmath>

#include "../utility/MyFunc.h"
#include "../utility/MyMat.h"

#include "ExpressionModuleSet.h"
#include "EEM.h"
#include "CoherenceBasedEEM.h"

using namespace std;
using namespace utility;
using namespace eem;

ExpressionModule::~ExpressionModule()
{
}

ExpressionModule::ExpressionModule(const ExpressionModule& e)
{
	id = e.id;
	eem = e.eem;
	activityProfile = e.activityProfile;
	moduleGenes = e.moduleGenes;
	seedGenes = e.seedGenes;
	expressionModuleCluster = e.expressionModuleCluster;
	minimunPvalue = e.minimunPvalue;
	PvalueCorrectedForMultipleTest = e.PvalueCorrectedForMultipleTest;
	numberOfMultipleTesting = e.numberOfMultipleTesting;
}

ExpressionModule::ExpressionModule(const string& id, EEM *eem, const MyMat& Exp)
{
	this->id = id;
	this->eem = eem;
	minimunPvalue = eem->getPvalue();
	PvalueCorrectedForMultipleTest = eem->getPvalue();
	numberOfMultipleTesting = 1;
	seedGenes = eem->getSeedGenes();
	moduleGenes = eem->getModuleGenes();
}


ExpressionModule::ExpressionModule(int irank, MPI_Comm comm)
{
	MyFunc::strRecv(id, irank, comm); // string
	eem = new CoherenceBasedEEM(irank, comm); // EEM*
	MyFunc::mstrdblRecv(activityProfile, irank, comm); // MyMap<string, double>
	MyFunc::vstrRecv(moduleGenes, irank, comm); // vector<string>
	MyFunc::vstrRecv(seedGenes, irank, comm); // vector<string>

	//recv(expressionModuleCluster); // ExpressionModuleSet* ; IMPOSSIBLE
	expressionModuleCluster = NULL;

	double tmp[3];
	MPI_Status status;
	MPI_Recv(tmp, 3, MPI_DOUBLE, irank, 2345, comm, &status);
	minimunPvalue = tmp[0];
	PvalueCorrectedForMultipleTest = tmp[1];
	numberOfMultipleTesting = (int) tmp[2];
}

void ExpressionModule::send_to(int irank, MPI_Comm comm)
{
	MyFunc::strSend(id, irank, comm); // string
	eem->send_to(irank, comm); // EEM*
	MyFunc::mstrdblSend(activityProfile, irank, comm); // MyMap<string, double>
	MyFunc::vstrSend(moduleGenes, irank, comm); // vector<string>
	MyFunc::vstrSend(seedGenes, irank, comm); // vector<string>

	//send(expressionModuleCluster); // ExpressionModuleSet* ; IMPOSSIBLE

	double tmp[3];
	tmp[0] = minimunPvalue;
	tmp[1] = PvalueCorrectedForMultipleTest;
	tmp[2] = numberOfMultipleTesting;
	MPI_Send(tmp, 3, MPI_DOUBLE, irank, 2345, comm);
}


string ExpressionModule::getId()
{
	return id;
}

double ExpressionModule::getPvalue()
{
	return PvalueCorrectedForMultipleTest;
}

const EEM* ExpressionModule::getEEM() const
{
	return eem;
}

string ExpressionModule::toString() const
{
	string str = id + "\t" + MyFunc::toString(PvalueCorrectedForMultipleTest) +
		"\t" + MyFunc::toString(numberOfMultipleTesting) +
		"\t" + getEEM()->toString();

	return str;
}

//--
// generate ExpressionModule
//--
ExpressionModule* ExpressionModule::getBetterExpressionModules(const ExpressionModule& e1, const ExpressionModule& e2)
{
	ExpressionModule* e3;
	if(e1.minimunPvalue > e2.minimunPvalue){
		e3 = new ExpressionModule(e1);//----------- generate ExpressionModule
		e3->minimunPvalue = e1.minimunPvalue;
	}else{
		e3 = new ExpressionModule(e2);//----------- generate ExpressionModule
		e3->minimunPvalue = e2.minimunPvalue;
	}
	e3->numberOfMultipleTesting = e1.numberOfMultipleTesting + e2.numberOfMultipleTesting;
	if(e3->minimunPvalue != DBL_MAX){
		double p = e3->minimunPvalue;
		if(p < 5){
			p = pow(10, -p);
			p = 1- pow(1-p, e3->numberOfMultipleTesting);
			p  = - log10(p);
		}else{
			p = p - log10( (double)e3->numberOfMultipleTesting );
		}
		e3->PvalueCorrectedForMultipleTest = p;
	}
	return e3;
		
}


namespace eem {
ostream &operator <<(ostream &stream, const ExpressionModule &em)
{
	stream << em.toString();
	return stream;
}
} // namespace eem
