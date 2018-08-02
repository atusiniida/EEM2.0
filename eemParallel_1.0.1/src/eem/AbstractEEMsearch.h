//
// AbstractEEMsearch.h
//
#ifndef AbstractEEMsearch_H
#define AbstractEEMsearch_H

#include "../utility/MyMap.h"
#include <vector>
#include <string>
#include <mpi.h>

using std::vector;
using std::string;
using std::MyMap;

#include "../utility/MyMat.h"
#include "../utility/StopWatch.h"

using utility::MyMat;
using utility::StopWatch;

#include "AbstractGeneSetAnalysis.h"
#include "EEMsearch.h"

namespace eem{
class EEM;
class ExpressionModuleSet;
///
/// inheritance relationships are:
/// CoherenceBasedEEMsearch->AbstractEEMsearch->AbstractGenesetAnalysis
///
class AbstractEEMsearch : public AbstractGeneSetAnalysis, public EEMsearch{

public:
	AbstractEEMsearch();
	AbstractEEMsearch(const AbstractEEMsearch& A);
	virtual ~AbstractEEMsearch();

	//--
	//
	//--
	MyMat  originalExp;
	int    itrForPvalue2Calculation; // iteration number for calculation precise p-values (Pvalue2)
	double Pvalue1Cutoff; // cutoff for approximate p-values (Pvalue1), used to find candidates, which are subjected to Pvalue2 calculation
	StopWatch stopWatch;

	MyMap<string, EEM*> eems;	// store EEM objects for each geneset. keys are seed  gene set ID
	vector<string> seeds; //  seed (input) gene set IDs
	vector<string> candidates; // IDs of gene set subjected to Pvalue2 calculation
	string errLog;
	string timeLog;
	MyMap<int, vector<int> > nullDistrbutionForRecycle;


	virtual void setGeneSets(const MyMap<string, vector<string> >& geneSets);
	virtual void setPvalue1Cutoff(double d);
	virtual void suppressPvalue1Cutoff();
	virtual void setItrForPvalue2Calculation(int i);
	virtual void recycleNullDistribution();

protected:
	// these are used for EEM search and to be parallelized
	virtual void findModuleGenes(MPI_Comm comm);
	virtual void calculatePvalue1(MPI_Comm comm);
	virtual void findCandidates(); // find candidates, which are subjected to Pvalue2 calculation, based on Pvalue1Cutoff
	virtual void calculatePvalue2(MPI_Comm comm);

private:
	void calculatePvalue2_inner(int myRank, int i, vector<int>& candidates_flg);

public:
	virtual void printResults(const string& outfile) throw ();
	virtual void printLog(const string& outfile) throw ();

protected:
	virtual void setEEM(MPI_Comm comm);

public:
	//
	// generate ExpressionModuleSet, ExpressionModule
	//
	virtual ExpressionModuleSet* getExpressionModuleSet(); // get result of EEM search
	virtual string getLog();

protected:
	virtual void writeTimeLog();

public:
	virtual void perform(MPI_Comm comm);
	virtual MyMap<string, double> getPvalues();

	// add method
	string name(){ return "AbstractEEMsearch";}

};
} // namespace eem
#endif // AbstractEEMsearch_H
