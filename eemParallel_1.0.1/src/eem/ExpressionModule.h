//
// ExpressionModule.h
//
#ifndef ExpressionModule_H
#define ExpressionModule_H

#include "../utility/MyMap.h"
#include <vector>
#include <string>
#include <ostream>
#include <mpi.h>

using std::MyMap;
using std::vector;
using std::string;
using std::ostream;

#include "../utility/MyMat.h"

using utility::MyMat;

//
// a ExpressionModule class object contains EEM search result for each gene set
//
namespace eem{
class ExpressionModuleSet;
class EEM;

class ExpressionModule{
private:
	ExpressionModule(){};
public:
	ExpressionModule(const ExpressionModule& e);
	ExpressionModule(const string& id, EEM *eem, const MyMat &Exp);
	virtual ~ExpressionModule();

	ExpressionModule(int irank, MPI_Comm comm);
	void send_to(int irank, MPI_Comm comm);

private:

	string id;
	EEM  *eem;
	MyMap<string, double> activityProfile; // activity module of expression modules
	vector<string>      seedGenes;
	vector<string>      moduleGenes;

	// for multiple hypothesis testing (used when EEM search is performed with multiple parameters)
	double minimunPvalue; //-log10 scale
	double PvalueCorrectedForMultipleTest;
	int    numberOfMultipleTesting;

	ExpressionModuleSet* expressionModuleCluster;// for ExpressionModuleClustering


public:
	string getId();
	double getPvalue();
	const EEM*   getEEM() const;

	string toString() const;

	///
	/// obtain better ExpressionModules based on p-values while correcting multiple hypothesis testing
	///
	static ExpressionModule* getBetterExpressionModules(const ExpressionModule& e1, const ExpressionModule& e2);
};

ostream &operator <<(ostream &stream, const ExpressionModule& em);

} // namespace eem

#endif // ExpressionModule_H
