//
// EEMsearch.h
//
#ifndef EEMsearch_H
#define EEMsearch_H

//#include <vector>
#include <string>
#include "../utility/MyMap.h"

using std::string;
using std::MyMap;

namespace eem{
class ExpressionModuleSet;

class EEMsearch{
public:
	EEMsearch();
	virtual ~EEMsearch();

	virtual void setPvalue1Cutoff(double d) =0;
	virtual void suppressPvalue1Cutoff() =0;

	////////////////virtual void perform(){ ;}//=0;

	virtual ExpressionModuleSet* getExpressionModuleSet() =0;
	virtual void printResults(const string& outfile) throw () =0;
	virtual void printLog(const string& outfile) throw () =0;
	virtual MyMap<string, double> getPvalues() =0;

	//////////////virtual void setGeneSets(const MyMap<string, vector<string> >& geneSets){ ;}//=0;
	//////////////virtual void setMaxGeneSetSize(int i){ ;}//=0;
	//////////////virtual void setMinGeneSetSize(int i){ ;}//=0;

	virtual void setItrForPvalue2Calculation(int i) =0;
	virtual void recycleNullDistribution() =0;
	virtual string getLog() =0;
};
} // namespace eem
#endif // EEMsearch_H
