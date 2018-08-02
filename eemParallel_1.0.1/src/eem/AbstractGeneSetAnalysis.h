//
// AbstractGeneSetAnalysis.h
//
#ifndef AbstractGeneSetAnalysis_H
#define AbstractGeneSetAnalysis_H

#include <vector>
#include "../utility/MyMap.h"
#include <string>
#include <mpi.h>

using std::vector;
using std::MyMap;
using std::string;

#include "GeneSetAnalysis.h"

namespace eem{

class AbstractGeneSetAnalysis : public GeneSetAnalysis{
public:
	AbstractGeneSetAnalysis();
	virtual ~AbstractGeneSetAnalysis();

	vector<string> allGenes;

	//for filtering of input gene sets based on their size
	int maxSeedGeneSize;// = 2000;
	int minSeedGeneSize;// = 10 ;
	MyMap<string, vector<string> > seedGeneSets;

	//--
	//
	//--
	virtual void setGeneSets(const MyMap<string, vector<string> >& geneSets);
	virtual void setMaxGeneSetSize(int i);
	virtual void setMinGeneSetSize(int i);

	virtual void perform(MPI_Comm comm);

};
} // namespace eem
#endif // AbstractGeneSetAnalysis_H
