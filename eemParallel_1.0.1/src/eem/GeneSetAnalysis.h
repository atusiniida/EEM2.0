//
// GeneSetAnalysis.h
//
#ifndef GeneSetAnalysis_H
#define GeneSetAnalysis_H

#include "../utility/MyMap.h"
#include <string>
#include <vector>
#include <mpi.h>

using std::MyMap;
using std::string;
using std::vector;

namespace eem{

class GeneSetAnalysis{
public:
	GeneSetAnalysis();
	virtual ~GeneSetAnalysis();

	virtual void perform(MPI_Comm comm) =0;
	virtual void setGeneSets(const MyMap<string, vector<string> >& geneSets) =0;
	virtual void setMaxGeneSetSize(int i) =0;
	virtual void setMinGeneSetSize(int i) =0;
};
} // namespace eem
#endif // GeneSetAnalysis_H
