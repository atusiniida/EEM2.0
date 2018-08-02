///
///EEM.h
///

#ifndef EEM_H
#define EEM_H

#include "../utility/MyMap.h"
#include <vector>
#include <string>
#include <mpi.h>

using std::MyMap;
using std::vector;
using std::string;

namespace eem{

///
///the EEM class is  an abstract  class of CoherenceBasedEEM
///

class EEM{
public:
	EEM();
	virtual ~EEM();

	virtual const vector<string>& getSeedGenes() const =0;
	virtual const vector<string>& getModuleGenes() const =0;
	virtual double getPvalue1() =0;
	virtual double getPvalue2() =0;
	virtual double getPvalue() =0;
	virtual void setPvalue2(double p) =0;
	virtual void setPvalue1(double p) =0;
	virtual void setPvalue(double p) =0;
	virtual void findModuleGenes() throw () =0;
	virtual void calculatePvalue1() =0;
	virtual void calculatePvalue2() =0;
	virtual void cutParent() =0;

	//
	// add interface method
	//
public:
	virtual string toString() const =0;
protected:
	virtual void findCenter() =0;
	virtual void refineCenter() throw () =0;
	virtual void refineCenter0()throw () =0;

	virtual vector<string> getGlobalCoherentGenes() =0;
	virtual void generateNullDist(int n) =0;

	virtual MyMap<int, double> getDensityDist() =0;

	// add method2 '14.3
public:
	virtual vector<int> get_densityDist() =0;
	virtual void set_densityDist(const vector<int>& dendist) =0;

public:
	virtual double  getAbsoluteRadius() =0;
	virtual double  getRelativeRadius() =0;
	virtual const vector<double>&  getCenter() const =0;
	virtual const string& getCenterGene() const =0;
	virtual void   setAbsoluteRadius(double r) =0;

	//
	// add method
	//
public:
	virtual string name(){ return "EEM";}
	//--
	// add method '2014.3
	//--
	virtual void setModuleGenes(const vector<string>& module_genes) =0;
	virtual void setCenter(const vector<double>& c) =0;
	virtual void setCenterGene(const string& center_gene) =0;

	virtual void setNullDist(const vector<int>& vint) =0;
	virtual vector<int> getNullDist() =0;

	virtual void setParentNullDist(const vector<int>& null_dist) =0;


	EEM(int irank, MPI_Comm comm);
	virtual void send_to(int irank, MPI_Comm comm)=0;
};
} // namespace eem
#endif // EEM_H
