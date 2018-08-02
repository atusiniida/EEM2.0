//
// CoherenceBasedEEM.h
//
#ifndef CoherenceBasedEEM_H
#define CoherenceBasedEEM_H

#include <vector>
#include <string>
#include "../utility/MyMap.h"

using std::vector;
using std::string;
using std::MyMap;

#include "EEM.h"

namespace eem{
//
//public class CoherenceBasedEEM implements EEM, Serializable{
//
class CoherenceBasedEEMsearch;

///
///an  CoherenceBasedEEM  class object performs an EEM seach for each geneset
///
///
    
class CoherenceBasedEEM : public EEM {
//protected:
public:
	CoherenceBasedEEM();
	CoherenceBasedEEM(CoherenceBasedEEMsearch *parent, const vector<string>& gene);
	virtual ~CoherenceBasedEEM();

protected:
	double absoluteRadius;
	double relativeRadius;
	vector<string> seedGenes;
	vector<string> moduleGenes;

	double Pvalue1; // approximate p-value (minus log scale)
	double Pvalue2; // precise p-value (minus log scale)
	double Pvalue;  // Pvalue2 if available, or Pvalue1

	int itrForPvalue2Calculation;
	int coreGeneSize;
	vector<double> center;

	//--
	// transient(Java)
	//--
	vector<int> nullDist;
	vector<int> densityDist;

	CoherenceBasedEEMsearch *parent;

	string centerGene;

public:
	virtual string toString() const;
protected:
	virtual void findCenter();
	virtual void refineCenter() throw ();
	virtual void refineCenter0()throw ();

public:
	virtual void findModuleGenes() throw ();

protected:
	virtual vector<string> getGlobalCoherentGenes();
	virtual void generateNullDist(int n);

public:
	virtual void calculatePvalue1();
	virtual void calculatePvalue2();

protected:
	virtual MyMap<int, double> getDensityDist();

	// add method2 '14.3
public:
	virtual vector<int> get_densityDist();
	virtual void set_densityDist(const vector<int>& dendist){ densityDist= dendist;}

public:
	//--
	// inline
	//--
	virtual const vector<string>& getModuleGenes() const {return moduleGenes;}
	virtual void  cutParent(){parent = NULL;}

	virtual double  getAbsoluteRadius(){return absoluteRadius;};
	virtual double  getRelativeRadius(){return relativeRadius;};
	virtual const vector<double>&  getCenter() const {return center;};
	virtual const string& getCenterGene() const {return centerGene;};
	virtual void   setAbsoluteRadius(double r){absoluteRadius = r;}

	virtual double  getPvalue(){return Pvalue;}
	virtual double  getPvalue1(){return Pvalue1;}
	virtual double  getPvalue2(){return Pvalue2;}
	virtual void  setPvalue2(double p){Pvalue2 = p;}
	virtual void  setPvalue1(double p){Pvalue1 = p;}
	virtual void  setPvalue(double p){Pvalue = p;}
	virtual const vector<string>& getSeedGenes() const {return seedGenes;}
	//--
	// add method '2014.3
	//--
	virtual void setModuleGenes(const vector<string>& module_genes){ moduleGenes=module_genes; }
	virtual void setCenter(const vector<double>& c){                 center=c; }
	virtual void setCenterGene(const string& center_gene){           centerGene=center_gene; }

	virtual void setNullDist(const vector<int>& vint){ nullDist = vint;}
	virtual vector<int> getNullDist(){ return nullDist;}

	virtual void setParentNullDist(const vector<int>& null_dist);

public:
	// add method
	virtual string name(){ return "CoherenceBasedEEM";}


	CoherenceBasedEEM(int irank, MPI_Comm comm);
	virtual void send_to(int irank, MPI_Comm comm);
};
} // namecpace eem
#endif // CoherenceBasedEEM_H
