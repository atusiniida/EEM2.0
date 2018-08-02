//
// CoherenceBasedEEMsearch.h
//
#ifndef CoherenceBasedEEMsearch_H
#define CoherenceBasedEEMsearch_H

#include <vector>
#include <string>

using std::vector;
using std::string;

#include "../utility/Dist.h"
#include "../utility/MyFunc.h"
#include "../utility/MyMat.h"

using utility::Dist;
using utility::MyFunc;
using utility::MyMat;

#include "AbstractEEMsearch.h"

namespace eem{
class DistConverter;
class CoherenceBssedEEM;
    
///
/// an CoherenceBasedEEMsearch class objects performs EEM analysis
///
/// inheritance relationships are:
/// CoherenceBasedEEMsearch->AbstractEEMsearch->AbstractGenesetAnalysis
///
class CoherenceBasedEEMsearch : public AbstractEEMsearch{
public:
	CoherenceBasedEEMsearch();
	virtual ~CoherenceBasedEEMsearch();


	MyMat Exp; // expression matrix
	Dist  Cor; // correlation half matrix

	// these are parameters for EEM search
	double absoluteRadius;
	double relativeRadius; // the parameter to be parallelized
	int coreGeneSize;

	DistConverter *distConverter;  // use for parameter conversion (relativeRadius <-> absoluteRadius)

	//--
	// inner class
	//--
	class DistFunc {
	public:
		DistFunc(){}
		DistFunc(Dist* cor){ Cor = cor;}
		virtual ~DistFunc(){}

		Dist* Cor;//------ 2013.11.21 pointer

		virtual double get(const vector<double>& a, const vector<double>& b){ return 1;}//=0;
		virtual double get(const string& a, const string& b){ return 1;}//=0;
	};
	class Correlation: public DistFunc {
		public:
			Correlation(){}
			Correlation(Dist* cor){ Cor = cor;}
			virtual ~Correlation(){}

			virtual double get(const vector<double>& a, const vector<double>& b){
				return 1.0 - MyFunc::pearsonCorrelationForNormalizedList(a,b);
			}
			virtual double get(const string& a, const string& b){
				return 1.0 - Cor->get(a,b);
			}
	};

	//-----------------
	DistFunc *distfunc;
	//-----------------

	CoherenceBasedEEMsearch(const MyMat& E);
	
protected:
	// prepare for EEM search
	// here a correlation matrix is calculated (to be parallelized) and parameter conversion performed
	virtual void setEEM(MPI_Comm comm);// generate CoherenceBasedEEM

public:
	virtual void setCoreGeneSize(int i);

	virtual void setCor(const Dist& D);
	// a correlation matrix is calculated (to be parallelized)
	virtual void calculateCor(MPI_Comm comm);

	// parameter setting and conversion (absoluteRadius <-> relativeRadius)
	virtual void setAbsoluteRadius(double r, MPI_Comm comm);
	virtual void setRelativeRadius(double r, MPI_Comm comm);

	virtual string toString();

};
} // namespace eem
#endif // CoherenceBasedEEMsearch_H
