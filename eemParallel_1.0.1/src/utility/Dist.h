//
// Dist.h
//
#ifndef Dist_H
#define Dist_H

#include <vector>
#include <string>
#include "MyMap.h"
#include <mpi.h>

using std::vector;
using std::string;
using std::MyMap;

#include "MyMat.h"

namespace utility{

class Dist{
public:
	Dist();
	/*type must be 'e'(euclideanDist), 'c' (pearsonCorrelation),
	 * or 'C' (pearsonCorrelationForNormalizedList)
	 */
	Dist(const MyMat& m, char type, MPI_Comm comm);

	virtual ~Dist();

private:
	static long long serialVersionUID;// = -1890276595674423636L;

	vector<vector<double> >    M;
	MyMap<string, int>  name2index;
	int n;
	vector<string> name;

	double diagonalElement;

	bool b_argCon;//use Dist(MyMat&, char) ?

public:
	const vector<string>& getNames() const;
	void setDiagonalElement(double d);
	double get(const string& s, const string& t);

	//static Dist readFromBinary(const string& infile) throw ();// FileNotFoundException, IOException, ClassNotFoundException

	//--
	// operator 2013.11.16
	//--
	Dist& operator=(const Dist& dist);

	//--
	// add method 2013.11.16
	//--
	bool use_argumentConstruct(){ return b_argCon;}
	void call_argumentConstruct(const MyMat& m, char type, MPI_Comm comm);
	void clear();//------------2013.11.20

private:
	void calculateM(const MyMat& m, MPI_Comm comm);
};
} // namespace utility
#endif // Dist_H
