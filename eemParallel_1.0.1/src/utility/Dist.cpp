//
// Dist.cpp
//
#include "Dist.h"

#include <cmath>
#include <sstream>

#include "MyFunc.h"
#include "MyException.h"

using namespace std;
using namespace utility;

Dist::Dist()
	: diagonalElement(0)
{
	diagonalElement= 0;
	b_argCon = false;
}
/*type must be 'e'(euclideanDist), 'c' (pearsonCorrelation),
* or 'C' (pearsonCorrelationForNormalizedList)
*/
Dist::Dist(const MyMat& m, char type, MPI_Comm comm)
{
	//Dist(m, type);
	call_argumentConstruct(m, type, comm);
}

Dist::~Dist()
{
}
//----
//
//----
const vector<string> &Dist::getNames() const
{
	return name;
}

void Dist::setDiagonalElement(double d)
{
	diagonalElement = d;	
}

double Dist::get(const string& s, const string& t)
{
	int i = name2index[s];
	int j = name2index[t];
	if(i > j){
		return M[i][j];
	}
	if(j > i){
		return M[j][i];
	}
	return diagonalElement;
}


//--
// operator 2013.11.16
//--
Dist& Dist::operator =(const Dist& dist)
{
	M = dist.M;
	name2index = dist.name2index;	
	n = dist.n;
	name = dist.name;

	diagonalElement = dist.diagonalElement;

	b_argCon = dist.b_argCon;//------------- add 2013.11.20

	return *this;
}

//--
// add method
//--
void Dist::call_argumentConstruct(const MyMat& m, char type, MPI_Comm comm)
{
	diagonalElement= 0;
	b_argCon = true;

	n = m.rowSize();
	name = m.getRowNames();
	int i;
	for(i =0; i< n; i++){
		name2index[name[i]] = i;
	}
	M.resize(n);
	for(i=0; i < n; i++){
		M[i].resize(i);
	}

	switch(type){
		/*case('c'):
			{
			MyMat copy_m = m;
			copy_m.normalizeRows();
			for(i = 0; i < n; i++){
				for(j = 0; j < i; j++){
					M[i][j] = MyFunc::pearsonCorrelationForNormalizedList(copy_m.getRow(i),copy_m.getRow(j));
			    }
			}
			setDiagonalElement(1.0);
			}
			break;*/

		case('C'):
			calculateM(m, comm);
			break;
		/*case('e'):
			for(i = 0; i < n; i++){
				for(j = 0; j < i; j++){
					M[i][j] = MyFunc::euclideanDist(m.getRow(i),m.getRow(j));
			    }
			}
			break;*/
		default:
			//cerr << "type must be 'c','C', or 'e'" << endl;
			cerr << "type must be 'C'"<< endl;
			break;
	}
}

void Dist::calculateM(const MyMat& m, MPI_Comm comm)
{
	int i, j;
	int numProcs, myRank;
	double stctime,edctime;//calculate time
	double stdtime,eddtime;//distribute time

	MPI_Comm_size(comm, &numProcs);
	MPI_Comm_rank(comm, &myRank);

	vector<int> minRow, maxRow, widthRow;
	//--
	// divide M
	//--
	MyFunc::minmaxRow( minRow, maxRow, widthRow, n, numProcs);

	// parameter check.
	stringstream ss;
	ss << "Dist::myRank=" << myRank << ", M row size=" << n
	   << ", minRow=" << minRow[myRank] << ", maxRow=" << maxRow[myRank] << ", widthRow=" << widthRow[myRank] << "\n";
	cerr << ss.str();

	MPI_Barrier(comm);
	stctime = MPI_Wtime();
	//--
	// parallelized
	//--
	for(i = minRow[myRank]; i < maxRow[myRank]; i++){// single => for(i = 0; i < n; i++)
		for(j = 0; j < i; j++){
			M[i][j] = MyFunc::pearsonCorrelationForNormalizedList(m.getRow(i),m.getRow(j));
		}
	}
	edctime = MPI_Wtime();
	MPI_Barrier(comm);
	stdtime = MPI_Wtime();
	//--
	// result distribute
	//--
	vector<double> vdbl(n*(n-1)/2, 0.0);
	for(i = minRow[myRank]; i < maxRow[myRank]; i++) {
		int i0 = i*(i-1)/2;
		for(j = 1; j < i; j++) {
			vdbl[i0 + j] = M[i][j];
		}
	}
	double *tmp = new double[vdbl.size()];
	MPI_Allreduce(&vdbl[0], tmp, vdbl.size(), MPI_DOUBLE, MPI_SUM, comm);
	for (i = 0; i < vdbl.size(); i++) {
		vdbl[i] = tmp[i];
	}
	delete[] tmp;
	vector<double>::iterator it = vdbl.begin();
	for(i = 0; i < n; i++) {
		for(j = 0; j < i; j++) {
			M[i][j] = *it;
			it++;
		}
	}
	eddtime = MPI_Wtime();
	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "Dist::calculateM: proc=" << myRank << "/" << numProcs
	   << ", total time:" << eddtime - stctime
	   << ", calculate time:" << edctime - stctime
	   << ", wait time:" << stdtime - edctime
	   << ", distribute time:" << eddtime - stdtime << "\n";
	cerr << ss.str();

	setDiagonalElement(1);
}

void Dist::clear()
{
	size_t nRow = M.size();
	for(size_t i=0; i < nRow; i++) M[i].clear();
	M.clear();

	name2index.clear();
	n = 0;
	name.clear();

	diagonalElement = 0.0;

	b_argCon=false;
}
