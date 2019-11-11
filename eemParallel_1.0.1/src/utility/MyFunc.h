//
// MyFunc.h
//
#ifndef MyFunc_H
#define MyFunc_H

#include <vector>
#include "MyMap.h"
#include <string>
#include <mpi.h>
#include <random>

using std::vector;
using std::MyMap;
using std::string;

namespace utility {

class MyFunc{
public:
	MyFunc();
	virtual ~MyFunc();

	static double mean(const vector<double> &v);
	static double sd(const vector<double> &v);
	static double var(const vector<double> &v);
	static double abs(const vector<double> &v);
	static double innerProduct(const vector<double> &v, const vector<double> &u);
	static double pearsonCorrelation(const vector<double> &v, const vector<double> &u);
	static double euclideanDist(const vector<double> &v, const vector<double> &u);

	static double pearsonCorrelationForNormalizedList(const vector<double> &v, const vector<double> &u);

	static double max(const vector<double> &v);
	static double min(const vector<double> &v);
	static double median(const vector<double> &v);
	static double mode(const vector<double> &v);

	static double percentile(const vector<double> &v, double d);
	static double upperQuartile(const vector<double> &v);
	static double lowerQuartile(const vector<double> &v);
	static double sum(const vector<double> &v);
	static vector<double> normalize( const vector<double> &v);
	static double tStatistic( const vector<double> &v, const vector<double> &u );

	static double ANOVAStatistics(const vector<vector<double> > &X);

	static string join(const string &d, const vector<string> &s);
	static void uniq(const vector<string> &O, vector<string>& result);
	static void isect(const vector<string> &A, const vector<string> &B, vector<string>& result);
	static void sample(const vector<string> &S, int n, vector<string>& result);

	static double calculatePvalueForSetOverlap(int t, int a, int b, int u);
	static void sortKeysByAscendingOrderOfValues(const MyMap<string, double>& org_map, vector<string>& result);
	static void sortKeysByDescendingOrderOfValues(const MyMap<string, double>& org_map, vector<string>& result);
	static vector<vector<string> > getAllCombination(const vector<string>& v, int n);

	static void split(const string &line, char delim, vector<string>& result);
	static void readGeneSetFromGmtFile(const string &gmtfile, MyMap<string, vector<string> >& GeneSetMap);

	// add function
	static bool contains(const vector<string>& vec_large, const vector<string>& vec_small);
	static double parseDouble(const string& str);
	static string toString(double val);
	static string toString(int val);
	static void erase(vector<string>& vstr, const string& s);

	// add function 2
	static void minmaxRow(vector<int>& minRow, vector<int>& maxRow, vector<int>& widthRow, int& size_M, int& numProcs);// use Dist, parallel M

	static void vstrPack(const vector<string>& vstr, vector<int>& vint, string& str);
	static void vstrUnpack(vector<string>& vstr, const vector<int>& vint, const string& str);
	static void vstrBcast(vector<string>& vstr, int irank, MPI_Comm comm);// Bcast for vector<string> 
	static void vdblBcast(vector<double>& vdbl, int irank, MPI_Comm comm);// Bcast for vector<double>
	static void vintBcast(vector<int>& vint, int irank, MPI_Comm comm);   // Bcast for vector<int>
	static void strBcast(string& str, int irank, MPI_Comm comm);// Bcast for string

	static void strSend(const string& str, int irank, MPI_Comm comm);
	static void strRecv(string& str, int irank, MPI_Comm comm);
	static void vstrSend(const vector<string>& vstr, int irank, MPI_Comm comm);
	static void vstrRecv(vector<string>& vstr, int irank, MPI_Comm comm);
	static void vdblSend(const vector<double>& vdbl, int irank, MPI_Comm comm);
	static void vdblRecv(vector<double>& vdbl, int irank, MPI_Comm comm);
	static void vintSend(const vector<int>& vint, int irank, MPI_Comm comm);
	static void vintRecv(vector<int>& vint, int irank, MPI_Comm comm);
	static void mstrdblSend(const MyMap<string,double>& mstrdbl, int irank, MPI_Comm comm);
	static void mstrdblRecv(MyMap<string,double>& mstrdbl, int irank, MPI_Comm comm);
};

} // namespace utility

#endif // MyFunc_H
