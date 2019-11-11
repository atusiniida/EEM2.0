//
// MyFunc.cpp
//
#include "MyFunc.h"

#include <set>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <random>

#include "HypergeometricDistributionImpl.h"

using namespace std;
using namespace utility;

MyFunc::MyFunc()
{
}
MyFunc::~MyFunc()
{
}

double MyFunc::mean(const vector<double> &v){
	double m = 0;
	for(int i = 0, n = v.size(); i < n; i++){
		m += v[i];
	}
	return m/v.size();
}
double MyFunc::sd(const vector<double> &v){
	return sqrt(var(v));
}
double MyFunc::var(const vector<double> &v){
	double var = 0;
	double m = mean(v);
	for(int i = 0, n = v.size(); i < n; i++){
		double tmp = v[i] - m;
		var += tmp * tmp;
	}
	var /= (double)v.size();
	return var;
}
double MyFunc::abs(const vector<double> &v){
	double tmp = 0;
	for(int i = 0, n = v.size(); i < n; i++){
		tmp += v[i] * v[i];
	}
	return sqrt(tmp);
}
double MyFunc::innerProduct(const vector<double> &v, const vector<double> &u){
	assert( v.size() == u.size() );
	double tmp = 0;
	for(int i = 0, n = v.size(); i < n; i++){
		tmp += v[i] * u[i];
	}
	return tmp;
}
double MyFunc::pearsonCorrelation(const vector<double> &v, const vector<double> &u){
	assert( v.size() == u.size() );
	double c = 0;
	double vm = mean(v);
	double um = mean(u);
	for(int i = 0, n = v.size(); i < n; i++){
		c +=  (v[i]-vm) * (u[i]-um);
	}
	c /= ((double)v.size())*sd(v)*sd(u);
	return c;
}
double MyFunc::euclideanDist(const vector<double> &v, const vector<double> &u){
	assert( v.size() == u.size() );
	double d = 0;
	for(int i = 0, n = v.size(); i < n; i++){
		double tmp = v[i]-u[i];
		d += tmp * tmp;
	}
	return sqrt(d);
}

double MyFunc::pearsonCorrelationForNormalizedList(const vector<double> &v, const vector<double> &u){
	assert( v.size() == u.size() );
	double c = 0;
	for(int i = 0, n = v.size(); i < n; i++){
		c +=  v[i]*u[i];
	}
	c /= (double)v.size();
	return c;
}

double MyFunc::max(const vector<double> &v){
	double  m = v[0];
	for(int i = 1, n = v.size(); i < n; i++){
		if(v[i] > m){
			m = v[i];
		}
	}
	return m;
}
double MyFunc::min(const vector<double> &v){
	double  m = v[0];
	for(int i = 1, n = v.size(); i < n; i++){
		if(v[i] <  m){
			m = v[i];
		}
	}
	return m;
}
double MyFunc::median(const vector<double> &v){
	vector<double> u(v);
	sort(u.begin(), u.end());
	int n = u.size();
	if(n%2 == 1){
		return u[(n-1)/2];
	}else{
		return (u[n/2] + u[(n/2)-1])/2;
	}
}
double MyFunc::mode(const vector<double> &v){
	vector<double> u(v);
	sort(u.begin(), u.end());
	double M = u[u.size()-1]; // max(v);
	double m = u[0];          // min(v);
	assert( m < M );
	double binNum = 1000;
	double delta = (M-m)/binNum;

	int count = 0;
	int countMax = 0;
	double  d = m+delta;
	double  dMax = d;
	for(int i = 0, n = u.size(); i < n; i++){
		double e = u[i];
		if(e <= d){
			count++;
		}else{
			if(count > countMax){
				dMax = d;
				countMax = count;
			}
			while(e > d){
				d += delta;
			}
			count = 1;
		}
	}
	dMax = dMax-(0.5*delta);
	return dMax;
}


double MyFunc::percentile(const vector<double> &v, double d){
	vector<double> u(v);
	sort(u.begin(), u.end());
	int n = u.size();
	if(d > 1.0){
		d /= 100;
	}
	int i = (int)floor(n*d)-1;
	if(i < 0){
		i = 0;
	}
	return u[i];
}
double MyFunc::upperQuartile(const vector<double> &v){
	return percentile(v,0.75);
}
double MyFunc::lowerQuartile(const vector<double> &v){
	return percentile(v,0.25);
}
double MyFunc::sum(const vector<double> &v){
	double s = 0;
	for(int i = 0, n = v.size(); i < n; i++){
		s += v[i];
	}
	return s;
}
vector<double> MyFunc::normalize( const vector<double> &v){
	vector<double> tmp(v);
	double m = mean(v);
	double s = sd(v);
	for(int i = 0, n = v.size(); i < n; i++){
		tmp[i] = (tmp[i] - m)/s;
	}
	return tmp;
}
double MyFunc::tStatistic( const vector<double> &v, const vector<double> &u ){
	double V  = sqrt(var(v)/v.size() + var(u)/u.size());
	return (mean(v) - mean(u))/V;
}


double MyFunc::ANOVAStatistics(const vector<vector<double> > &X){
	double SStotal = 0, SSbetween = 0;
	vector<double> all;
	for(int i = 0, n = X.size(); i < n; i++) {
		all.insert(all.end(), X[i].begin(), X[i].end());
	}
	double meanAll = mean(all);
	for(int i = 0, n = X.size(); i < n; i++) {
		double tmp = mean(X[i])-meanAll;
		SSbetween += X[i].size()* tmp * tmp;
		for(int j = 0, m = X[i].size(); j < m; j++){
			tmp = X[i][j]-meanAll;
			SStotal += tmp * tmp;
		}
	}
	double SSwithin = SStotal -SSbetween;
	double MSbetween = SSbetween/((double)X.size()-1);
	double MSwithin = SSwithin/((double)all.size() - X.size());
	return  MSbetween/MSwithin;
}

string MyFunc::join(const string &d, const vector<string> &s){
	string S(s[0]);
	for(int i = 1, n = s.size(); i < n; i++){
		S += d;
		S += s[i];
	}
	return S;
}

void MyFunc::uniq(const vector<string>& O, vector<string>& result){
	result.assign(O.begin(), O.end());
	sort(result.begin(), result.end());
	result.erase(unique(result.begin(), result.end()), result.end());
}

void MyFunc::isect(const vector<string> &A, const vector<string> &B, vector<string>& result){
	result.clear();
	vector<string> Asort(A), Bsort(B);
	sort(Asort.begin(), Asort.end());
	sort(Bsort.begin(), Bsort.end());
	set_intersection(Asort.begin(), Asort.end(), Bsort.begin(), Bsort.end(),
			 back_inserter(result));
}

void MyFunc::sample(const vector<string> &S, int n, vector<string>& result){
	result.assign(S.begin(), S.end());
#ifdef HAVE_SHUFFLE
	static std::random_device seed_gen;
	static std::mt19937 engine(seed_gen());
	shuffle(result.begin(), result.end(), engine);
#else
	random_shuffle(result.begin(), result.end());
#endif
	result.erase(result.begin()+n, result.end());
}

/*  Set <?> A , B, Background;
 *   calculatePvalueForSetOverlap( Background.size(), A.size(), B.size(), (isect(A,B)).size);
 */
double MyFunc::calculatePvalueForSetOverlap(int t, int a, int b, int u){
	if(!(t >  0 && a > 0 && b > 0 && u > 0  && t >= a && t >= b && a >= u && b >= u)){
		//throw new  ArithmeticException("calculatePvalueForSetOverlap: input error!");
		cerr << "calculatePvalueForSetOverlap: input error!\n";
		exit(1);
	}
	int s;
	int l;
	if(a>b){
		l = a;
		s = b;
	}else{
		l = b;
		s = a;
	}
	HypergeometricDistributionImpl H(t, l, s);
	return  H.upperCumulativeProbability(u);

}

void MyFunc::sortKeysByAscendingOrderOfValues(const MyMap<string,double> &org_map, vector<string>& result){
	typedef MyMap<string, double> org_map_t;
	typedef map<double, vector<string> > new_map_t;
	new_map_t new_map;
	for (org_map_t::const_iterator it = org_map.begin(); it != org_map.end(); it++) {
		double val = it->second;
		string key(it->first);
		new_map[val].push_back(key);
	}
	result.clear();
	for (new_map_t::iterator it2 = new_map.begin(); it2 != new_map.end(); it2++) {
		vector<string> &s = it2->second;
		result.insert(result.end(), s.begin(), s.end());
	}
}

void MyFunc::sortKeysByDescendingOrderOfValues(const MyMap<string,double> &org_map, vector<string>& result){
	sortKeysByAscendingOrderOfValues(org_map, result);
	reverse(result.begin(), result.end());
}

vector<vector<string> > MyFunc::getAllCombination(const vector<string> &v, int n) {
	vector<vector<string> > V;
	if (v.size() <= n) return V;
	set<vector<string> > seen;
	int i;
	for (i = 0; i < v.size(); i++) {
		vector<string> tmp;
		tmp.push_back(v[i]);
		if (seen.find(tmp) != seen.end()) {
			V.push_back(tmp);
			seen.insert(tmp);
		}
	}

	int j, k;
	for (i = 1; i < n; i++) {
		vector<vector<string> > TMP;
		for (j = 0; j < V.size(); j++) {
			for (k = 0; k < v.size(); k++) {
				vector<string> tmp(V[j]);
				tmp.push_back(v[k]);
				sort(tmp.begin(), tmp.end());
				if (seen.find(tmp) != seen.end()) continue;
				if (find(V[j].begin(), V[j].end(), v[k]) == V[j].end()) {
					TMP.push_back(tmp);
					seen.insert(tmp);
				}
			}
		}
		V = TMP;
	}
	return V;
}

void MyFunc::split(const string &line, char delim, vector<string>& result) {
	result.clear();
	istringstream iss(line);
	string tmp;
	while (getline(iss, tmp, delim)) {
		result.push_back(tmp);
	}
}

void MyFunc::readGeneSetFromGmtFile(const string &gmtfile, MyMap<string, vector<string> >& GeneSetMap) {
	ifstream inputStream(gmtfile.c_str());
	vector<string> str;
	GeneSetMap.clear();
	string line;
	while (getline(inputStream, line)) {
		if (line[0] == '#') continue;
		split(line, '\t', str);
		if (str.size() < 3) continue;
		string id = str[0];
		//str = MyFunc::uniq(vector<string>(str.begin()+2, str.end()));
		//GeneSetMap[id] = vector<string>(str);
		MyFunc::uniq(vector<string>(str.begin()+2, str.end()), GeneSetMap[id]);
	}
	inputStream.close();
	if (GeneSetMap.empty()) {
		//throw new DataFormatException("readGeneSetFromGmtFile: file format is wrong!");
		cerr << "readGeneSetFromGmtFile: file format is wrong!" << endl;
		exit(1);
	}
}


//--
// add function
//--
bool MyFunc::contains(const vector<string>& vec_l, const vector<string>& vec_s)
{
	vector<string>::const_iterator it = vec_s.begin();

	for(; it != vec_s.end(); it++){
		if( find(vec_l.begin(), vec_l.end(), *it) == vec_l.end() ){
			return false;
		}
	}
	return true;
}
double MyFunc::parseDouble(const string& str)
{
	double val;

	stringstream ss;
	ss << str;
	ss >> val;

	return val;
}
string MyFunc::toString(double val)
{
	string str;

	stringstream ss;
	ss << setprecision(17) << val;
	ss >> str;

	return str;
}
string MyFunc::toString(int val)
{
	string str;

	stringstream ss;
	ss << val;
	ss >> str;

	return str;
}
void MyFunc::erase(vector<string>& vstr, const string& s)
{
	vector<string>::iterator it;

	it = vstr.begin();
	while (it != vstr.end()) {
		if(*it == s) {
			it = vstr.erase(it);
		} else {
			it++;
		}
	}
}

//--
// add function 2 
//--
//
// divide M[i][j] for parallel
//
void MyFunc::minmaxRow(vector<int>& minRow, vector<int>& maxRow, vector<int>& widthRow, int& size_M, int& numProcs)
{
	//--------------------------
	// use Dist, divide M[i][j] 
	//--------------------------

	////--
	//// flat divide
	////--
	////widthRow = size_M / numProcs;
	////
	////minRow = widthRow * irank;
	////if(irank != numProcs-1){
	////	maxRow = minRow + widthRow;
	////}else{
	////	maxRow = size_M;
	////}

	//--
	// parameter init
	//--
	minRow.resize(numProcs);
	maxRow.resize(numProcs);
	widthRow.resize(numProcs);

	vector<double> vDiv;
	vDiv.resize(numProcs);
	//--
	// divide number
	//--
	for(int irank=0; irank < numProcs; irank++){
		// vDiv[irank] = 1 + (irank*2);
		vDiv[irank] = 1 + irank + sqrt(double(irank*(irank+1)));
	}
	//--
	// widthRow
	//--
	int maxRank = numProcs-1;
	for(int irank=maxRank; irank >= 0; irank--){
		int difW(0);
		for(int i=maxRank; i > irank; i--) difW += widthRow[i];

		widthRow[irank] = (size_M - difW) / vDiv[irank];
	}
	//--
	// minRow,maxRow
	//--
	for(int irank=maxRank; irank >= 0; irank--){
		int difW(0);
		for(int i=maxRank; i >= irank; i--) difW += widthRow[i];

		minRow[irank] = size_M - difW;
	}
	for(int irank=0; irank < numProcs; irank++){
		maxRow[irank] = widthRow[irank] + minRow[irank];
	}
}

//
// pack vector<string> into a single string and index
//
void MyFunc::vstrPack(const vector<string>& vstr, vector<int>& vint, string& str)
{
	vint.clear();
	str.clear();
	vector<string>::const_iterator it;
	int isum = 0;
	for (it = vstr.begin(); it != vstr.end(); it++) {
		isum += it->size();
		vint.push_back(isum);
		str += *it;
	}
}
//
// unpack vector<string> from a single string and index
//
void MyFunc::vstrUnpack(vector<string>& vstr, const vector<int>& vint, const string& str)
{
	vstr.clear();
	int istart = 0;
	vector<int>::const_iterator it;
	for (it = vint.begin(); it != vint.end(); it++) {
		int iend = *it;
		string strtmp(str, istart, iend-istart);
		vstr.push_back(strtmp);
		istart = iend;
	}
}

//
// Bcast for vector<string>
//
void MyFunc::vstrBcast(vector<string>& vstr, int irank, MPI_Comm comm)
{
	int myRank;
	MPI_Comm_rank(comm, &myRank);

	vector<int> vint;
	string str;
	if(irank == myRank) vstrPack(vstr, vint, str);

	vintBcast(vint, irank, comm);
	strBcast(str, irank, comm);

	if(irank != myRank) vstrUnpack(vstr, vint, str);
}
//
// Bcast for vector<double>
//
void MyFunc::vdblBcast(vector<double>& vdbl, int irank, MPI_Comm comm)
{
	int myRank;
	MPI_Comm_rank(comm, &myRank);

	int vecsize;
	if(irank == myRank){
		vecsize= vdbl.size();
	}
	MPI_Bcast(&vecsize, 1, MPI_INT, irank, comm);//vector size

	if(irank != myRank) vdbl.resize(vecsize);

	MPI_Bcast(&vdbl[0], vecsize, MPI_DOUBLE, irank, comm);
}
//
// Bcast for vector<int>
//
void MyFunc::vintBcast(vector<int>& vint, int irank, MPI_Comm comm)
{
	int myRank;
	MPI_Comm_rank(comm, &myRank);

	int vecsize;
	if(irank == myRank){
		vecsize= vint.size();
	}
	MPI_Bcast(&vecsize, 1, MPI_INT, irank, comm);//vector size

	if(irank != myRank) vint.resize(vecsize);

	MPI_Bcast(&vint[0], vecsize, MPI_INT, irank, comm);
}
//
// Bcast for string
//
void MyFunc::strBcast(string& str, int irank, MPI_Comm comm)
{
	int myRank;
	MPI_Comm_rank(comm, &myRank);

	int strsize;
	if(irank == myRank) strsize=str.size();

	MPI_Bcast(&strsize, 1, MPI_INT, irank, comm);//string size

	char* tmpch = new char[strsize + 1];
	if(irank == myRank){
		strcpy(tmpch, str.c_str() );
	}
	MPI_Bcast(tmpch, strsize+1, MPI_CHAR, irank, comm);

	if(irank != myRank) str=tmpch;

	delete []tmpch;
}
//
// Send & Recv for string
//
void MyFunc::strSend(const string& str, int irank, MPI_Comm comm)
{
	const char *p = str.c_str();
	int len = str.size();
	MPI_Send(const_cast<char*>(p), len, MPI_CHAR, irank, 1002, comm);
}
void MyFunc::strRecv(string& str, int irank, MPI_Comm comm)
{
	int len;
	MPI_Status status;
	MPI_Probe(irank, 1002, comm, &status);
	MPI_Get_count(&status, MPI_CHAR, &len);
	char *p = new char[len + 1];
	MPI_Recv(p, len, MPI_CHAR, irank, 1002, comm, &status);
	p[len] = '\0';
	str = p;
	delete[] p;
}
//
// Send & Recv for vector<string>
//
void MyFunc::vstrSend(const vector<string>& vstr, int irank, MPI_Comm comm)
{
	vector<int> vint;
	string str;
	vstrPack(vstr, vint, str);
	vintSend(vint, irank, comm);
	strSend(str, irank, comm);
}
void MyFunc::vstrRecv(vector<string>& vstr, int irank, MPI_Comm comm)
{
	vector<int> vint;
	string str;
	vintRecv(vint, irank, comm);
	strRecv(str, irank, comm);
	vstrUnpack(vstr, vint, str);
}
//
// Send & Recv for vector<double>
//
void MyFunc::vdblSend(const vector<double>& vdbl, int irank, MPI_Comm comm)
{
	int len = vdbl.size();
	MPI_Send(const_cast<double*>(&vdbl[0]), len, MPI_DOUBLE, irank, 1004, comm);
}
void MyFunc::vdblRecv(vector<double>& vdbl, int irank, MPI_Comm comm)
{
	int len;
	MPI_Status status;
	MPI_Probe(irank, 1004, comm, &status);
	MPI_Get_count(&status, MPI_DOUBLE, &len);
	vdbl.resize(len);
	MPI_Recv(&vdbl[0], len, MPI_DOUBLE, irank, 1004, comm, &status);
}
//
// Send & Recv for vector<int>
//
void MyFunc::vintSend(const vector<int>& vint, int irank, MPI_Comm comm)
{
	int len = vint.size();
	MPI_Send(const_cast<int*>(&vint[0]), len, MPI_INT, irank, 1006, comm);
}
void MyFunc::vintRecv(vector<int>& vint, int irank, MPI_Comm comm)
{
	int len;
	MPI_Status status;
	MPI_Probe(irank, 1006, comm, &status);
	MPI_Get_count(&status, MPI_INT, &len);
	vint.resize(len);
	MPI_Recv(&vint[0], len, MPI_INT, irank, 1006, comm, &status);
}
//
// Send & Recv for map<string,double>
//
void MyFunc::mstrdblSend(const MyMap<string,double>& mstrdbl, int irank, MPI_Comm comm)
{
	MyMap<string,double>::const_iterator it;
	vector<string> vstr;
	vector<double> vdbl;
	for (it = mstrdbl.begin(); it != mstrdbl.end(); it++) {
		vstr.push_back(it->first);
		vdbl.push_back(it->second);
	}
	vstrSend(vstr, irank, comm);
	vdblSend(vdbl, irank, comm);
}
void MyFunc::mstrdblRecv(MyMap<string,double>& mstrdbl, int irank, MPI_Comm comm)
{
	vector<string> vstr;
	vector<double> vdbl;
	vstrRecv(vstr, irank, comm);
	vdblRecv(vdbl, irank, comm);
	mstrdbl.clear();
	int len = vstr.size();
	for (int i = 0; i < len; i++) {
		mstrdbl[vstr[i]] = vdbl[i];
	}
}
