//
// MyMat.h
//
#ifndef MyMat_H
#define MyMat_H

#include <vector>
#include "MyMap.h"
#include <string>

using std::vector;
using std::MyMap;
using std::string;

namespace utility{
class MyMat{
public:
	MyMat();
	MyMat(const string &infile) throw ();// throws IOException, DataFormatException
	MyMat(const MyMat& m);
	virtual ~MyMat();

private:
	static long long serialVersionUID;// = -6849257508515306759L;

protected:
	vector<vector<double> > M;
	MyMap<string, int> colname2index;
	MyMap<string, int> rowname2index;
	int ncol;
	int nrow;
	vector<string> colname;
	vector<string> rowname;

public:
	double get(int i, int j) const;
	const vector<string>& getRowNames() const;
	const vector<string>& getColNames() const;

	int rowSize() const;
	int colSize() const;

	const vector<double>& getRow(int i) const;
	const vector<double>& getRow(const string &i) const;

	void normalizeRows();

	vector<double> getRowMeans(const vector<string> &rows);

	//--
	// operator 2013.11.16
	//--
	MyMat& operator=(const MyMat& mat);
};
} // namespace utility

#endif // MyMat_H
