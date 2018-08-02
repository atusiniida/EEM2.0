//
// MyMat.cpp
//
#include "MyMat.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "MyFunc.h"
#include "MyException.h"

using namespace std;
using namespace utility;

MyMat::MyMat()
{
}
MyMat::~MyMat()
{
}

MyMat::MyMat(const string &infile) throw ()
{
	ifstream ifs(infile.c_str());
	if (!ifs) {
		cerr << "file opening failed " << infile << endl;
		// TODO: throw() or exit()
		exit(EXIT_FAILURE);
	}

	vector<string> lines;
	string line;
	nrow = -1;
	while( getline(ifs,line) ){// != null){
		if(line[0] == '#'){
			continue;
		}
		lines.push_back(line);
		nrow++;
	}

	int i;

	int l = 0;
	line = lines[0];

	char delim = '\t';
	vector<string> str;
	MyFunc::split(line, delim, str);

	ncol = str.size()-1;
	for(i = 1; i < str.size(); i++){
		if(find(colname.begin(), colname.end(), str[i]) != colname.end() ){
			cerr << "MyMat: file format is wrong (colnames are not unique)! \n" << endl;
		}
		colname.push_back(str[i]);
	}

	M.resize(nrow);
	for(size_t ii=0; ii < nrow; ii++) M[ii].resize(ncol);

	for(l=1; l <= nrow; l++){
		line = lines[l];

		MyFunc::split(line, delim, str);

		if(str.size() != ncol+1){
			cerr << "MyMat: file format is wrong (column size is inconsistent) !\n" << endl;
		}
		if( find(rowname.begin(), rowname.end(), str[0]) != rowname.end() ){
			cerr << "MyMat: file format is wrong (rownames are not unique)! \n" << endl;
		}
		rowname.push_back(str[0]);

		for(i=1; i < str.size(); i++){
			try{
				M[l-1][i-1] = MyFunc::parseDouble(str[i]);

			}catch(MyException e){// NumberFormatException e
				cerr << "MyMat: file format is wrong (an element is not Double)! \n" << endl;
			}
		 }
	}
	for(i=0;i<ncol;i++){
		colname2index[colname[i]] = i;
	}
	for(i=0;i<nrow;i++){
		rowname2index[rowname[i]] = i;
	}

}

MyMat::MyMat(const MyMat& m)
{
	nrow = m.nrow;
	ncol = m.ncol;
	colname2index = m.colname2index;
	rowname2index = m.rowname2index;
	colname = m.colname;
	rowname = m.rowname;

	int i,j;
	M.resize(nrow);//M = new double[nrow][ncol];
	for(i=0; i < nrow; i++){
		M[i].resize(ncol);
		for(j=0;j<ncol;j++){
			M[i][j] = m.M[i][j];
		}
	}
}


double MyMat::get(int i, int j) const
{
	if(i >= nrow || j >= ncol){
		cerr << "MyMat::get(int i, int j) " << endl;
		throw MyException();
	}
	return M[i][j];
}

const vector<string>& MyMat::getRowNames() const
{
	return rowname;
}

const vector<string>& MyMat::getColNames() const
{
	return colname;
}

int MyMat::rowSize() const
{
	return nrow;
}

int MyMat::colSize() const
{
	return ncol;
}


const vector<double>& MyMat::getRow(int i) const
{
	return M[i];
}

const vector<double>& MyMat::getRow(const string &i) const
{
	static vector<double> tmp;
	MyMap<string, int>::const_iterator it;
	it = rowname2index.find(i);
	if (it == rowname2index.end()) return tmp;
	return M[ it->second ];
}

void MyMat::normalizeRows()
{
	int i,j;
	for(i=0; i<nrow; i++){
		vector<double> tmp = getRow(i);
		double tmp2 = MyFunc::mean(tmp);
		double tmp3 = MyFunc::sd(tmp);

		for(j=0; j<ncol; j++){
			M[i][j] = (M[i][j] - tmp2)/tmp3;
		}
	}
}

vector<double> MyMat::getRowMeans(const vector<string> &rows)
{
	vector<double> v;
	int i,j;
	for(j=0; j<ncol; j++){
		double tmp = 0.0;
		for(i=0; i < rows.size(); i++){
		  tmp += get(rowname2index[rows[i]], j);
		}
		if(rows.size() > 0) {
			v.push_back(tmp/(double)rows.size()) ;
		} else {
			v.push_back(0.0);
		}
	}
	return v;
}

//--
// operator 2013.11.16
//--
MyMat& MyMat::operator =(const MyMat &mat)
{
	M = mat.M;
	colname2index = mat.colname2index;
	rowname2index = mat.rowname2index;
	ncol = mat.ncol;
	nrow = mat.nrow;
	colname = mat.colname;
	rowname = mat.rowname;

	return *this;
}
