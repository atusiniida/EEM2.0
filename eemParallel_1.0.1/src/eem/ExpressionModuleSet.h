//
// ExpressionModuleSet.h
//
#ifndef ExpressionModuleSet_H
#define ExpressionModuleSet_H

#include "../utility/MyMap.h"
#include <vector>
#include <string>

using std::MyMap;
using std::vector;
using std::string;

#include "ExpressionModule.h"

//
// class ExpressionModuleSet implements Serializable
//
namespace eem{
///
/// an ExpressionModuleSet class object is a container of ExpressionModule class objects
///
class ExpressionModuleSet{
public:
	ExpressionModuleSet();
	virtual ~ExpressionModuleSet();

private:

	MyMap<string, ExpressionModule*> expressionModules;
	vector<string> ids;

public:
	// add expressionModule
	// if an expressionModule with same ID exists, set the better one.
	void add(ExpressionModule* expMod);

	ExpressionModule* get(const string& id);
	ExpressionModule* get(int i);

	const vector<string>& getIds() const;

	// void writeToFile(const string& file) throw ();
};
} // namespace eem
#endif // ExpressionModuleSet_H
