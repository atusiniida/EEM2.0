//
// HyperGeometricDistributionImpl.cpp
//
#include "HypergeometricDistributionImpl.h"

#include <string>
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;
using namespace utility;

HypergeometricDistributionImpl::HypergeometricDistributionImpl() {
}

HypergeometricDistributionImpl::~HypergeometricDistributionImpl() {
}

HypergeometricDistributionImpl::HypergeometricDistributionImpl(int populationSize,
			       int numberOfSuccesses, int sampleSize) {
	if (numberOfSuccesses > populationSize) {
		cerr << "number of successes (" << numberOfSuccesses
		     << ") must be less than or equal to population size ("
		     << populationSize << ")\n";
		exit(1);
	}
	if (sampleSize > populationSize) {
		cerr << "sample size (" << sampleSize
		     << ") must be less than or equal to population size ("
		     << populationSize << ")\n";
		exit(1);
	}
	setPopulationSize(populationSize);
	setSampleSize(sampleSize);
	setNumberOfSuccesses(numberOfSuccesses);
}

void HypergeometricDistributionImpl::setNumberOfSuccesses(int num) {
	if(num < 0){
		cerr << "number of successes must be non-negative (" << num << ")\n";
		exit(1);
	}
	numberOfSuccesses = num;
}

void HypergeometricDistributionImpl::setPopulationSize(int size) {
	if(size <= 0){
		cerr << "population size must be positive (" << size << ")\n";
		exit(1);
	}
	populationSize = size;
}

void HypergeometricDistributionImpl::setSampleSize(int size) {
	if (size < 0) {
		cerr << "sample size must be positive (" << size << ")\n";
		exit(1);
	}
	sampleSize = size;
}

int HypergeometricDistributionImpl::getNumberOfSuccesses() {
	return numberOfSuccesses;
}

int HypergeometricDistributionImpl::getPopulationSize() {
	return populationSize;
}

int HypergeometricDistributionImpl::getSampleSize() {
	return sampleSize;
}

double HypergeometricDistributionImpl::probability(int x) {
	double ret;

	int n = getPopulationSize();
	int m = getNumberOfSuccesses();
	int k = getSampleSize();

	int domain[2] = {
		getLowerDomain(n, m, k),
		getUpperDomain(m, k)
	};
	if(x < domain[0] || x > domain[1]){
		ret = 0.0;
	} else {
		ret = probability(n, m, k, x);
	}

	return ret;
}

double HypergeometricDistributionImpl::cumulativeProbability(int x) {
	double ret;

	int n = getPopulationSize();
	int m = getNumberOfSuccesses();
	int k = getSampleSize();

	int domain[2] = {
		getLowerDomain(n, m, k),
		getUpperDomain(m, k)
	};
	if (x < domain[0]) {
		ret = 0.0;
	} else if(x >= domain[1]) {
		ret = 1.0;
	} else {
		ret = innerCumulativeProbability(domain[0], x, 1, n, m, k);
	}

	return ret;
}

double HypergeometricDistributionImpl::upperCumulativeProbability(int x) {
	double ret;

	int n = getPopulationSize();
	int m = getNumberOfSuccesses();
	int k = getSampleSize();

	int domain[2] = {
		getLowerDomain(n, m, k),
		getUpperDomain(m, k)
	};
	if (x < domain[0]) {
		ret = 1.0;
	} else if(x > domain[1]) {
		ret = 0.0;
	} else {
		ret = innerCumulativeProbability(domain[1], x, -1, n, m, k);
	}

	return ret;
}

static int max(int m, int n) {
	return m < n ? n : m;
}

static int min(int m, int n) {
	return m < n ? m : n;
}

int HypergeometricDistributionImpl::getLowerDomain(int n, int m, int k) {
	return max(0, m - (n - k));
}

int HypergeometricDistributionImpl::getUpperDomain(int m, int k){
	return min(k, m);
}

//static double binomialCoefficient(int n, int k) {
//}
//static double binomialCoefficientDouble(int n, int k) {
//}

static double binomialCoefficientLog(int n, int k) {
	if (n < k) {
		cerr << "invalid parameters order in binomial (" << n << ", " << k << ")\n";
		exit(1);
	}
	if (n < 0) {
		cerr << "negative parameter in binomial (" << n << ")\n";
		exit(1);
	}

	if ((n == k) || (k == 0)) return 0;
	if ((k == 1) || (k == n - 1)) return log((double)n);

	//if (n < 67) return log(binomialCoefficient(n,k));
	//if (n < 1030) return log(binomialCoefficientDouble(n, k));

	if (k > n / 2) return binomialCoefficientLog(n, n - k);

	double logSum = 0;

	// n!/(n-k)!
	for (int i = n - k + 1; i <= n; i++) {
		logSum += log((double)i);
	}

	// divide by k!
	for (int i = 2; i <= k; i++) {
		logSum -= log((double)i);
	}

	return logSum;
}

double HypergeometricDistributionImpl::probability(int n, int m, int k, int x) {
	return exp(binomialCoefficientLog(m, x) +
		   binomialCoefficientLog(n - m, k - x) -
		   binomialCoefficientLog(n, k));
}

double HypergeometricDistributionImpl::innerCumulativeProbability(
	int x0, int x1, int dx, int n, int m, int k)
{
	double ret = probability(n, m, k, x0);
	while (x0 != x1) {
		x0 += dx;
		ret += probability(n, m, k, x0);
	}
	return ret;
}
