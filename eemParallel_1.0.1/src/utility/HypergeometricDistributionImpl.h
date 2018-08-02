//
// HyperGeometricDistributionImpl.cpp
//
#ifndef HyperGeometricDistributionImpl_H
#define HyperGeometricDistributionImpl_H

namespace utility {

class HypergeometricDistributionImpl {
private:
	int numberOfSuccesses;
	int populationSize;
	int sampleSize;

public:
	HypergeometricDistributionImpl();
	virtual ~HypergeometricDistributionImpl();
	HypergeometricDistributionImpl(int populationSize,
				       int numberOfSuccesses, int sampleSize);
	void setNumberOfSuccesses(int num);
	void setPopulationSize(int size);
	void setSampleSize(int size);
	int getNumberOfSuccesses();
	int getPopulationSize();
	int getSampleSize();
	double probability(int x);
	double cumulativeProbability(int x);
	double upperCumulativeProbability(int x);

private:
	int getLowerDomain(int n, int m, int k);
	int getUpperDomain(int m, int k);
	double probability(int n, int m, int k, int x);
	double innerCumulativeProbability(
		int x0, int x1, int dx, int n, int m, int k);
};

} // namespace utility

#endif // HyperGeometricDistributionImpl_H
