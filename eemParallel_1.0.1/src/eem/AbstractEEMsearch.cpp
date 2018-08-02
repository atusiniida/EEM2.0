//
// AbstractEEMsearch.cpp
//
#include "AbstractEEMsearch.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include "../utility/MyFunc.h"
#include "../utility/MyException.h"

#include "EEM.h"
#include "ExpressionModule.h"
#include "ExpressionModuleSet.h"

using namespace std;
using namespace utility;
using namespace eem;

AbstractEEMsearch::AbstractEEMsearch()
{
	maxSeedGeneSize = 2000;
	minSeedGeneSize = 10;
	itrForPvalue2Calculation = 500;
	Pvalue1Cutoff = 1.0;

	//stopWatch = new StopWatch();
	//seedGeneSets = new HashMap<String, List<String>>();
	//eems = new HashMap<String, EEM>();
	//seeds = new ArrayList<String>();
	//candidates  = new ArrayList<String>();

	errLog = "";
	timeLog = "";

	nullDistrbutionForRecycle.clear();
}
AbstractEEMsearch::~AbstractEEMsearch()
{
	//TODO : change safe delete
	MyMap<string, EEM*>::iterator it;
	for(it=eems.begin(); it != eems.end(); it++){
		EEM *eem = it->second;
		delete eem;
	}
}

AbstractEEMsearch::AbstractEEMsearch(const AbstractEEMsearch& A)
{
	originalExp = A.originalExp;
	allGenes = A.allGenes;
	maxSeedGeneSize = A.maxSeedGeneSize;
	minSeedGeneSize = A.minSeedGeneSize;
	itrForPvalue2Calculation = A.itrForPvalue2Calculation;
	Pvalue1Cutoff = A.Pvalue1Cutoff;

	seedGeneSets = A.seedGeneSets;

	errLog = "";
	timeLog = "";
}



void AbstractEEMsearch::setGeneSets(const MyMap<string, vector<string> >& geneSets)
{
	MyMap<string, vector<string> >::const_iterator it;
	for(it=geneSets.begin(); it != geneSets.end(); it++){

		vector<string> tmp;
		MyFunc::isect(it->second, allGenes, tmp);

		if(tmp.size() < minSeedGeneSize || tmp.size() > maxSeedGeneSize){
			stringstream ss;
			ss << it->first << ": seed geneset size (" << tmp.size() << ") is out of range!\n";
			cerr << ss.str();

			ss.str("");
			ss.clear(stringstream::goodbit);
			ss << it->first << ": seed geneset size (" << tmp.size() << ") is out of range!\n";
			errLog += ss.str();

			continue;
		}
		seedGeneSets[it->first] = tmp;
	}
}

void AbstractEEMsearch::setPvalue1Cutoff(double d)
{
	Pvalue1Cutoff = d;
}

void AbstractEEMsearch::suppressPvalue1Cutoff()
{
	Pvalue1Cutoff = 0.0;
}

void AbstractEEMsearch::setItrForPvalue2Calculation(int i) 
{
	itrForPvalue2Calculation= i;
}

void AbstractEEMsearch::recycleNullDistribution()
{
	nullDistrbutionForRecycle.clear();// = new HashMap<Integer, List<Integer>>();
}

// parallelized
void AbstractEEMsearch::findModuleGenes(MPI_Comm comm)
{
	cerr << "Findng module gene..." << endl;

	double stctime,edctime;//calculate time
	double stdtime,eddtime;//distribute time

	int i;
	int n = eems.size();
	vector<string> new_seeds;// = new ArrayList<String>();

	vector<int> seeds_flg(seeds.size(), 0);

	int numProcs, myRank;
	MPI_Comm_size(comm, &numProcs);
	MPI_Comm_rank(comm, &myRank);

	MPI_Barrier(comm);
	stctime = MPI_Wtime();

	for(i=myRank; i < seeds.size(); i+=numProcs){// single => (i=0; i < seeds.size(); i++)
		string s = seeds[i];
		try{
			eems[s]->findModuleGenes();// EEM*

			stringstream ss;
			ss << s << "(" << i+1 << "/" << n << "): succeed! " << eems[s]->getModuleGenes().size() << "/" << eems[s]->getSeedGenes().size() << "\n";
			cerr << ss.str();

			seeds_flg[i] = myRank + 1;
		}
		catch (MyException &err) {
			stringstream ss;
			ss << s << "(" << i+1 << "/" << n << "): failed!";

			cerr << err.what_message(ss.str()) << endl;

			errLog += s +  ": unable to find module genes!\n"; 
			continue;
		}
	}
	edctime = MPI_Wtime();
	MPI_Barrier(comm);
	stdtime = MPI_Wtime();
	{
		int *ibuf = new int[seeds.size()];
		MPI_Allreduce(&seeds_flg[0], ibuf, seeds.size(), MPI_INT, MPI_SUM, comm);
		for (i = 0; i < seeds.size(); i++) seeds_flg[i] = ibuf[i];
		delete[] ibuf;
	}

	for (i = 0; i < seeds.size(); i++) {
		int irank = seeds_flg[i] - 1;
		if (irank < 0) continue;
		string s = seeds[i];
		new_seeds.push_back(s);

		vector<string> vModule;//eem->moduleGenes
		vector<double> vCenter;//eem->center
		string         sCenGen;//eem->centerGene
		EEM *eemp = eems[s];

		if(irank == myRank){
			vModule = eemp->getModuleGenes();
			vCenter = eemp->getCenter();
			sCenGen = eemp->getCenterGene();
		}
		//MPI communication
		MyFunc::vstrBcast(vModule, irank, comm);
		MyFunc::vdblBcast(vCenter, irank, comm);
		MyFunc::strBcast(sCenGen, irank, comm);

		if(irank != myRank){
			eemp->setModuleGenes(vModule);
			eemp->setCenter(vCenter);
			eemp->setCenterGene(sCenGen);
		}
	}
	eddtime = MPI_Wtime();

	seeds = new_seeds;

	if(seeds.empty()){
		//throw new MyException("findModuleGenes: Any module has passed!");
		throw MyException();// TODO
	}

	stringstream ss;
	ss << "AbstractEEMsearch::findModuleGenes: proc=" << myRank << "/" << numProcs
	   << ", total time:" << eddtime - stctime
	   << ", calculate time:" << edctime - stctime
	   << ", wait time:" << stdtime - edctime
	   << ", distribute time:" << eddtime - stdtime << "\n";
	cerr << ss.str();
}

// parallelized
void AbstractEEMsearch::calculatePvalue1(MPI_Comm comm)
{
	double stctime,edctime;//calculate time
	double stdtime,eddtime;//distribute time

	if(Pvalue1Cutoff == 0.0){
		return;
	}else{
		cerr << "Calculating approximate P values to filter out non-significant genesets..." << endl;
		int i;
		int n = seeds.size();
		vector<string> new_seeds;// = new ArrayList<String>();

		vector<int> seeds_flg(seeds.size(), 0);

		int numProcs, myRank;
		MPI_Comm_size(comm, &numProcs);
		MPI_Comm_rank(comm, &myRank);

		MPI_Barrier(comm);
		stctime = MPI_Wtime();

		for(i=myRank; i < seeds.size(); i+=numProcs){ // single process -> (i=0; i < seeds.size(); i++)
			string s = seeds[i];
			try{
				eems[s]->calculatePvalue1();

				stringstream ss;
				ss << s << "(" << i+1 << "/" << n << "): " << "succeed! " << eems[s]->getPvalue1() << "\n";
				cerr << ss.str();

				seeds_flg[i] = myRank + 1;
			}
			catch (MyException &err) {
				stringstream ss;
				ss << s << "(" << i << "/" << n << "): failed!";

				cerr << err.what_message(ss.str()) << endl;

				errLog +=  s +  ": unable to calculate a P value  based on hypergeometric distribution!\n";
				continue;
			}
		}
		edctime = MPI_Wtime();

		MPI_Barrier(comm);
		stdtime = MPI_Wtime();
		{
			int *ibuf = new int[seeds.size()];
			MPI_Allreduce(&seeds_flg[0], ibuf, seeds.size(), MPI_INT, MPI_SUM, comm);
			for (i = 0; i < seeds.size(); i++) seeds_flg[i] = ibuf[i];
			delete[] ibuf;
		}

		for (i = 0; i < seeds.size(); i++) {
			int irank = seeds_flg[i] - 1;
			if (irank < 0) continue;
			string s = seeds[i];
			new_seeds.push_back(s);

			double valP[2];
			EEM *eemp = eems[s];

			if(irank == myRank){
				valP[0] = eemp->getPvalue();
				valP[1] = eemp->getPvalue1();
			}
			MPI_Bcast(valP, 2, MPI_DOUBLE, irank, comm);

			if(irank != myRank){
				eemp->setPvalue(valP[0]);
				eemp->setPvalue1(valP[1]);
			}
		}
		eddtime = MPI_Wtime();

		seeds = new_seeds;

		if(seeds.empty()){
			//throw new MyException(" calculatePvalue1: Any module has passed!");
			throw MyException();// TODO
		}
		stringstream ss;
		ss << "AbstractEEMsearch::calculatePvalue1: proc=" << myRank << "/" << numProcs
		   << ", total time:" << eddtime - stctime
		   << ", calculate time:" << edctime - stctime
		   << ", wait time:" << stdtime - edctime
		   << ", distribute time:" << eddtime - stdtime << "\n";
		cerr << ss.str();

	}//end if(Pvalue1Cutoff == 0.0)
}

void AbstractEEMsearch::findCandidates()
{
	if(Pvalue1Cutoff == 0.0){
		candidates = seeds;// candidates = new ArrayList<String>(seeds);
	}else{
		for(size_t i=0; i < seeds.size(); i++){
			string s = seeds[i];
			if(eems[s]->getPvalue1() >= Pvalue1Cutoff){
				candidates.push_back(s);
			}
		}
		if(candidates.empty()){
			//throw MyException("findCandidates: Any module has passed!");
			cerr << "findCandidates: Any module has passed!" << endl;
		}
	}
}

void AbstractEEMsearch::calculatePvalue2_inner(int myRank, int i, vector<int>& candidates_flg)
{
	int n = candidates.size();
	string s;
	try{
		s = candidates[i];
		eems[s]->calculatePvalue2();

		stringstream ss;
		ss << s << "(" << i+1 << "/"  << n << "): " << "succeed! " << eems[s]->getPvalue2() << "\n";
		cerr << ss.str();

		candidates_flg[i] = myRank + 1;
	}
	catch (MyException &err) {
		// err.printStackTrace();
		stringstream ss;
		ss << s << "(" << i+1 << "/" << n << "): failed!";
		cerr << err.what_message( ss.str() ) << endl;

		errLog +=  s +  ": unable to calculate a P value based on arandomization test!\n"; 
		//continue;
	}
}

// parallelized
void AbstractEEMsearch::calculatePvalue2(MPI_Comm comm)
{
	cerr << "Calculating accurate P values..." << endl;

	int i;
	vector<string> new_candidates;// = new ArrayList<String>();

	vector<int> candidates_flg(candidates.size(), 0);

	int numProcs, myRank;
	MPI_Comm_size(comm, &numProcs);
	MPI_Comm_rank(comm, &myRank);

	MPI_Barrier(comm);
	double stctime = MPI_Wtime();

#ifndef SWITCHING_NUM_PROCS
#  define SWITCHING_NUM_PROCS 4
#endif
	if(numProcs < SWITCHING_NUM_PROCS || candidates.size() < numProcs) { /* STATIC TASK DISTRIBUTION */
		for(i=myRank; i < candidates.size(); i+=numProcs){ // single process -> (i=0; i < n; i++){
			calculatePvalue2_inner(myRank, i, candidates_flg);
		}
	} else { /* DYNAMIC TASK DISTRIBUTION */
		if(myRank == 0) { /* Parent process */
			int next = numProcs - 1;
			int num_finished = 0;
			int finished_self = 1; // Rank 0 is NON_WORKING_HOST; just distribute jobs
			// if(numProcs < SWITCHING_NUM_PROCS) {
			// 	finished_self = 0; // Rank 0 is WORKING_HOST
			// }
			while(num_finished < numProcs - 1 || !finished_self) {
				MPI_Status status;
				int flag = 0;
				if (num_finished < numProcs - 1) {
					/* Check if any process has finished the assigned task */
					MPI_Iprobe(MPI_ANY_SOURCE, 0, comm, &flag, &status);
				}
				if(flag) { /* Distribute new task or send finish-sign */
					int irank = status.MPI_SOURCE;
					int todo;
					/* Receive which job has finished */
					MPI_Recv(&todo, 1, MPI_INT, irank, 0, comm, &status);
					if(next < candidates.size()) { /* if any task left */
						todo = next; /* next task */
						next++;
					} else { /* no more task */
						todo = -1; /* finish-sign */
						num_finished++;
					}
					MPI_Send(&todo, 1, MPI_INT, irank, 0, comm);
				} else if (!finished_self) { /* Assign task to myself */
					if (next < candidates.size()) { /* if any task left */
						stringstream ss;
						ss << next+1 << " is assigned to local-rank 0\n";
						cerr << ss.str();
						calculatePvalue2_inner(myRank, next, candidates_flg);
						next++;
					} else { /* no more task */
						finished_self = 1;
						//cerr << "local-rank 0: finished\n";
					}
				}
			}
		} else { /* Child processes */
			int todo = myRank - 1;
			while(todo >= 0) {
				stringstream ss;
				MPI_Status status;
				ss << todo+1 << " is assigned to local-rank " << myRank << "\n";
				cerr << ss.str();
				calculatePvalue2_inner(myRank, todo, candidates_flg);
				/* tell parent that the assigned task is done */
				MPI_Send(&todo, 1, MPI_INT, 0, 0, comm);
				/* receive new task */
				MPI_Recv(&todo, 1, MPI_INT, 0, 0, comm, &status);
			}
			//cerr << "local-rank " << myRank << ": finished\n";
		}
	}
	double edctime = MPI_Wtime();

	MPI_Barrier(comm);
	double stdtime = MPI_Wtime();
	{
		int *ibuf = new int[candidates.size()];
		MPI_Allreduce(&candidates_flg[0], ibuf, candidates.size(), MPI_INT, MPI_SUM, comm);
		for (i = 0; i < candidates.size(); i++) candidates_flg[i] = ibuf[i];
		delete[] ibuf;
	}

	for (i = 0; i < candidates.size(); i++) {
		int irank = candidates_flg[i] - 1;
		if (irank < 0) continue;
		string s = candidates[i];
		new_candidates.push_back(s);

		double valP[2];
		vector<int> nDist, denDist;
		vector<string> vModule;//eem->moduleGenes
		vector<double> vCenter;//eem->center
		string         sCenGen;//eem->centerGene
		EEM *eemp = eems[s];

		if(irank == myRank){
			valP[0] = eemp->getPvalue();
			valP[1] = eemp->getPvalue2();
			nDist = eemp->getNullDist();
			denDist= eemp->get_densityDist();
			vModule= eemp->getModuleGenes();
			vCenter= eemp->getCenter();
			sCenGen= eemp->getCenterGene();
		}
		MPI_Bcast(valP,  2, MPI_DOUBLE, irank, comm);
		//
		// MyFunc:: *** Bcast 
		//
		MyFunc::vintBcast(nDist, irank, comm);
		MyFunc::vintBcast(denDist, irank, comm);
		MyFunc::vstrBcast(vModule, irank, comm);
		MyFunc::vdblBcast(vCenter, irank, comm);
		MyFunc::strBcast(sCenGen, irank, comm);

		if(irank != myRank){
			eemp->setPvalue(valP[0]);
			eemp->setPvalue2(valP[1]);
			eemp->setNullDist(nDist);
			eemp->set_densityDist(denDist);
			eemp->setModuleGenes(vModule);
			eemp->setCenter(vCenter);
			eemp->setCenterGene(sCenGen);

			eemp->setParentNullDist(nDist);
		}
	}
	double eddtime = MPI_Wtime();

	stringstream ss;
	ss << "AbstractEEMsearch::calculatePvalue2: proc=" << myRank << "/" << numProcs
	   << ", total time:" << eddtime - stctime
	   << ", calculate time:" << edctime - stctime
	   << ", wait time:" << stdtime - edctime
	   << ", distribute time:" << eddtime - stdtime << "\n";
	cerr << ss.str();

	candidates = new_candidates;

	if(candidates.empty()){
		//throw MyException("calculatePvalue2: Any module has passed!");
		cerr << "calculatePvalue2: Any module has passed!" << endl;
	}
}

void AbstractEEMsearch::printResults(const string& outfile) throw ()
{
	MyMap<string, double> p;

	MyMap<string, EEM*>::iterator it;
	for(it=eems.begin(); it != eems.end(); it++){
		p[it->first] = it->second->getPvalue();
	}
	// tmp <= p
	vector<string> tmp;
	MyFunc::sortKeysByDescendingOrderOfValues( p, tmp );

	// tmp => outfile
	ofstream ofs;
	ofs.open(outfile.c_str());
	size_t i;
	for( i=0; i < tmp.size(); i++ ){
		string id = tmp[i];
		EEM *eemp = eems[id];

		stringstream ss;
		ss << id << "\t" << eemp->getModuleGenes().size() << "/" << eemp->getSeedGenes().size()
		   << "\t" << eemp->getPvalue() << "\t" << MyFunc::join("\t", eemp->getModuleGenes()) << "\n";

		ofs << ss.str();
	}
	ofs.close();
}

void AbstractEEMsearch::printLog(const string& outfile) throw ()
{
	ofstream ofs;

	ofs.open(outfile.c_str(), ios::out|ios::app);
	ofs << getLog() << endl;
	ofs.close();
}

void AbstractEEMsearch::setEEM(MPI_Comm comm)
{
	//
	// throw new UnsupportedOperationException();
	//
	cerr << "Unsupported Operation  Exception : void AbstractEEMsearch::setEEM()" << endl;
}
//--
// generate ExpressionModuleSet, ExpressionModule
//--
ExpressionModuleSet* AbstractEEMsearch::getExpressionModuleSet()
{
	ExpressionModuleSet* tmp = new ExpressionModuleSet();// generate ExpressionModuleSet

	MyMap<string, EEM*>::iterator it;

	for(it=eems.begin(); it != eems.end(); it++){
		if(it->second->getPvalue() >= 0.0){
			it->second->cutParent();
			tmp->add(new ExpressionModule(it->first, it->second, originalExp));// generate ExpressionModule
		}
	}
	return tmp;
}

string AbstractEEMsearch::getLog()
{
	return timeLog + "\n" + this->name() + "\n" + errLog;
}

void AbstractEEMsearch::writeTimeLog()
{
	timeLog = "started at " +  stopWatch.getStartDate().toString() + "\n" +
			"finished at " + stopWatch.getStopDate().toString() + "\n" +
			"It took " + stopWatch.toString() + "\n";
}

void AbstractEEMsearch::perform(MPI_Comm comm)
{
	try{
		setEEM(comm);
		findModuleGenes(comm);
		calculatePvalue1(comm);
		findCandidates();
		calculatePvalue2(comm);

	}catch(MyException& err){
		cerr << err.what() << endl;

		stopWatch.stop();
		writeTimeLog();
	}
}

MyMap<string, double> AbstractEEMsearch::getPvalues()
{
	MyMap<string, double> P;

	MyMap<string, EEM*>::iterator it;
	for(it=eems.begin(); it != eems.end(); it++){

		if(it->second->getPvalue() >= 0.0){
			P[it->first]= it->second->getPvalue();
		}
	}

	return P;
}
