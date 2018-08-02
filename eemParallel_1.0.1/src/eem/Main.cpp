//
// Main.cpp (eemParallel)
//
// Define Macro "HAVE_UNORDERED_MAP" for (C++11 std) unordered_map
// Define Macro "VC" for Visual C++
//

#ifdef VC
#include "wgetopt.h"
#else
#include <getopt.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <mpi.h>

#include "../utility/MyFunc.h"
#include "../utility/MyMat.h"
#include "../utility/Dist.h"

#include "CoherenceBasedEEMsearch.h"
#include "ExpressionModuleSet.h"

using namespace std;
using namespace utility;
using namespace eem;

#ifndef VC
static void usage(const string& cmdline, const struct option* options)
{
	cerr << cmdline << endl;
	cerr << "Options:" << endl;
	while (options->name) {
		cerr << "--" << options->name;
		if (options->val) {
			cerr << ", -" << char(options->val);
		}
		if (options->has_arg == required_argument) {
			cerr << " arg";
		}
		//
		cerr << endl;
		options++;
	}
}
#endif

static int extract_relrad(const char* optarg, vector<double>& relrads)
{
	int n_relrad = 0;
	for (const char *p = optarg; *p != '\0'; p++) {
		if (*p == ',') {
			n_relrad++;
		}
	}
	n_relrad++;
	stringstream ss;
	ss.str(string(optarg));
	for (int i = 0; i < n_relrad; i++) {
		double tmp;
		ss >> tmp;
		ss.ignore(256,',');
		relrads.push_back(tmp);
	}
	return n_relrad;
}

int main(int argc, char* argv[])
{
	//--
	// MPI initialize
	//--
	MPI_Init(&argc,&argv);

	MPI_Barrier(MPI_COMM_WORLD);
	double stctime= MPI_Wtime();

	int numProcs, myRank;
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	//--
	// options
	//--
	bool hasOption_R = false;
	int n_relrad = 0;
	vector<double> relrads;
	bool hasOption_p = false;
	double pcut = 0.0;
	bool hasOption_i = false;
	int itr = 0;
	bool hasOption_o = false;
	string outfile;
	bool hasOption_m = false;
	int mingeneset = 0;
	bool hasOption_M = false;
	int maxgeneset = 0;
	bool hasOption_l = false;
	string logfile;
	bool hasOption_c = false;
	int coregenesize = 0;

#ifndef VC
	static struct option options[] = {
		{"relrad",     required_argument, 0, 'R'}, // multiple relativeRadius should be given l like  "-R 0.05,0.10,0.15"
		{"pcut",       required_argument, 0, 'p'},
		{"itr",        required_argument, 0, 'i'},
		{"outfile",    required_argument, 0, 'o'},
		{"mingeneset", required_argument, 0, 'm'},
		{"maxgeneset", required_argument, 0, 'M'},
		{"log",        required_argument, 0, 'l'},
		{0,            0,                 0, 0  }
	};
#endif

	while (1) {
		int option_index = 0;

		int c = getopt_long(argc, argv, "R:p:i:o:m:M:l:",
						options, &option_index);

		if (c == -1)
			break;

		switch (c) {
		case 'R':
			hasOption_R = true;
			n_relrad = extract_relrad(optarg, relrads);
			if (n_relrad <= 0) return 1;
			break;
		case 'p':
			hasOption_p = true;
			pcut = atof(optarg);
			break;
		case 'i':
			hasOption_i = true;
			itr = atoi(optarg);
			break;
		case 'o':
			hasOption_o = true;
			outfile = string(optarg);
			break;
		case 'm':
			hasOption_m = true;
			mingeneset = atoi(optarg);
			break;
		case 'M':
			hasOption_M = true;
			maxgeneset = atoi(optarg);
			break;
		case 'l':
			hasOption_l = true;
			logfile = string(optarg);
			break;
		case 'c':
			hasOption_c = true;
			coregenesize = atoi(optarg);
		case '?':
			break;
		default:
			usage(string(argv[0]) + " [options] expfile gmtfile", options);
			return 1;
		}
	}

	if (argc - optind != 2) {
		usage(string(argv[0]) + " [options] expfile gmtfile", options);
		return 1;
	}

	argv += optind;

	//--
	// split communicator for multiple relrad's
	//--
	MPI_Comm comm;
	if (n_relrad == 0) n_relrad = 1;
	if (numProcs % n_relrad != 0) {
		cerr << "ERROR: num of MPI procs must be multiple of num of relrad's" << endl;
		return 1;
	}
	int numProcsLocal = numProcs / n_relrad;
	int color = myRank / numProcsLocal;
	int key = myRank % numProcsLocal;
	MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm);
	int myRankS, numProcsS;
	MPI_Comm_size(comm, &numProcsS);
	MPI_Comm_rank(comm, &myRankS);
	stringstream ss;
	ss << "global-rank: " << myRank << '/' << numProcs
	   << ", local-rank: " << myRankS << '/' << numProcsS << '/' << color << '\n';
	cerr << ss.str();

	//--
	// main process
	//--
	string matfile(argv[0]);

	MyMat mymat(matfile); // read expression matrix
	CoherenceBasedEEMsearch CE(mymat);

	MyMap<string, vector<string> > GeneSetMap;
	MyFunc::readGeneSetFromGmtFile(argv[1], GeneSetMap);
	CE.setGeneSets(GeneSetMap);
	if(hasOption_c){
		CE.setCoreGeneSize(coregenesize);
	}
	if(hasOption_p){
		CE.setPvalue1Cutoff(pcut);
	}
	if(hasOption_i){
		CE.setItrForPvalue2Calculation(itr);
	}
	if(hasOption_m){
		CE.setMinGeneSetSize(mingeneset);
	}
	if(hasOption_M){
		CE.setMaxGeneSetSize(maxgeneset);
	}
	if(hasOption_R){
		CE.setRelativeRadius(relrads[color], comm);
	}
	// parallelized for multiple  RelativeRadius parameters
	CE.perform(comm);

	double edctime= MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
	double stgtime= MPI_Wtime();

	//--
	// re-group procs to gather data in root procs to rank 0
	//--
	MPI_Comm_free(&comm);
	color = myRank % numProcsLocal;
	key = myRank / numProcsLocal;
	MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm);
	MPI_Comm_size(comm, &numProcsS);
	MPI_Comm_rank(comm, &myRankS);
	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "global-rank: " << myRank << '/' << numProcs
	   << ", local-rank: " << myRankS << '/' << numProcsS << '/' << color << '\n';
	cerr << ss.str();

	ExpressionModuleSet *ems = CE.getExpressionModuleSet();

	//--
	// gather result to rank 0
	//--
	if(color == 0){
		if(myRankS > 0) { /* send em */
			cerr << "sending results to rank=0" << endl;
			vector<string> ids = ems->getIds();
			vector<string>::iterator it;
			int len = ids.size();
			MPI_Send(&len, 1, MPI_INT, 0, 1234, comm);
			for(it = ids.begin(); it != ids.end(); it++) {
				ems->get(*it)->send_to(0, comm);
			}
		} else { /* recv em */
			int irank;
			for(irank = 1; irank < numProcsS; irank++) {
				cerr << "receiving results from rank=" << irank*numProcsLocal << endl;
				int len;
				MPI_Status status;
				MPI_Recv(&len, 1, MPI_INT, irank, 1234, comm, &status);
				for(int i = 0; i < len; i++) {
					ExpressionModule *em = new ExpressionModule(irank, comm);
					ems->add(em);
				}
			}
		}
	}

	double edgtime= MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
	double stotime= MPI_Wtime();

	//--
	// output
	//--
	if(myRank == 0) {
		if(hasOption_o){
			ofstream os(outfile.c_str());
			vector<string> ids = ems->getIds();
			vector<string>::iterator it;
			for(it = ids.begin(); it != ids.end(); it++) {
				os << *ems->get(*it) << endl;
			}
			os.close();
		}else{
			vector<string> ids = ems->getIds();
			vector<string>::iterator it;
			for(it = ids.begin(); it != ids.end(); it++) {
				cout << *ems->get(*it) << endl;
			}
		}
	}
	if(hasOption_l){
		stringstream ss;
		ss << logfile << '.' << myRank;
		ofstream os(ss.str().c_str());
		os << CE.getLog();
		os.close();
	}

	MPI_Comm_free(&comm);

	double edotime= MPI_Wtime();
	ss.str(""); ss.clear(stringstream::goodbit);
	ss << "main: rank=" << myRank << "/" << numProcs
	   << ", total time=" << edotime-stctime
	   << ", calculation time=" << edctime-stctime
	   << ", wait time=" << stgtime-edctime
	   << ", gather time=" << edgtime-stgtime
	   << ", wait time=" << stotime-edgtime
	   << ", output time=" << edotime-stotime << "\n";
	cerr << ss.str();
	//--
	// finalize
	//--
	MPI_Finalize();

	return 0;
}
