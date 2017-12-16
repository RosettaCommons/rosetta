// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file external/calibur/Clustering.hh
/// @author SC Li & YK Ng (kalngyk@gmail.com)

#ifndef external_calibur_Clustering_HH
#define external_calibur_Clustering_HH

#include <vector>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <assert.h>
#ifndef PYROSETTA
#include <time.h>
#endif
#if !defined(__WIN32__) && !defined(WIN32)
#include <sys/resource.h>
#endif

#include <utility/pointer/owning_ptr.hh>
#include <protocols/cluster/calibur/pdb_util.hh>
#include <protocols/cluster/calibur/SimPDB.hh>
#include <protocols/cluster/calibur/PreloadedPDB.hh>

namespace protocols {
namespace cluster {
namespace calibur {

class AdjacentList; // fwd
typedef utility::pointer::shared_ptr< AdjacentList > AdjacentListOP;

enum EST_THRESHOLD_MODE {
	PERCENT_EDGES,   // % quantile pairwise distance (default)
	MIN_AVG_DIST_BASED, // t = a * (min avg dist) + b
	MOST_FREQ_BASED, // use (2/3)*(most frequently occuring distance)
	ROSETTA,   // ROSETTA's method
	SAMPLED_ROSETTA, // ROSETTA's method using sampled decoys
	USER_SPECIFIED   // user supplied
};

// Instead of all these statics, we're going to have clustering objects actually
// own their configuration.

class Stru
{
public:
	std::vector<double> calpha_vector_;
	std::vector<double> signature_; //signature
	SimPDBOP mPDB;
	Stru(SimPDBOP pdb, int len, bool use_sig);
	~Stru();
	double dist(double x, double y, double z, double const *zz);
	double dist(double x, double y, double z);
};

typedef utility::pointer::shared_ptr< Stru > StruOP;

class Clustering
{
public:
	bool pref_use_scud;
	bool pref_use_sig;
	bool pref_ref_use_first_few_decoys;

	EST_THRESHOLD_MODE EST_THRESHOLD;
	bool FILTER_MODE;
	double xPercentile;
	bool autoAdjustPercentile;
	double xFactor;

	std::string mInputFileName;   // file which contains all PDB filenames
	StringVec names_; // all decoy (file) names
	std::vector< StruOP > PDBs_;  // all decoy PDBs
	int n_pdbs_;   // will be set to mPDBs->size()
	int mLen;      // #residues

	double THRESHOLD;  // clustering threshold. most important parameter

	// - = - = - = - = - = - = - = - = - = - = - = - = - = -
	// for auxiliary grouping

	double CLU_RADIUS;  // cluster radius for auxClustering()
	std::vector< int > mCluCen; // cluster centers found using auxClustering()
	std::vector< int >** auxiliary_clusters_; // cluster elements
	std::vector< int > mCen;
	std::vector< double > mD2C;    // distance from decoy in auxCluster to CluCen
	std::vector< int > num_neighbors_;   // the number of neighbors of each decoy
	double bestClusMargin; // size(bestClus) -origsuze(2ndClus) /size(bestClus)
	int bestClusSize;

	// - = - = - = - = - = - = - = - = - = - = - = - = - = -
	// for clustering

	std::vector< AdjacentListOP > adjacent_lists_; // lists of all neighbors
	std::vector< double > references_;   // for {lower,upper}bounds through references

	int mFinalDecoy;
	std::vector< AdjacentListOP > mFinalClusters;
	std::vector< int > remaining_list_;
	std::vector< int > remaining_list_index_;
	int mRemainingSize;

	// RMSD storage, apparently?
	// Why aren't these temporaries
	// I think for preallocation, use multiple times.
	bool spaceAllocatedForRMSD;
	double *result_coords;
	// storage for rmsfit_() computation
	double *coord1;
	double *coord2;

	//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
	// METHODS
	//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

	Clustering(); // do nothing
	void initialize( std::string const & filename, double threshold);
	void reinitialize( StringVec const &, std::vector< StruOP > const &, double const threshold );
	void cluster();
	void showClusters(int);
	void getPDBs( StringVec &, std::vector< StruOP > &, std::vector<int> const &, int);


	//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

	// for reading decoys from input files
	void readDecoyNames();
	void readDecoys( StringVec const &, std::vector< StruOP > const &);
	std::vector< StruOP > readDecoys( StringVec const &);
	//void refilterDecoys(std::vector<char *>*, std::vector<Stru *>*);

	// - = - = - = - = - = - = - = - = - = - = - = -
	// for finding threshold

	void getThresholdAndDecoys();
	double get_threshold(std::vector< std::vector< double > > const & nbors, int, int, int, int, double, double);
	void get_neighbor_list(std::vector< StruOP > const &, StringVec const &, int, std::vector< std::vector< double > > & nbors );
	//double ** get_neighbor_list(std::vector<Stru *> *, std::vector<Stru *> *, int);
	void get_neighbor_list(std::vector< StruOP > const &, int, double &, double &, std::vector< std::vector< double > > & nbors );
	void estimateDist( StringVec const &, int, int, double,
		double &, double &, double &, double &);
	void estimateDist(std::vector< StruOP > const &, double,
		double &, double &, double &, double &);

	std::vector< std::string > getRandomDecoyNames( StringVec const &, int );
	void destroyRandomDecoys( StringVec const &, std::vector< StruOP > const &);


	// - = - = - = - = - = - = - = - = - = - = - = -
	// for clustering

	void auxClustering();
	void buildAdjacentLists();
	void listAdjacentLists();

	void findLargestClusters();
	int findDecoyWithMostNeighbors();
	void removeDecoys(AdjacentList & adj);

	// - = - = - = - = - = - = - = - = - = - = - = -

	void initRef(int* index);
	void refBound(int i, int j, double& lower, double& upper);
	//bool find(int which, std::vector<int> *elements); // too slow

	// - = - = - = - = - = - = - = - = - = - = - = -

	void realignDecoys(int ref);
	void superimposeAndReplace(std::vector<double> & coor1, std::vector<double> & coor2);
	double eucD(int i, int j);

	// - = - = - = - = - = - = - = - = - = - = - = -

	// methods for rmsd computation
	double estD( int i, int j );
	double estD( Stru const & a, Stru const & b );
	double trueD( Stru const & a, Stru const & b );
	double trueD( int i, int j );
	// storage for RMSD() computation
	void allocateSpaceForRMSD(int len);
};

}
}
}

#endif
