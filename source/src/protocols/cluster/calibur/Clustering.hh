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
#include <iosfwd>
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

#include <core/types.hh>

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
	std::vector<core::Real> calpha_vector_;
	std::vector<core::Real> signature_; //signature
	SimPDBOP mPDB;
	Stru(SimPDBOP pdb, int len, bool use_sig);
	~Stru();
	core::Real dist(core::Real x, core::Real y, core::Real z, core::Real const *zz);
	core::Real dist(core::Real x, core::Real y, core::Real z);
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
	core::Real xPercentile;
	bool autoAdjustPercentile;
	core::Real xFactor;

	std::string mInputFileName;   // file which contains all PDB filenames
	StringVec names_; // all decoy (file) names
	std::vector< StruOP > PDBs_;  // all decoy PDBs
	int n_pdbs_;   // will be set to mPDBs->size()
	int mLen;      // #residues

	core::Real THRESHOLD;  // clustering threshold. most important parameter

	// - = - = - = - = - = - = - = - = - = - = - = - = - = -
	// for auxiliary grouping

	core::Real CLU_RADIUS;  // cluster radius for auxClustering()
	std::vector< int > mCluCen; // cluster centers found using auxClustering()
	std::vector< int >** auxiliary_clusters_; // cluster elements
	std::vector< int > mCen;
	std::vector< core::Real > mD2C;    // distance from decoy in auxCluster to CluCen
	std::vector< int > num_neighbors_;   // the number of neighbors of each decoy
	core::Real bestClusMargin; // size(bestClus) -origsuze(2ndClus) /size(bestClus)
	int bestClusSize;

	// - = - = - = - = - = - = - = - = - = - = - = - = - = -
	// for clustering

	std::vector< AdjacentListOP > adjacent_lists_; // lists of all neighbors
	std::vector< core::Real > references_;   // for {lower,upper}bounds through references

	int mFinalDecoy;
	std::vector< AdjacentListOP > mFinalClusters;
	std::vector< int > remaining_list_;
	std::vector< int > remaining_list_index_;
	int mRemainingSize;

	// RMSD storage, apparently?
	// Why aren't these temporaries
	// I think for preallocation, use multiple times.
	bool spaceAllocatedForRMSD;
	core::Real *result_coords;
	// storage for rmsfit_() computation
	core::Real *coord1;
	core::Real *coord2;

	//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
	// METHODS
	//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

	Clustering(); // do nothing
	void initialize( std::string const & filename, core::Real threshold);
	void reinitialize( StringVec const &, std::vector< StruOP > const &, core::Real const threshold );
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
	core::Real get_threshold(std::vector< std::vector< core::Real > > const & nbors, int, int, int, int, core::Real, core::Real);
	void get_neighbor_list(std::vector< StruOP > const &, StringVec const &, int, std::vector< std::vector< core::Real > > & nbors );
	//core::Real ** get_neighbor_list(std::vector<Stru *> *, std::vector<Stru *> *, int);
	void get_neighbor_list(std::vector< StruOP > const &, int, core::Real &, core::Real &, std::vector< std::vector< core::Real > > & nbors );
	void estimateDist( StringVec const &, int, int, core::Real,
		core::Real &, core::Real &, core::Real &, core::Real &);
	void estimateDist(std::vector< StruOP > const &, core::Real,
		core::Real &, core::Real &, core::Real &, core::Real &);

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
	void refBound(int i, int j, core::Real& lower, core::Real& upper);
	//bool find(int which, std::vector<int> *elements); // too slow

	// - = - = - = - = - = - = - = - = - = - = - = -

	void realignDecoys(int ref);
	void superimposeAndReplace(std::vector<core::Real> & coor1, std::vector<core::Real> & coor2);
	core::Real eucD(int i, int j);

	// - = - = - = - = - = - = - = - = - = - = - = -

	// methods for rmsd computation
	core::Real estD( int i, int j );
	core::Real estD( Stru const & a, Stru const & b );
	core::Real trueD( Stru const & a, Stru const & b );
	core::Real trueD( int i, int j );
	// storage for RMSD() computation
	void allocateSpaceForRMSD(int len);
};

}
}
}

#endif
