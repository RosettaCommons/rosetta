// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/kalngyk/InitCluster.cc
/// @author SC Li & YK Ng (kalngyk@gmail.com)

//#define _DEBUG_NBORLIST_

#include <vector>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <assert.h>
#include <time.h>
#ifndef __WIN32__
#include <sys/resource.h>
#endif

//#include <protocols/cluster/calibur/SimpPDB.hh>
#include <numeric/random/random.hh>
#include <protocols/cluster/calibur/rmsd.hh>
#include <protocols/cluster/calibur/Clustering.hh>
#include <protocols/cluster/calibur/AdjacentList.hh>
#include <protocols/cluster/calibur/InitCluster.hh>

#include <utility/exit.hh> // for runtime_assert

namespace protocols {
namespace cluster {
namespace calibur {

ADJ_LIST_MODE AdjacentList::list_mode_;

// _MATRIX_MODE_LIMIT_ is determined by the RAM of the system
#define _MATRIX_MODE_LIMIT_ 13000

#define _SHOW_PERCENTAGE_COMPLETE_

#define _OVER_RMSD_ 9999999

#define REFERENCE_SIZE 6
#define RANDOM_DECOY_SIZE_FOR_FILTERING 101
#define RANDOM_DECOY_SIZE_FOR_THRESHOLD 101
#define NUM_TRIALS_FOR_THRESHOLD 16
#define DEFAULT_PERCENTILE_FOR_THRESHOLD 10
#define MAX_PERCENTILE_FOR_THRESHOLD 50
#define MIN_PERCENTILE_FOR_THRESHOLD 3

//============================== Structure ====================================
Stru::Stru(SimPDBOP pdb, int len, bool use_sig):
	calpha_vector_( pdb->calpha_vector_ ),
	mPDB( pdb )
{
	if ( use_sig ) {
		signature_.resize( len );
		//double cx = 0;
		//double cy = 0;
		//double cz = 0;
		//compute the centroid
		/*for (int i=0; i < len; i++)
		{
		cx += calpha_vector_[3*i];
		cy += calpha_vector_[3*i+1];
		cz += calpha_vector_[3*i+2];
		}*/
		//cx /= len;
		//cy /= len;
		//cz /= len;
		for ( int i=0; i<len; i++ ) {
			//signature_[i] = dist(cx, cy, cz, calpha_vector_+3*i);
			// for the new SimpPDB.cc (it became slower so I reverted back)
			signature_[i] = dist(calpha_vector_[3*i], calpha_vector_[3*i+1], calpha_vector_[3*i+2]);
		}
	}
}

Stru::~Stru() {}

double
Stru::dist(double x, double y, double z, double *zz)
{
	double xd=x-zz[0];
	double yd=y-zz[1];
	double zd=z-zz[2];
	return sqrt(xd*xd+yd*yd+zd*zd);
}

double
Stru::dist(double x, double y, double z)
{
	return sqrt(x*x+y*y+z*z);
}
//============================= Structure =====================================

//============================ Main object ====================================
Clustering::Clustering() :
	pref_use_scud( false ),
	pref_use_sig( false ),
	pref_ref_use_first_few_decoys( false ),
	EST_THRESHOLD( ROSETTA ), // an arbitrary default
	FILTER_MODE( false ),
	xPercentile( 0 ),
	autoAdjustPercentile( false ),
	xFactor( 0 ),
	mInputFileName( "" ),
	//names_
	//PDBs_
	n_pdbs_( 0 ),
	mLen( 0 ),
	THRESHOLD( 0 ),
	CLU_RADIUS( 0 ),
	// mCluCen
	auxiliary_clusters_( nullptr ),
	//mCen( nullptr ),
	//mD2C( nullptr ),
	//num_neighbors_( nullptr ),
	bestClusMargin( 0 ),
	bestClusSize( 0 ),
	//adjacent_lists_
	//references_( nullptr ),
	mFinalDecoy( 0 ),
	//mFinalClusters
	//remaining_list_
	//remaining_list_index_
	mRemainingSize( 0 ),
	spaceAllocatedForRMSD( false ),
	result_coords( nullptr ),
	coord1( nullptr ),
	coord2( nullptr )
{
	mLen = 0;
	spaceAllocatedForRMSD = false;
	bestClusMargin = 1.; // should be a value that will not trigger re-cluster
#ifdef _SPICKER_SAMPLING_
	FILTER_MODE = false;
#endif
	std::cout.setf(std::ios::fixed);
	std::cout << std::setprecision(5);
}

/**
* Main aim:
* (1) read decoys into names_ and PDBs_.
* (2) decide THRESHOLD, CLUS_RADIUS
* (3) initialize data structure needed for auxiliary grouping, and
* (4) ...for clustering.
*/
void
Clustering::initialize( std::string const & filename, double const threshold )
{
	mInputFileName = filename;
	THRESHOLD = threshold;

	// Read decoys and get threshold - = - = - = - = - = - = -

	if ( THRESHOLD != -1 ) {
		EST_THRESHOLD = USER_SPECIFIED;
	}

	std::cout << "Filtering " << (FILTER_MODE? "on": "off") << std::endl;
	std::cout << "Signature mode " << (pref_use_sig? "on": "off") << std::endl;
	std::cout << "Using chains ";
	for ( char * c = SimPDB::chains; *c; c++ ) {
		std::cout << "'" << *c << "'" << (*(c+1)=='\0'? "": ", ");
	}
	std::cout << " in PDB files" << std::endl;

	getThresholdAndDecoys();

	CLU_RADIUS = THRESHOLD / 2.0 - 0.00001;

	std::cout << "Initialized " << n_pdbs_ << " decoys." << std::endl;

	// Realign decoys if needed - = - = - = - = - = - = - = -

	if ( THRESHOLD > 2.5 ) pref_use_scud = true;

	std::cout << "Use of SCUD's bound " << (pref_use_scud? "on": "off") << std::endl;

	if ( pref_use_scud ) realignDecoys(0);

	if ( EST_THRESHOLD == USER_SPECIFIED ) {
		std::cout << "Using user specified threshold: " << THRESHOLD << std::endl;
	}
	// Initialize auxiliary decoys grouping - = - = - = - = -

	// mCluCen contains the cluster centers. Each cluster center is referenced
	// by an integer x, which is the decoy's index in PDBs_.
	// Elements of the cluster is in auxiliary_clusters_[x];
	// Since most decoys are going to end up not being cluster centers, for
	// most i, auxiliary_clusters_[i] will be NULL.
	mCluCen.clear(); // to be certain = new std::vector<int>(0);
	auxiliary_clusters_ = new std::vector<int>*[n_pdbs_];
	for ( int i=0; i < n_pdbs_; i++ ) {
		auxiliary_clusters_[i] = nullptr;
	}
	mD2C.resize( n_pdbs_ );
	mCen.resize( n_pdbs_ );

	// Initialize clustering - = - = - = - = - = - = - = - = -

	adjacent_lists_.resize(n_pdbs_);

#ifdef _ADD_LITE_MODE_
	if (AdjacentList::list_mode_ != LITE) {
#endif
	// auto-switch between LIST and MATRIX mode
	if ( n_pdbs_ > _MATRIX_MODE_LIMIT_ ) {
		std::cout << "Using LIST mode" << std::endl;
		AdjacentList::list_mode_ = LIST;
	} else {
		std::cout << "Using MATRIX mode" << std::endl;
		AdjacentList::list_mode_ = MATRIX;
	}
#ifdef _ADD_LITE_MODE_
	}
	else {
		std::cout << "Using LITE mode" << std::endl;
		AdjacentList::list_mode_ = LITE;
	}
#endif

	for ( int i = 0; i < n_pdbs_; i++ ) {
		if ( AdjacentList::list_mode_ == MATRIX ) {
			adjacent_lists_[i] = AdjacentListOP( new AdjacentList(i, n_pdbs_) );
		} else {
			adjacent_lists_[i] = AdjacentListOP( new AdjacentList() );
			adjacent_lists_[i]->which_ = i;
		}
	}
}


/**
* Useful for when we want to redo clustering with a set of already
* initialized PDBs.
*/
void
Clustering::reinitialize(
	std::vector< std::string > const & nNames,
	std::vector< StruOP > const & nPDBs,
	double const threshold
) {
	THRESHOLD = threshold;

	// Read decoys and get threshold - = - = - = - = - = - = -

	CLU_RADIUS = THRESHOLD / 2.0 - 0.00001;
	names_ = nNames;
	for ( auto const & n : nPDBs ) PDBs_.push_back( n );
	//PDBs_ = nPDBs;
	int _n_pdbs_ = n_pdbs_;
	n_pdbs_ = PDBs_.size();

	// Initialize auxiliary decoys grouping - = - = - = - = -

	// remove old groupings
	for ( int i=0; i < _n_pdbs_; i++ ) {
		if ( auxiliary_clusters_[i] ) {
			delete auxiliary_clusters_[i];
		}
	}
	delete [] auxiliary_clusters_;

	// initialize new ones
	mCluCen.clear(); // = new std::vector<int>(0);
	auxiliary_clusters_ = new std::vector<int>*[n_pdbs_];
	for ( int i=0; i < n_pdbs_; i++ ) {
		auxiliary_clusters_[i] = NULL;
	}
	mD2C.resize( n_pdbs_ );
	mCen.resize( n_pdbs_ );

	// Initialize clustering - = - = - = - = - = - = - = - = -

	// delete old...
	adjacent_lists_.clear();

	// ...create new
	adjacent_lists_.resize(n_pdbs_);
	for ( int i=0; i < n_pdbs_; i++ ) {
		if ( AdjacentList::list_mode_ == MATRIX ) {
			adjacent_lists_[i] = AdjacentListOP(new AdjacentList(i, n_pdbs_));
		} else {
			adjacent_lists_[i] = AdjacentListOP(new AdjacentList());
			adjacent_lists_[i]->which_ = i;
		}
	}

	// Delete cluster data - = - = - = - = - = - = - = - = -

	num_neighbors_.clear();
}


/**
* Get these: names_, PDBs_, mLen
*/
void
Clustering::getThresholdAndDecoys()
{
	readDecoyNames(); // results in names_

	// decide min max target cluster sizes - = - = - = -

#ifdef _USE_FIX_CLUSTER_SIZES_
	// good for when comparing with cluster_info_silent (ROSETTA)
	unsigned int minClusterSize = 15;
	unsigned int maxClusterSize = 75;
	unsigned int targetClusterSize = 45;
#else
	unsigned int minClusterSize = names_.size()/100;
	unsigned int maxClusterSize = names_.size()/15; // INFLUENTIAL
	unsigned int targetClusterSize = (maxClusterSize + minClusterSize)/2; // INFLUENTIAL
#endif
	if ( maxClusterSize > (names_.size()-1) ) {
		maxClusterSize = names_.size()-1;
	}
	if ( targetClusterSize > (names_.size()-1) ) {
		targetClusterSize = names_.size()-1;
	}

	// decide the min max thresholds - = - = - = - = -

	double minDist, maxDist, mostFreqDist, xPercentileDist;

	estimateDist(names_,
		NUM_TRIALS_FOR_THRESHOLD,
		names_.size() > 2*RANDOM_DECOY_SIZE_FOR_THRESHOLD?
		RANDOM_DECOY_SIZE_FOR_THRESHOLD: (names_.size()/2),
		xPercentile,
		minDist,
		maxDist,
		mostFreqDist,
		xPercentileDist);

	if ( EST_THRESHOLD == MOST_FREQ_BASED ) {
		THRESHOLD = minDist + xFactor * (mostFreqDist-minDist) ;
		std::cout << "Finding threshold using most frequent distance" << std::endl;
		std::cout << "Threshold = " << THRESHOLD << "  ( " << minDist
			<< " + " << xFactor << "(" << mostFreqDist
			<< " - " << minDist << ") )" << std::endl;
	}

	if ( EST_THRESHOLD == PERCENT_EDGES ) {
		THRESHOLD = xPercentileDist;
		std::cout << "Finding threshold to keep " << xPercentile << "% edges" << std::endl;
		std::cout << "Threshold = " << THRESHOLD << std::endl;
	}

	if ( EST_THRESHOLD == MIN_AVG_DIST_BASED ) {
		// Uses a set of randomly chosen decoys to estimate the threshold
		// as in cluster_info_silent (i.e. rosetta)
		//
		// For each sampled decoy, we get its average distance to other decoys
		// Then, take the least of such average distance, m
		// Threshold is then computed as "a*m + b", where a and b are some
		// constant we decide by trial and error.
		//
		std::cout << "Finding threshold using average distance. "
			<< "This may take a while..." << std::endl;
		clock_t start = clock();

		int numDecoys = names_.size() > 101? 101: names_.size();
		std::vector< std::string > names = getRandomDecoyNames(names_, numDecoys);
		std::vector< StruOP > decoys = readDecoys(names);

		std::vector< std::vector< double > > nbors;
		get_neighbor_list(decoys, names_, names_.size()-1, nbors);

		destroyRandomDecoys(names, decoys);

		// find the minimum of average distances
		double min_avg_dist = _OVER_RMSD_;
		for ( int i=0; i < numDecoys; i++ ) {
			double sum_of_dist = 0;
			for ( unsigned int j=0; j < names_.size(); j++ ) {
				sum_of_dist += nbors[i][j];
			}
			double avg_dist = sum_of_dist / names_.size();
			if ( avg_dist < min_avg_dist ) {
				min_avg_dist = avg_dist;
			}
		}

		THRESHOLD = xFactor * min_avg_dist + minDist;

		double elapsed = (clock() - start)/(double)CLOCKS_PER_SEC;
		std::cout << "Minimum average distance = " << min_avg_dist
			<< ". Found in " << elapsed << " s" << std::endl;
		std::cout << "Threshold = " << THRESHOLD << "  (" << minDist
			<< " + " << xFactor << " * " << min_avg_dist << ")" << std::endl;
	}

	if ( EST_THRESHOLD == SAMPLED_ROSETTA ) {
		// Uses a set of randomly chosen decoys to estimate the threshold
		// as in cluster_info_silent (i.e. rosetta)
		//
		// The belief is that the largest cluster should contain a set X
		// number of decoys, within a threshold range (t1, t2).
		// This is what ROSETTA does:
		// 1. Build sorted lists of decoys, that is
		//   decoy1: neighbor1, neighbor2, ..., neighborX
		//   decoy2: neighbor1, neighbor2, ..., neighborX
		//    ...
		// 2. If there is a list such that X neighbors can be contained
		// within a threshold t, t1 <= t <= t2, then we are done.
		//
		// In the following, only the sampled decoys get their lists built.
		//
		// Required parameters:
		//   mRandomDecoys
		//   {min,max,target}ClusterSize
		//   {min,max}Threshold
		//
		std::cout << "Finding threshold using random decoys (max cluster size = "
			<< maxClusterSize << ")" << std::endl;
		clock_t start = clock();

		int numDecoys = names_.size() > 101? 101: names_.size();
		std::vector<std::string> names = getRandomDecoyNames(names_, numDecoys);
		std::vector< StruOP > decoys = readDecoys(names);

		std::vector< std::vector< double > > nbors;
		get_neighbor_list(decoys, names_, maxClusterSize, nbors);

		destroyRandomDecoys(names, decoys);

		double minThreshold = minDist;
		double maxThreshold = (minDist + maxDist)/2;

		THRESHOLD = get_threshold(nbors,
			numDecoys,   // # of lists to build
			maxClusterSize, // length of list
			minClusterSize,
			targetClusterSize,
			minThreshold,
			maxThreshold );

		double elapsed = (clock() - start)/(double)CLOCKS_PER_SEC;
		std::cout << "Threshold = " << THRESHOLD
			<< ". Found in " << elapsed << " s" << std::endl;
	}

	// read decoys - = - = - = - = - = - = - = - = -

	std::vector<std::string> randNames;
	std::vector< StruOP > randDecoys;
	if ( FILTER_MODE ) {
		int numDecoys = names_.size() > 2*RANDOM_DECOY_SIZE_FOR_FILTERING?
			RANDOM_DECOY_SIZE_FOR_FILTERING: (names_.size()/2);
		randNames = getRandomDecoyNames(names_, numDecoys);
		randDecoys = readDecoys(randNames);
	}

	clock_t start = clock();
	readDecoys(randNames, randDecoys); // results in PDBs_
	double elapsed = (clock() - start)/(double)CLOCKS_PER_SEC;
	std::cout << "Decoys read in " << elapsed << " s" << std::endl;

	if ( FILTER_MODE ) {
		destroyRandomDecoys(randNames, randDecoys);
	}

	// find threshold using ROSETTA mode - = - = - = - = - = - = - = - = -
	// (this has to be done with the full decoys) - = - = - = - = - = - =

	if ( EST_THRESHOLD == ROSETTA ) {
		// This is exactly ROSETTA's way of getting threshold.
		// We add this in for comparison with the output of ROSETTA.
		//
		// Note that targetClusterSize is guaranteed since we pass the
		// minimum and maximum distance between decoys --obtained from
		// get_neighbor_list-- to getThreshold().
		//
		// Required parameters:
		//   PDBs_
		//   {min,max,target}ClusterSize
		//
		//   {min,max}Threshold are automatically decided accurately
		//   This is so that targetClusterSize will always pass
		//
		std::cout << "Finding threshold with ROSETTA's method (max cluster size = "
			<< maxClusterSize << ")" << std::endl;
		clock_t start = clock();

		double minThreshold, maxThreshold;

		std::vector< std::vector< double > > nbors;
		get_neighbor_list(PDBs_,
			maxClusterSize,
			minThreshold,
			maxThreshold,
			nbors );

		THRESHOLD = get_threshold(nbors,
			PDBs_.size()-1,
			maxClusterSize,
			minClusterSize,
			targetClusterSize,
			minThreshold,
			maxThreshold);

		double elapsed = (clock() - start)/(double)CLOCKS_PER_SEC;
		std::cout << "Threshold = " << THRESHOLD
			<< ". Found in " << elapsed << " s" << std::endl;
	}
}


/**
* Repeatedly estimate distances for numTrials trials, and get the average
* values of
* 1. min of pairwise distance,
* 2. max of pairwise distance,
* 3. most frequent distance,
* 4. distance t at which xPercent of pairwise distances are below t,
* from all trials. The number of random decoys to use in each trial
* is stated in randDecoySize.
*/
void
Clustering::estimateDist( std::vector< std::string > const & allNames,
	int numTrials,
	int randDecoySize,
	double xPercent,
	double & minDist,
	double & maxDist,
	double & mostFreqDist,
	double & xPercentileDist
) {
	StringVec randNames;
	std::vector< StruOP > randDecoys;

	std::cout << "Estimating threshold range...";

	randNames = getRandomDecoyNames(allNames, randDecoySize);
	randDecoys = readDecoys(randNames);

	estimateDist(randDecoys,
		xPercent,
		minDist,
		maxDist,
		mostFreqDist,
		xPercentileDist);

	destroyRandomDecoys(randNames, randDecoys);

#ifdef _ANALYZE_RANDOM_SAMPLES_
	std::cout << std::endl
		 << "sampling 1: dist [" << *minDist << "," << *maxDist
		 << "], most freq = " << *mostFreqDist << ", "
		 << xPercent << " percentile = " << *xPercentileDist
		 << std::endl;
#elif defined(_SHOW_PERCENTAGE_COMPLETE_)
	printf("\rEstimating threshold range... %4.1f%%\r", 100. /numTrials);
	fflush(stdout);
#endif

	std::vector< double > maxDists( numTrials - 1 );
	std::vector< double > minDists( numTrials - 1 );
	std::vector< double > mostFreqDists( numTrials - 1 );
	std::vector< double > xPercentileDists( numTrials - 1 );
	for ( int i=0; i < numTrials-1; i++ ) {
		randNames = getRandomDecoyNames(allNames, randDecoySize);
		randDecoys = readDecoys(randNames);
		estimateDist(randDecoys,
			xPercent,
			minDists[i],
			maxDists[i],
			mostFreqDists[i],
			xPercentileDists[i]);

		destroyRandomDecoys(randNames, randDecoys);

#ifdef _ANALYZE_RANDOM_SAMPLES_
		std::cout << "sampling " << (i+2) << ": dist [" << minDists[i] << ","
			 << maxDists[i] << "], most freq = "
			 << mostFreqDists[i] << ", " << xPercent
			 << " percentile = " << xPercentileDists[i] << std::endl;
#elif defined(_SHOW_PERCENTAGE_COMPLETE_)
		printf("Estimating threshold range... %4.1f%%\r", 100.*(i+2)/numTrials);
		fflush(stdout);
#endif
	}

#if !defined(_ANALYZE_RANDOM_SAMPLES_) && !defined(_SHOW_PERCENTAGE_COMPLETE_)
	std::cout << " done" << std::endl;
#elif defined(_SHOW_PERCENTAGE_COMPLETE_)
	std::cout << "Estimating threshold range... complete" << std::endl;
#endif

	for ( int i=0; i < numTrials-1; i++ ) {
		if ( minDists[i] < minDist ) {
			minDist = minDists[i];
		}
		if ( maxDists[i] > maxDist ) {
			maxDist = maxDists[i];
		}
		mostFreqDist += mostFreqDists[i];
		xPercentileDist += xPercentileDists[i];
	}
	mostFreqDist /= numTrials;
	xPercentileDist /= numTrials;

	std::cout << "Estimated threshold range = ["
		<< minDist << ", " << maxDist << "]" << std::endl
		<< "Most frequent distance = " << mostFreqDist
		<< std::endl
		<< xPercent << " percent distances below " << xPercentileDist
		<< std::endl;
}


/**
* Read from the input file the names of decoys
*/
void
Clustering::readDecoyNames()
{
	switch (filetype(mInputFileName)) {
	PreloadedPDB * pdbs;
	case SILENT_FILE :
		pdbs = new PreloadedPDB();
		pdbs->loadSilentFile(mInputFileName); // Preload PDBs from silent file
		SimPDB::preloadedPDB = pdbs; // Attach the preloaded PDBs to SimPDB
		SimPDB::preloadPDB = true;
		names_ = *pdbs->names_;
		std::cout << "Read " << names_.size() << " decoy names" << std::endl;
		return;
	case PDB_LIST :
		if ( SimPDB::preloadPDB ) {
			unsigned int numDecoys = num_lines_in_file(mInputFileName);
			if ( numDecoys > PreloadedPDB::ADVISED_THRESHOLD ) {
				std::cout << "More than " << PreloadedPDB::ADVISED_THRESHOLD
					<< " decoys: not preloading PDBs" << std::endl;
				SimPDB::preloadPDB = false;
			}
		}
		if ( !SimPDB::preloadPDB ) { // either switched off by user, or by above
			break;
		}
		pdbs = new PreloadedPDB();
		pdbs->loadPDBFromList(mInputFileName); // Preload PDBs from pdb list
		SimPDB::preloadedPDB = pdbs; // Attach the preloaded PDBs to SimPDB
		std::cout << "Before names_ assignment" << std::endl;
		names_ = *pdbs->names_;
		std::cout << "Read " << names_.size() << " decoy names" << std::endl;
		return;
	default :
		break;
	}

	std::ifstream input(mInputFileName);
	if ( !input ) {
		std::cerr << "Cannot find file \"" << mInputFileName << "\"" << std::endl;
		exit(0);
	}

	char buf[400];
	names_.clear(); // porbably already empty
	while ( !input.eof() ) {
		input.getline(buf, 400);
		char* token = strtok(buf, " ");
		if ( !token ) continue;
		char* name = new char[strlen(token)+1];
		strcpy(name, token);
		names_.push_back(name);
	}
	input.close();
	std::cout << "Read " << names_.size() << " decoy names" << std::endl;
}


void
Clustering::allocateSpaceForRMSD( int len ) {
	if ( spaceAllocatedForRMSD ) {
		return;
	}
	result_coords = new double[3*len];
	spaceAllocatedForRMSD = true;
}


/**
* Read from the input files the decoys and return them
*/
std::vector< StruOP >
Clustering::readDecoys( std::vector<std::string> const & decoynames ) {
	SimPDBOP firstPDB( new SimPDB( decoynames[0] ) );
	mLen = firstPDB->num_residue_;

	allocateSpaceForRMSD(mLen);

	std::vector< StruOP > decoys( decoynames.size(), nullptr );
	decoys[0] = StruOP( new Stru( firstPDB, mLen, pref_use_sig ) );//std::shared_ptr<Stru>( new Stru(firstPDB, mLen) );

	for ( unsigned int i=1; i < decoynames.size(); i++ ) {
		decoys[i] = StruOP( new Stru( SimPDBOP( new SimPDB( decoynames[i], mLen ) ), mLen, pref_use_sig ) );
	}
	return decoys;
}


/**
* Read from the input files all the decoys, filtering when necessary
*/
void
Clustering::readDecoys(
	std::vector< std::string > const & randonames_,
	std::vector< StruOP > const & randomDecoys
) {
	std::vector< std::string > newNames;
	std::vector< StruOP > newDecoys;
#ifdef _SPICKER_SAMPLING_
	// Spicker samples decoys at a fixed interval delta
	// such that exactly 13000 decoys are sampled
	std::cout << "Sampling 13000 decoys" << std::endl;
	double delta = names_.size()/(double)13000;
	int sampled_decoy_id = 1;
	if (delta < 1) delta = 1;
#endif

	std::cout << "Reading decoys with signature mode " << (pref_use_sig?"on":"off")
		<< " and filter mode " << (FILTER_MODE? "on": "off") << std::endl;

	int size = names_.size();
	for ( int i=0; i < size; i++ ) {
#ifdef _SPICKER_SAMPLING_
		if ((i+1) < sampled_decoy_id * delta) {
			continue;
		}
		sampled_decoy_id++;
#endif
#ifdef _SHOW_PERCENTAGE_COMPLETE_
		printf("Read %4.1f%%\r", 100.*i/size);
		fflush(stdout);
#endif
		// read in the decoy's PDB
		std::string dName = names_[i];
		StruOP s;
		if ( mLen == 0 ) {
			SimPDBOP aPDB( new SimPDB(dName));
			mLen = aPDB->num_residue_;
			s = StruOP( new Stru(aPDB, mLen, pref_use_sig));
			allocateSpaceForRMSD(mLen);
		} else {
			s = StruOP( new Stru( SimPDBOP( new SimPDB(dName, mLen)), mLen, pref_use_sig) );
		}
		bool isOutlier = false;
		if ( FILTER_MODE ) { // then we shall decide whether to include s
			int randomDecoysSize = randomDecoys.size();
			isOutlier = true;
			for ( int j=0; j < randomDecoysSize; j++ ) {
				if ( dName == randonames_[j] ) {
					continue;
				}
				if ( pref_use_sig && estD(*s,*randomDecoys[j]) > 2*THRESHOLD ) {
					continue;
				}
				if ( trueD(*s,*randomDecoys[j]) > 2*THRESHOLD ) {
					continue;
				}
				// s is near to a random decoy
				isOutlier = false;
				break;
			}
		}
		if ( !isOutlier ) {
			newNames.push_back(dName);
			newDecoys.push_back( StruOP( new Stru( *s ) ) );
		}
	}
	std::cout << "Read " << newNames.size() << " decoys.";
	if ( FILTER_MODE ) {
		std::cout << " Filtered " << (names_.size() - newNames.size())
			<< " outlier decoys.";
	}
	std::cout << std::endl;
	names_ = newNames;
	PDBs_  = newDecoys;
	n_pdbs_ = PDBs_.size();
	if ( !n_pdbs_ ) {
		std::cout << "No decoy read. Exiting..." << std::endl;
		exit(0);
	}
}


/**
* Main call. Find all neighbors within threshold.
*/
void
Clustering::cluster()
{
	clock_t start;
	double elapsed;

	std::cout << "Auxiliary clustering...";
#ifdef _SHOW_PERCENTAGE_COMPLETE_
	printf("\r");
#endif
	start = clock();
	initRef(NULL);
	auxClustering(); // Cluster to speed-up computation of neighbors
	elapsed = (clock() - start)/(double)CLOCKS_PER_SEC;
#ifdef _SHOW_PERCENTAGE_COMPLETE_
	std::cout << "\nAuxiliaryClustering...";
#endif
	std::cout << " completed in " << elapsed << " s" << std::endl;

	std::cout << "Finding decoys neighbors...";
	//start = clock();
#ifndef __WIN32__
#ifndef PYROSETTA
	_get_elapsed(1);
#endif
#endif
	buildAdjacentLists(); // Find all the neighbors for each decoy
	//elapsed = (clock() - start)/(double)CLOCKS_PER_SEC;
#ifndef __WIN32__
#ifndef PYROSETTA
	elapsed = _get_elapsed(0);
	std::cout << "\nFinding all neighbors completed in " << elapsed << " s" << std::endl;
#else
	std::cout << " completed" << std::endl;
#endif
#else
	std::cout << " completed" << std::endl;
#endif

	//listAdjacentLists();

	std::cout << "Finding and removing largest clusters...";
	start = clock();
	findLargestClusters(); // Find largest clusters and remove. Recurse.
	elapsed = (clock() - start)/(double)CLOCKS_PER_SEC;
	std::cout << " completed in " << elapsed << " s" << std::endl;

#ifdef _ADD_LITE_MODE_
	return;
#endif

	/**
	* If the first cluster is not much bigger than the second cluster,
	* two independent clusters may exist. Hence if this difference is
	* below a certain margin we may want to do something about it.
	*
	* We want an almost constant margin so we compare the sizes of the two
	* clusters prior to the removal of decoys from the second. If we use
	* the size of the clusters after the removal, we will have to consider
	* the fact that: at higher thresholds the overlap between clusters will
	* will increase, and the margin will increase likewise.
	*
	* That the second cluster is not originally large due to its similarity
	* to the first cluster (i.e. they are two almost equal center from the
	* same cluster) is more or less guaranteed from the fact that, if this
	* is the case, then most of its decoys would've been removed during the
	* clustering and it would have lost out to some other cluster that
	* has less overlap with the best cluster. Hence the presence of at least
	* one such "other cluster* should be ensured.
	*/
	AdjacentList const & bestClus = *mFinalClusters[0];
	bestClusSize = bestClus.num_neighbors_;
	if ( mFinalClusters.size() > 1 ) {
		AdjacentList const & nextClus = *mFinalClusters[1];
		int nextSize = num_neighbors_[nextClus.which_];
		bestClusMargin = (bestClusSize - nextSize) /(double) bestClusSize;
		std::cout << std::endl << std::endl<<"Largest 2 cluster centers: "
			<< names_[bestClus.which_] << "(" << bestClusSize << "), "
			<< names_[nextClus.which_] << "(" << nextSize
			<< "). Margin = " << (bestClusMargin*100) << "%" << std::endl;
	}
}


void
Clustering::showClusters(int numClus)
{
#ifdef _ADD_LITE_MODE_
	if (AdjacentList::list_mode_ == LITE) {
		std::cout << "Best decoy: " << names_[mFinalDecoy] << "\t"
			 << adjacent_lists_[mFinalDecoy]->num_neighbors_ << std::endl;
		return;
	}
#endif

	int size = mFinalClusters.size();
	if ( size > numClus ) size = numClus;
	for ( int i=0; i<size; i++ ) {

		// For each cluster in the final clusters...

		// (Note: it's possible for the proper number of clusters, numClus, to be
		// smaller. There could be empty ones at the end of mFinalClusters, I guess,
		// which we skip.)

		AdjacentList const & clu = *mFinalClusters[i];
		//std::cout << clu->which_ << "'" << names_[clu->which_]
		std::cout << std::endl << "cluster = " << i << "; center = " << names_[clu.which_]
			<< "; n_decoy_members = " << clu.num_neighbors_ << "; members = ";
		uint num_neighbors_this_cluster = clu.num_neighbors_;
		std::vector<int> const & neigh = clu.neigh;

		for ( uint j=0; j < num_neighbors_this_cluster; j++ ) { // for each neighboring decoy...
			// If decoy is greater than names_.size(), just move on.
			// I don't yet know what could cause this.

			// output its decoy number and name...
			uint decoy = neigh[j];

			if ( decoy > names_.size() ) {
				std::cout << "Warning: attempted to access the " << decoy << " element of names_ (" << names_.size() << ")!" << std::endl;
				continue;
			}

			//std::cout << decoy << "'" << names_[decoy] << "'";
			// ...and its distance
			//if (AdjacentList::list_mode_ == LIST)
			// std::cout << (*(clu->dist))[j] << ", ";
			//else
			// std::cout << clu->getD(decoy) << ", ";
		}
		std::cout << std::endl << std::endl;
	}
}


/**
* Get the names and Strus of the (first num_elements) indices in list into
* Names and PDBs respectively.
*/
void
Clustering::getPDBs(std::vector< std::string> & Names, std::vector< StruOP > & PDBs, std::vector<int> const & list, int num_elements)
{
	for ( int i=0; i < num_elements; i++ ) {
		PDBs.push_back( PDBs_[ list[ i ] ] );
		Names.push_back( names_[ list[ i ] ] );
	}
}


void
Clustering::listAdjacentLists()
{
	std::cout << "Listing all adjacentLists" << std::endl;
	for ( int i=0; i < n_pdbs_; i++ ) {
		AdjacentList const & clu = *adjacent_lists_[i];
		std::cout << i << "," << clu.which_ << "," << names_[clu.which_]
			<< " " << clu.num_neighbors_ << ": ";
		int size2 = clu.num_neighbors_;
		for ( int j=0; j < size2; j++ ) { // for each neighboring decoy...
			// output its decoy number and name...
			std::cout << clu.neigh[j] << "'" << names_[clu.neigh[j]] << "'";
			// ...and its distance
			if ( AdjacentList::list_mode_ == LIST ) {
				std::cout << clu.dist[j] << ", ";
			} else {
				std::cout << clu.getD(clu.neigh[j]) << ", ";
			}
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}


//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// Codes for building the AdjacentLists
//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

/**
* A simplistic clustering to help us find the adjacent lists --- this is one
* of Calibur's strategy.
* In the output, each decoy either belongs to a cluster, or is a cluster
* center.
*
* For each decoy from 0... n_pdbs_,
*   If it is not within CLUS_RADIUS of any currently registered cluster center
*  Add it as a cluster center
*   Else add it to the cluster
*/
void
Clustering::auxClustering()
{
	double lower, upper, upper_scud;
	for ( int i=0; i < n_pdbs_; i++ ) { // for each decoy
#ifdef _SHOW_PERCENTAGE_COMPLETE_
		printf("Auxiliary clustering... %4.1f%%\r", 100.*i/n_pdbs_);
		fflush(stdout);
#endif
		bool clustered = false;
		int numc = mCluCen.size();
		for ( int c = 0; c < numc; c++ ) { // for each center
			int cen = mCluCen[c];

			refBound(i, cen, lower, upper);
			if ( lower > CLU_RADIUS ) continue;

			if ( pref_use_sig && estD(i, cen) > CLU_RADIUS ) continue;

			if ( upper <= CLU_RADIUS ) {
				mD2C[ i ] = upper;
				mCen[ i ] = cen;
				clustered = true;
				auxiliary_clusters_[cen]->push_back(i);
				break;
			}

			if ( pref_use_scud ) {
				upper_scud = eucD(i, cen);
				if ( upper_scud <= CLU_RADIUS ) {
					mD2C[ i ] = upper_scud;
					mCen[ i ] = cen;
					clustered = true;
					auxiliary_clusters_[cen]->push_back(i);
					break;
				}
			}

			double d = trueD(i, cen);
			adjacent_lists_[i]->add(cen, d, false);
			adjacent_lists_[cen]->add(i, d, false);
			if ( d <= CLU_RADIUS ) {
				mD2C[ i ] = d;
				mCen[ i ] = cen;
				clustered = true;
				auxiliary_clusters_[cen]->push_back(i);
				break;
			}
		}
		if ( !clustered ) {
			/**
			* Declare i as a cluster center.
			* Important: itself must be added to its cluster.
			*/
			auxiliary_clusters_[ i ] = new std::vector<int>(0);
			auxiliary_clusters_[ i ]->push_back(i);
			mCluCen.push_back(i);
			mD2C[ i ] = 0;
			mCen[ i ] = i;
		}
	}
}


/** too slow
bool
Clustering::find(int which, std::vector<int> *elements)
{
int len=elements->size();
for (int i=0; i<len; i++)
{
if ((*elements)[i] == which) return true;
}
return false;
}
*/

/**
* Subroutine for cluster().
* Finds the neighbors for each decoy i within threshold and add these
* neighbors in adjacent_lists_[i]
*/
void
Clustering::buildAdjacentLists()
{
	double lower, upper;
	//, upper_scud;
	int numc = mCluCen.size();
	//std::cout << "Number of decoys=" << n_pdbs_
	//  << ", number of clusters=" << numc << std::endl;

	//=================================================================
	// We assume that every element belongs to exactly one cluster,
	// so we loop through all the clusters and determine if the cluster
	// element should be added to the neighbor set of some decoy.
	//=================================================================
#ifdef _SHOW_PERCENTAGE_COMPLETE_
	printf("\r");
#endif
	for ( int c=0; c < numc; c++ ) { // for each cluster center
		double d, _d;
		int cen = mCluCen[c];
#ifdef _SHOW_PERCENTAGE_COMPLETE_
		printf("Finding decoys' neighbors... completed %4.1f%%\r", 100.*c/numc);
		fflush(stdout);
#endif
		for ( int i=0; i < n_pdbs_; i++ ) { // for each decoy
			std::vector<int>* elements = auxiliary_clusters_[cen];
			int size = elements->size();

			//=========================================================
			// Case 1
			// Under this condition, every element in the currrent
			// cluster is a neighbor of each other.
			// Do not compute distance.
			// Set dist=-1 to indicate that it's added at this stage.
			//=========================================================
			if ( 2*CLU_RADIUS <= THRESHOLD ) {
				if ( mCen[ i ] == cen ) { // if i is an element of the cluster
#ifdef _ADD_LITE_MODE_
					if (AdjacentList::list_mode_ == LITE) {
						adjacent_lists_[i]->add(size);
						continue;
					}
#endif
					for ( int n=0; n < size; n++ ) {
						int elem = (*elements)[n];
						//if (elem != i) // do not add self
						adjacent_lists_[i]->add(elem, -1., true);
					}
					continue;
				}
			}

			//=========================================================
			// Case 2 and Case 3
			// (Use signature estimation for pruning)
			// If this condition is fulfilled none of the element in
			// the cluster is a neighbor of the current decoy
			//=========================================================

			refBound(i, cen, lower, upper);
			if ( lower - CLU_RADIUS > THRESHOLD ) {
				continue;
			}
			if ( upper + CLU_RADIUS <= THRESHOLD ) {
#ifdef _ADD_LITE_MODE_
				if (AdjacentList::list_mode_ == LITE) {
					adjacent_lists_[i]->add(size);
					continue;
				}
#endif
				for ( int n=0; n < size; n++ ) {
					adjacent_lists_[i]->add((*elements)[n], -6, true);
				}
				continue;
			}

			if ( pref_use_scud && eucD(i,cen) + CLU_RADIUS <= THRESHOLD ) {
#ifdef _ADD_LITE_MODE_
				if (AdjacentList::list_mode_ == LITE) {
					adjacent_lists_[i]->add(size);
					continue;
				}
#endif
				for ( int n=0; n < size; n++ ) {
					adjacent_lists_[i]->add((*elements)[n], -6, true);
				}
				continue;
			}

			if ( pref_use_sig && estD(i, cen) - CLU_RADIUS > THRESHOLD ) {
				continue;
			}

			//=========================================================
			// Find exact distance to center
			//=========================================================
			d = trueD(i, cen);
			adjacent_lists_[i]->add(cen, d, false);
			adjacent_lists_[cen]->add(i, d, false);

			//=========================================================
			// If this condition is fulfilled none of the element in
			// the cluster is a neighbor of the current decoy
			//=========================================================
			if ( d - CLU_RADIUS > THRESHOLD ) continue;

			//=========================================================
			// If this condition is fulfilled every element in
			// the cluster is a neighbor of the current decoy
			// When r <= d/2, this is also implied by
			// d <= THRESHOLD/2. or
			// d <= r
			// Both of which are more strict than d + r <= THRESHOLD
			//=========================================================
			if ( d + CLU_RADIUS <= THRESHOLD ) {
#ifdef _ADD_LITE_MODE_
				if (AdjacentList::list_mode_ == LITE) {
					adjacent_lists_[i]->add(size);
					continue;
				}
#endif
				for ( int n=0; n < size; n++ ) {
					adjacent_lists_[i]->add((*elements)[n], -2, true);
				}
				continue;
			}

			//========================================================
			// Consider each cluster element individually
			//========================================================
			for ( int j=0; j < size; j++ ) {
				int e = (*elements)[j];

				if ( mD2C[e]+d <= THRESHOLD ) {
#ifdef _ADD_LITE_MODE_
					if (AdjacentList::list_mode_ == LITE) {
						adjacent_lists_[i]->add(1);
						continue;
					}
#endif
					adjacent_lists_[i]->add(e, -3, true);
				} else {
					if ( pref_use_sig && estD(i, e) > THRESHOLD ) {
						continue;
					}

					refBound(i, e, lower, upper);
					if ( upper <= THRESHOLD ) {
#ifdef _ADD_LITE_MODE_
						if (AdjacentList::list_mode_ == LITE) {
							 adjacent_lists_[i]->add(1);
							 continue;
						}
#endif
						adjacent_lists_[i]->add(e, -4, true);
						continue;
					}

					if ( pref_use_scud && eucD(i, e) <= THRESHOLD ) {
#ifdef _ADD_LITE_MODE_
						if (AdjacentList::list_mode_ == LITE) {
							 adjacent_lists_[i]->add(1);
							 continue;
						}
#endif
						adjacent_lists_[i]->add(e, -7, true);
						continue;
					}

					if ( lower > THRESHOLD ) {
						continue;
					}

					_d = trueD(i, e);
					adjacent_lists_[e]->add(i, _d, false);
					adjacent_lists_[i]->add(e, _d, false);
					if ( _d <= THRESHOLD ) {
#ifdef _ADD_LITE_MODE_
						if (AdjacentList::list_mode_ == LITE) {
							 adjacent_lists_[i]->add(1);
							 continue;
						}
#endif
						adjacent_lists_[i]->add(e, _d, true);
					}
				}
			}
		}
	}
}

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// Codes for finding RMSD
//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

double
Clustering::estD(int i, int j)
{
	std::vector< double > const & sig1 = PDBs_[i]->signature_;
	std::vector< double > const & sig2 = PDBs_[j]->signature_;
	double rev = 0;
	for ( int c=0; c < mLen; c++ ) {
		double v = sig1[c]-sig2[c];
		rev += v*v;
	}
	return sqrt(rev/mLen);
}

double
Clustering::estD(Stru const & a, Stru const & b)
{
	std::vector< double > const & sig1 = a.signature_;
	std::vector< double > const & sig2 = b.signature_;
	double rev=0;
	for ( int c=0; c<mLen; c++ ) {
		double v = sig1[c]-sig2[c];
		rev += v*v;
	}
	return sqrt(rev/mLen);
}

double
Clustering::trueD(int i, int j)
{
	std::cout.flush();
	if ( i == j ) return 0;
	double d = adjacent_lists_[i]->getD(j);
	if ( d < _OVER_RMSD_ && d >= 0 ) return d;
	std::vector<double> coor1 = PDBs_[i]->calpha_vector_;
	std::vector<double> coor2 = PDBs_[j]->calpha_vector_;
	double rmsd;
	rmsd = fast_rmsd(coor1, coor2, mLen);
	if ( rmsd != rmsd ) { // crazy RMSD
		rmsd = RMSD(coor1, coor2, mLen);
	}
	return (double) rmsd;
}

double
Clustering::eucD(int i, int j)
{
	std::cout.flush();
	if ( i == j ) return 0;
	double d = adjacent_lists_[i]->getD(j);
	if ( d < _OVER_RMSD_ && d >= 0 ) return d;
	std::vector<double> coor1 = PDBs_[i]->calpha_vector_;
	std::vector<double> coor2 = PDBs_[j]->calpha_vector_;
	double rev=0;
	int l3=mLen*3;
	for ( int k=0; k<l3; k++ ) {
		double d=coor1[k]-coor2[k];
		rev+=d*d;
	}
	return sqrt(rev/mLen);
}


/**
* Realign decoys so that their point-wise distance to decoy "ref" is minimized
*/
void
Clustering::realignDecoys(int ref)
{
	std::cout << "Realigning decoys...";
	clock_t start = clock();
	for ( int i=1; i < n_pdbs_; i++ ) {
		superimposeAndReplace(PDBs_[ref]->calpha_vector_, PDBs_[i]->calpha_vector_);
	}
	double elapsed = (clock() - start)/(double)CLOCKS_PER_SEC;
	std::cout << " completed in " << elapsed << " s" << std::endl;
}


void
Clustering::superimposeAndReplace(std::vector<double> & coor1, std::vector<double> & coor2)
{
	double R[3][3];
	RMSD(coor1, coor2, mLen, R);
	rotate(coor2, mLen, R, result_coords);
	for ( int i=0; i < 3*mLen; i++ ) {
		coor2[i] = result_coords[i];
	}
}


double
Clustering::trueD(Stru const & a, Stru const & b)
{
	std::vector<double> coor1 = a.calpha_vector_;
	std::vector<double> coor2 = b.calpha_vector_;
	double rmsd = 0;
	rmsd = fast_rmsd(coor1, coor2, mLen);
	if ( rmsd != rmsd ) { // crazy RMSD
		rmsd = RMSD(coor1, coor2, mLen);
	}
	return rmsd;
}


//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// Codes for recursively finding the largest clusters
//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

/**
* Find from all AdjacentLists the one that is longest
*/
int
Clustering::findDecoyWithMostNeighbors()
{
	int largest = remaining_list_[0];
	int currentSize = adjacent_lists_[largest]->num_neighbors_;
	for ( int i=1; i < mRemainingSize; i++ ) {
		int next_index = remaining_list_[i];
		int next_size = adjacent_lists_[next_index]->num_neighbors_;
		if ( next_size > currentSize ) {
			currentSize = next_size;
			largest = next_index;
		}
	}
	return largest;
}


/**
* Remove all the decoys in adj from all AdjacentList of all decoys.
*/
void
Clustering::removeDecoys(AdjacentList & adj)
{
	// remove from the remainingList
	std::vector<int> const & neigh = adj.neigh;
	int size = adj.num_neighbors_;
	int this_list = adj.which_;
	for ( int i=0; i<size; i++ ) {
		int to_remove = neigh[i]; // for each decoy to_remove

		// remove "to_remove" from the RemainingList
		int index = remaining_list_index_[to_remove];
		remaining_list_[index] = remaining_list_[mRemainingSize-1];
		remaining_list_index_[remaining_list_[index]] = index;
		mRemainingSize--;

		// remove "to_remove" from the AdjacentList of all its neighbors
		AdjacentListOP const & _neighbors_of_to_remove = adjacent_lists_[to_remove];
		std::vector<int> const & neighbors_of_to_remove = _neighbors_of_to_remove->neigh;
		int size2 = _neighbors_of_to_remove->num_neighbors_;
		for ( int j=0; j<size2; j++ ) {
			// should not remove to_remove from current list. it is used as the final cluster :D
			if ( neighbors_of_to_remove[j] != this_list ) {
				adjacent_lists_[neighbors_of_to_remove[j]]->remove(to_remove);
			}
		}
	}
}


void
Clustering::findLargestClusters()
{
	// already there mFinalClusters = std::shared_ptr< std::vector<AdjacentList *> >( new std::vector<AdjacentList *>(0) );
	mFinalClusters.clear();

	// We use a list of the remaining decoys (that has not been removed)
	// to enumerate only the decoys that are still available.
	remaining_list_.resize(n_pdbs_);
	remaining_list_index_.resize(n_pdbs_); // map decoy index to remaining_list_ index
	mRemainingSize = n_pdbs_;

	num_neighbors_.resize( n_pdbs_ );

	for ( int i=0; i<n_pdbs_; i++ ) {
		remaining_list_[i] = i;
		remaining_list_index_[i] = i;
		// Remember the original size of the cluster
		num_neighbors_[i] = adjacent_lists_[i]->num_neighbors_;
	}

#ifdef _ADD_LITE_MODE_
	if (AdjacentList::list_mode_ == LITE) {
		mFinalDecoy = findDecoyWithMostNeighbors();
		return;
	}
#endif

	while ( mRemainingSize > 1 ) {
		int largest = findDecoyWithMostNeighbors();
		AdjacentListOP adj(adjacent_lists_[largest]);
		mFinalClusters.push_back(adj); //std::shared_ptr< AdjacentList *>(adj));
		removeDecoys(*adj);
	}
}


/*
void
Clustering::findFinalRMSD()
{
int size = mFinalClusters->size();
for(int i=0; i < size; i++)
{
AdjacentList* adj=(*mFinalClusters)[i];
int w = adj->which_;
std::vector<int>* neigh = adj->neigh;
int num_neigh = adj->num_neighbors_;
for (int j=0; j < num_neigh; j++)
{
adj->add((*neigh)[j], trueD(w, (*neigh)[j]), false);
}
}
}
*/



//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// Random decoys finding codes
//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

/**
* Randomly permute the numbers, and select the first size decoys
*/
std::vector<std::string >
Clustering::getRandomDecoyNames(std::vector<std::string> const & srcnames, int size ) {
	int totalsize = srcnames.size();
	std::vector< int > randomArray( totalsize );
	for ( int i = 0; i < totalsize; i++ ) {// create an array of 0,...,totalsize-1
		randomArray[i] = i;
	}

	// swap the first size elements with other elements from the whole array
	for ( int i = 0; i < size; i++ ) {
		// find an index j (i<j<totalsize) to swap with the element i
		int j = i + numeric::random::rg().uniform() * (totalsize-i);
		int t = randomArray[j];
		randomArray[j] = randomArray[i];
		randomArray[i] = t;
	}
	// copy the first randomDecoysSize elements into a smaller array
	std::vector<std::string> names;
	for ( int i = 0; i < size; i++ ) {
		names.push_back(srcnames[randomArray[i]]);
	}
	return names;
}


/**
* Destroy the random decoys
*/
void
Clustering::destroyRandomDecoys(std::vector<std::string> const & /*names*/, std::vector< StruOP > const & /*decoys*/)
{
	/*
	delete names; // these names are shared

	// decoys should be read using readDecoys(), and should be safe to delete
	int size = decoys->size();
	for (int i = 0; i < size; i++)
	{
	delete (*decoys)[i];
	}
	delete decoys;
	*/
}


//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -
// Threshold finding codes
//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

/**
* For each decoy c in decoys, build a list of listLength nearest
* neighbors of c, taken from the list nborsCandidates.
* FIXME there is no way to tell if a decoy in sampleDecoys is the same
*    as a decoy in nborsCandidates, and hence a decoy can be inserted
*    into its own list.
*/
void
Clustering::get_neighbor_list(std::vector< StruOP > const & decoys,
	StringVec const & nborsCandidates,
	int listLength,
	std::vector< std::vector< double > > & nbors
) {
	// AMW TODO: at some point we should be using a proper 2Darray.
	// (we want to preserve the previous zero-indexing, so FArrays are out.)
	// For now, a correctly initialized vector of vectors will do.

	// For each decoy i find its neighbors: nbors[0],nbors[1],...,nbors[N-1]
	double r;
	std::string dName;
	int N = decoys.size();   // decoys to build neighbor list of
	int M = nborsCandidates.size(); // #candidates to test to fill list
	if ( listLength > M ) listLength = M;

	nbors.resize(N);
	for ( int i=0; i < N; i++ ) {
		nbors[i].resize( listLength, 1000.0 );
	}
	for ( int i=0; i < M; i++ ) { // for each candidate...
		dName = nborsCandidates[i];
		StruOP a( new Stru( SimPDBOP( new SimPDB( dName, mLen ) ), mLen, pref_use_sig ) );
		for ( int j=0; j < N; j++ ) {// ...insert it into each candidate list.
			std::vector< double > & nl = nbors[j];
			Stru const & b = *decoys[j];

			// FIXME need to change trueD to store computed RMSDs in a cache.
			//    this cache should be used to pre-fill the AdjacentLists
			r = trueD(*a, b);

			if ( r < nl[listLength-1] ) { // nbors[I] should be altered
				nl[listLength-1] = r;
				// double r up the list nl till it's larger
				for ( int k=listLength-1; k > 0 && nl[k-1] > r; k-- ) {
					// swap nl[k] with nk[k-1]
					nl[k] = nl[k-1];
					nl[k-1] = r;
				}
			}
		}
	}
}


/**
* Builds neighbor lists using the input decoys.
* For use in ROSETTA mode.
* Same as the earlier get_neighbor_list but nborsCandidates is taken to be
* the same as decoys.
* Also, minDist and maxDist are computed so that they can be used to
* guarantee targetClusterSize is used in getThreshold().
*/
void
Clustering::get_neighbor_list(std::vector< StruOP > const & decoys,
	int listLength,
	double & minDist,
	double & maxDist,
	std::vector< std::vector< double > > & nbors
) {
	// AMW: curious problems in the integration test using thres_finder 3 (ROSETTA)
	// may stem from issues in this function.
	// AMW: this seems fine. What about the one that removes?

	// For each decoy i find its neighbors: nbors[0],nbors[1],...,nbors[N-1]
	double r;
#ifdef _DEBUG_NBORLIST_
	int ** nborsIndex, * ni;
#endif
	//char * dName;
	int N = decoys.size();   // decoys to build neighbor list of
	if ( listLength > N ) listLength = N;
	maxDist = 0;
	minDist = _OVER_RMSD_;

	nbors.resize( N );
#ifdef _DEBUG_NBORLIST_
	nborsIndex = (int **) calloc(N, sizeof(int *));
#endif
	for ( int i=0; i < N; i++ ) {
		nbors[i].resize( listLength, 1000.0 );
#ifdef _DEBUG_NBORLIST_
		nborsIndex[i] = (int *) calloc(listLength, sizeof(double));
#endif
	}
	for ( int i=0; i < N; i++ ) {
		Stru const & a = *decoys[i];
		for ( int j=i+1; j < N; j++ ) {
			Stru const & b = *decoys[j];
			r = trueD(a, b);
			if ( r < minDist ) {
				minDist = r;
			}
			if ( r > maxDist ) {
				maxDist = r;
			}
			std::vector< double > & nl = nbors[i];
#ifdef _DEBUG_NBORLIST_
			ni = nborsIndex[i];
#endif
			if ( r < nl[listLength-1] ) { // nbors[I] should be altered
				nl[listLength-1] = r;
#ifdef _DEBUG_NBORLIST_
				ni[listLength-1] = j;
#endif
				// double r up the list nl till it's larger
				for ( int k=listLength-1; k > 0 && nl[k-1] > r; k-- ) {
					// swap nl[k] with nk[k-1]
					nl[k] = nl[k-1];
					nl[k-1] = r;
#ifdef _DEBUG_NBORLIST_
					ni[k] = ni[k-1];
					ni[k-1] = j;
#endif
				}
			}
			nl = nbors[j];
#ifdef _DEBUG_NBORLIST_
			ni = nborsIndex[j];
#endif
			if ( r < nl[listLength-1] ) { // nbors[I] should be altered
				nl[listLength-1] = r;
#ifdef _DEBUG_NBORLIST_
				ni[listLength-1] = i;
#endif
				// double r up the list nl till it's larger
				for ( int k=listLength-1; k > 0 && nl[k-1] > r; k-- ) {
					// swap nl[k] with nk[k-1]
					nl[k] = nl[k-1];
					nl[k-1] = r;
#ifdef _DEBUG_NBORLIST_
					ni[k] = ni[k-1];
					ni[k-1] = i;
#endif
				}
			}
		}
	}

#ifdef _DEBUG_NBORLIST_
	for ( uint i=0; i < names_.size(); i++ ) {
		std::cout << i << "'" << names_[i];
		ni = (int *)nborsIndex[i];
		for ( int j=0; j < listLength; j++ ) {
			int index = ni[j];
			std::cout << " " << index << "'" << names_[index];
		}
		std::cout << std::endl;
	}
#endif
}

static int
__cmp(const void *a, const void *b)
{
	return *((double *)a) >= *((double *)b);
}

/**
* Find:
* (1) min, max distances
* (2) most frequently occuring distance
* (3) the x percentile distance. using the x percentile
*  distance is like claiming that there are only x percent of edges
*  between decoys, and that these edges are the nearest ones.
* Note that if pairwise distances follows Gaussian distribution then case (2)
* should be the same as case (3) at x=50%.
*/
void
Clustering::estimateDist(std::vector< StruOP > const & decoys, double xPercent,
	double & minDist, double & maxDist,
	double & mostFreqDist, double & xPercentileDist
) {
	//char * dName;
	double r;
	int numdecoys = decoys.size();
	int numofpairs = numdecoys * (numdecoys-1) / 2;
	// Has to be a raw array for later qsort.
	double * alldists = new double[numofpairs];
	int k = 0;
	maxDist = 0;
	minDist = _OVER_RMSD_;

	for ( int i=0; i < numdecoys; i++ ) {
		Stru const & a = *decoys[i];
		for ( int j=i+1; j < numdecoys; j++ ) {
			Stru const & b = *decoys[j];
			r = trueD(a, b);
			alldists[k] = r;
			k++;
			if ( r < minDist && r >= 0 ) { // be paranoid about crazy RMSD
				minDist = r;
			}
			if ( r > maxDist && r < _OVER_RMSD_ ) { // crazy RMSD
				maxDist = r;
			}
		}
	}

	// Break the interval [minDist,maxDist] into numofbins bins and find
	// the frequencies of the distances within the range of each bin
	int numofbins = numdecoys;
	double bin_size = (maxDist - minDist)/numofbins;
	std::vector< int > bins( numofbins );
	for ( int & bin : bins ) {
		bin = 0;
	}
	for ( int i=0; i < numofpairs; i++ ) {
		int bin_to_put = (int)floor((alldists[i] - minDist) / bin_size);
		if ( bin_to_put >= numofbins || bin_to_put < 0 ) { // crazy RMSD
			continue;
		}
		bins[bin_to_put]++;
	}

	//std::cout << *minDist << "," << *maxDist << ", bin size=" << bin_size << std::endl;
	//for (int i=0; i < numofbins; i++)
	//{
	// std::cout << "bin " << i << "," << (*minDist + (i+0.5)*bin_size)
	//   << ": " << bins[i] << std::endl;
	//}

	// Find the peak of the distribution: Find n consecutive bins which
	// adds up to the highest frequency and return the center bin of these.

	// decide n and the length of one of n's side
	int n = numofbins / 10;
	if ( n % 2 == 0 ) { // use odd
		n++;
	}
	int n_half = (n-1)/2;

	int sum_of_n_bins = 0;
	for ( int i=0; i < n-1; i++ ) { // add up n-1 bins first
		sum_of_n_bins += bins[i];
	}

	int highest_frequency = 0;
	int bin_of_highest = 0;
	for ( int i = n_half; i < numofbins - n_half; i++ ) {
		sum_of_n_bins += bins[i + n_half]; // add the n-th bin
		if ( sum_of_n_bins > highest_frequency ) {
			highest_frequency = sum_of_n_bins;
			bin_of_highest = i;
		}
		sum_of_n_bins -= bins[i - n_half]; // remove the first bin
	}

	mostFreqDist = (minDist + (bin_of_highest + 0.5) * bin_size);

	// Find the x percentile distance through qsort
	qsort(alldists, numofpairs, sizeof(double), __cmp);
	xPercentileDist = alldists[(int)(xPercent/100*numofpairs)];

	if ( autoAdjustPercentile ) { // Use smaller thresholds for larger data sets

		xPercentile = 100./sqrt(sqrt(names_.size()));
		if ( xPercent > MAX_PERCENTILE_FOR_THRESHOLD ) {
			xPercent = MAX_PERCENTILE_FOR_THRESHOLD;
		}
		if ( xPercent < MIN_PERCENTILE_FOR_THRESHOLD ) {
			xPercent = MIN_PERCENTILE_FOR_THRESHOLD;
		}
		xPercentileDist = alldists[(int)(xPercent/100*numofpairs)];
	}

	delete [] alldists;
}


/**
* Find a decoy with the smallest distance x of N nearest neighbors, where
* minThres <= x <= maxThres.
*
* Hence at threshold x, no other cluster has more elements than N.
*
* We try for all sizes N in [minSize, maxSize] starting from size.
* The direction that N is searched depends on the starting size
* If there is no N within these set sizes such that minThres <= x <= maxThres,
* we use either minSize or maxSize, depending on which is nearer.
*
* In which case,
*
* For example, given the case
*     ----------------------------------maxThres
*    |       |
*     ----------------------------------minThres
* threshold   |   smallest x   ______|
*  ^    |   at N _______/   |^x_1
*  |    |  ________/     |
*  |    |____/        |
*    N = minSize     N = maxSize
*  ---> N
*         ...we would use x=x_1 and N=maxSize.
*/
double
Clustering::get_threshold(
	std::vector< std::vector< double > > const & nbors,
	int numLists,
	int listLength,
	int minSize,
	int size,
	double minThres,
	double maxThres
) {
	int t = size;   // current number of nearest neighbors to try
	int last_t = 0; // the last number of nearest neighbors found suitable
	int maxSize = listLength;

	// I think we should initialize best_rmsd because there is no
	// true guarantee that size is between minSize and maxSize.
	// (So loop might never execute.)
	double rmsd, best_rmsd = _OVER_RMSD_, threshold;

	/** Look for a number t where:
	* 1. minSize <= t <= maxSize
	* 2. minThres <= smallest C(t) among all decoys <= maxThres
	* where C(t) for a decoy is the max rmsd of the decoy's t nearest neighbors
	*
	* At every stage below, no cluster of best_rmsd radies has size above t.
	*/
	while ( t <= maxSize && t >= minSize ) {
		/** smallest C(t) among all decoys */
		best_rmsd = _OVER_RMSD_;
		for ( int i=0; i < numLists; i++ ) {
			rmsd = nbors[i][t-1]; // C(t) of decoy t
			if ( rmsd < best_rmsd ) {
				best_rmsd = rmsd;  // smallest C(t) so far
			}
		}

		if ( best_rmsd < minThres ) { // rmsd too small, need to increase target...
			if ( last_t == t+1 ) break; // ...but, we reduced it just now!
			last_t = t;
			t++;
		} else if ( best_rmsd > maxThres ) { // rmsd too large, need to decrease t
			if ( last_t == t-1 ) break; // ...but, we increased it just now!
			last_t = t;
			t--;
		} else {
			std::cout << "Cluster size " << t << " => threshold=" << best_rmsd
				<< ". [" << minThres << "," << maxThres << "]" << std::endl;
			break; // rmsd within threshold. Yahoo
		}
	}

	if ( t > maxSize ) t = maxSize;
	else if ( t < minSize ) t = minSize;

	if ( best_rmsd > maxThres || best_rmsd < minThres ) {
		std::cout << "No cluster size results in threshold within (" << minThres
			<< "," << maxThres << ")." << std::endl
			<< "Using threshold=" << best_rmsd
			<< ". Clusters' size will be at most " << t << std::endl;
	}

	threshold = best_rmsd;
	assert (t <= maxSize && t >= minSize);

	return threshold;
}


/**
* For each decoy j, and for each of the decoys i in index[], find and
* cache dist(i,j). index[] is of size REFERENCE_SIZE, and is initially
* set to decoys 1, 2, ..., REFERENCE_SIZE-1.
* This mNuPDBs_* distance is stored in references_
*/
void
Clustering::initRef( int* index ) {
	bool delete_later = false;
	if ( index == nullptr ) {
		delete_later = true;
		index = new int[REFERENCE_SIZE]; // a small memory leak here
		if ( pref_ref_use_first_few_decoys ) {
			for ( int i=0; i < REFERENCE_SIZE; i++ ) {
				index[i] = i;
			}
		} else {
			int step = n_pdbs_/REFERENCE_SIZE-1;
			int x = 0;
			for ( int i=0; i < REFERENCE_SIZE; i++ ) {
				index[i] = x;
				x += step;
			}
		}
	}

	//int numPDB = names_->size();
	references_.resize( REFERENCE_SIZE * n_pdbs_ );

	for ( int j=0; j < n_pdbs_; j++ ) {
		for ( int i=0; i < REFERENCE_SIZE; i++ ) {
			double d = trueD(index[i],j);
			adjacent_lists_[index[i]]->add(j, d, false);
			adjacent_lists_[j]->add(index[i], d, false);
			references_[j*REFERENCE_SIZE+i] = d; //trueD(index[i],j);
		}
	}

	if ( delete_later ) {
		delete [] index;
	}
}


void
Clustering::refBound( int i, int j, double& lower, double& upper )
{
	lower = 0;
	upper = _OVER_RMSD_;
	uint index_first  = i*REFERENCE_SIZE;
	uint index_second = j*REFERENCE_SIZE;
	for ( int k=0; k < REFERENCE_SIZE; k++ ) {
		double diff = fabs( references_[ index_first + k ] - references_[ index_second + k ] );
		double sum  = references_[ index_first + k ] + references_[ index_second + k ] ;
		if ( diff > lower ) lower = diff;
		if ( sum  < upper ) upper = sum;
	}
}

}
}
}

