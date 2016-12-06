// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/kalngyk/calibur.cc
/// @brief  A tool for clustering decoys
/// @author YK Ng & SC Li (kalngyk@gmail.com)
/// @author Jared Adolf-Bryfogle (testing, cleanup, etc.)

#include <devel/init.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
//#include <basic/options/keys/cluster.OptionKeys.gen.hh>

#include <math.h>
#include <stdlib.h>

//Calibur Includes - Note this is EXTERNAL and is written outside of the coding conventions for SPEED.
// If you are intrepid, please convert these to Rosetta proper if possible!
#include <calibur/InitCluster.hh>
#include <calibur/SimpPDB.hh>

#include <utility/string_constants.hh>
#include <utility/string_util.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_1GRP_KEY( String, input, pdb_list ) // 1st default param
OPT_1GRP_KEY( Integer, res, start ) // r
OPT_1GRP_KEY( Integer, res, end   ) // r
OPT_1GRP_KEY( String,  res, chains ) // c
OPT_1GRP_KEY( Real,    strategy, thres ) // 2nd default param
OPT_1GRP_KEY( Boolean, strategy, nofilter ) // n
OPT_1GRP_KEY( Integer, strategy, thres_finder ) // t

int main(int argc, char** argv)
{
	try {
		NEW_OPT( input::pdb_list, "A file specifying the decoys. Each line of the file"
			" specifies a path (relative to the working directory) to a decoy's PDB file (-in:file:l also works here)",
			"" );
		NEW_OPT( res::start, "(Optional) starts from this Residues C-alpha atom (instead of starting"
			" from the first C-alpha atom)", 1);
		NEW_OPT( res::end, "(Optional) ends at this Residues C-alpha atom (instead of ending at"
			" the last C-alpha atom)", 0);
		NEW_OPT( res::chains, "(Optional) specifies the chains to be used."
			" By default, we scan for all alhpanumerics and empty chain, i.e. 'A', 'C', or unspecified", "");
		NEW_OPT( strategy::thres_finder, "(Optional) specifies the threshold finding strategy."
			" This should be one of 0, 1, 2, 3. (default strategy: 0) "
			" 0: threshold results in only x% of \"edges\" between decoys;"
			" default x=100/sqrt(sqrt(#decoys)) ||"
			" 1: threshold = min dist + x * (most frequent dist - min dist);"
			" default x=0.666667 (=2/3) ||"
			" 2: threshold = min dist + x * min(the avarage dist of decoys from a"
			" decoy; default x=0.666667 (=2/3) ||"
			" 3: find threshold using ROSETTA's method (auto-detect parameters).  May result in no clustering",
			0);
		NEW_OPT( strategy::nofilter, "(Optional) disables the filtering of outlier"
			" decoys.", false);
		NEW_OPT( strategy::thres, "(Optional) specifies a doubleing point number x."
			" which is used differently according to the threshold strategy specified."
			" x is ignored if the strategy does not use it."
			" x is used as the threshold if no threshold strategy is specified.",
			-1.0);

		devel::init(argc, argv);

		Clustering * ic = new Clustering();

		//JAB - thres_finder 4 does not work.  I'm disabling it instead of trying to fix it.

		// FIXME Make this user-definable
		SimPDB::preloadPDB = true;

		// Handle res::types
		std::string chains = utility::ALPHANUMERICS+" ";
		SimPDB::chains = strdup(chains.c_str());
		if ( option[res::chains].user() && strcmp(option[res::chains]().c_str(), "") ) {
			SimPDB::chains = strdup( option[res::chains]().c_str() );
		}

		// Handle res::{start,end}
		SimPDB::s_residue = 1;
		if ( option[res::start].user() ) {
			SimPDB::s_residue = option[res::start]();
			if ( SimPDB::s_residue < 1 ) {
				SimPDB::s_residue = 1;
			}
			if ( SimPDB::s_residue > LONGEST_CHAIN ) {
				SimPDB::s_residue = LONGEST_CHAIN;
			}
		}
		SimPDB::e_residue = LONGEST_CHAIN;
		if ( option[res::end].user() ) {
			SimPDB::e_residue = option[res::end]();
			if ( SimPDB::e_residue < SimPDB::s_residue ) {
				SimPDB::e_residue = SimPDB::s_residue;
			}
			if ( SimPDB::e_residue > LONGEST_CHAIN ) {
				SimPDB::e_residue = LONGEST_CHAIN;
			}
		}
		cout << "Using C-alphas #" << SimPDB::s_residue << "-";
		if ( SimPDB::e_residue == LONGEST_CHAIN ) {
			cout << "end" << endl;
		} else {
			cout << "#" << SimPDB::e_residue << endl;
		}

		// Handles strategy::nofilter
		Clustering::FILTER_MODE = (!option[strategy::nofilter]());

		// Handles strategy::thres_finder
		Clustering::EST_THRESHOLD = PERCENT_EDGES;
		bool strategy_specified = false;
		if ( option[strategy::thres_finder].user() ) {
			strategy_specified = true;
			switch (option[strategy::thres_finder]())
					{
					case 0 :
						Clustering::EST_THRESHOLD = PERCENT_EDGES;
						break;
					case 1 :
						Clustering::EST_THRESHOLD = MOST_FREQ_BASED;
						break;
					case 2 :
						Clustering::EST_THRESHOLD = MIN_AVG_DIST_BASED;
						break;
					case 3 :
						Clustering::EST_THRESHOLD = ROSETTA;
						break;
					case 4 :
						//Clustering::EST_THRESHOLD = SAMPLED_ROSETTA;
						utility_exit_with_message("Threshold strategy 4 no longer supported");
						break;
					default : break;
					}
		}

		char * filename;
		if ( option[input::pdb_list].user() ) {
			filename = strdup( option[input::pdb_list]().c_str() );
		} else if ( option[in::file::l].user() ) {
			filename = strdup( utility::to_string(option[in::file::l]()).c_str() );
		} else {
			utility_exit_with_message( "Missing -l or -pdb_list (please run with -help for usage)");
		}
		// Handles strategy::thres
		double threshold = -1;
		if ( option[strategy::thres].user() ) {
			double c = option[strategy::thres]();
			if ( !strategy_specified ) {
				threshold = c; // use supplied value as threshold
			} else { // use supplied value to guide the threshold finding strategy
				switch (Clustering::EST_THRESHOLD)
						{
						case MOST_FREQ_BASED:
						case MIN_AVG_DIST_BASED :
							Clustering::xFactor = c;
							break;
						case PERCENT_EDGES :
							Clustering::autoAdjustPercentile = false;
							Clustering::xPercentile = c;
							break;
						default :
							break; // ignore
						}
			}
		}

		/*
		cout << "filter mode=" << Clustering::FILTER_MODE << endl;
		cout << "strategy="    << Clustering::EST_THRESHOLD << endl;
		cout << "filename="    << filename << endl;
		cout << "xfactor="     << Clustering::xFactor << endl;
		cout << "xpercent="    << Clustering::xPercentile << endl;
		cout << "threshold="   << threshold << endl;
		cout << "residues="    << SimPDB::s_residue << "," << SimPDB::e_residue << endl;
		exit(0);
		*/

		ic->initialize(filename, threshold);
		ic->cluster();

		double acceptMargin = 0.15;
		if ( ic->bestClusMargin < acceptMargin ) {
			cout << "Best cluster larger than 2nd best cluster by only "
				<< (ic->bestClusMargin*100) << "% (<"
				<< (acceptMargin*100) << "%)" << endl
				<< "Two possible clusters could be present." << endl
				<< "Starting refined clustering..." << endl;

			// create new PDBs and Names out of the elements in the best two
			// clusters

			// first get the lists
			vector<AdjacentList *> * finalClusters = ic->mFinalClusters;
			vector<char *>* Names = new vector<char *>(0);
			vector<Stru *>* PDBs = new vector<Stru *>(0);

			// then add elements into them
			AdjacentList* clus;
			clus = (*finalClusters)[1];
			ic->getPDBs(Names, PDBs, clus->neigh, clus->mNumNeigh);
			clus = (*finalClusters)[0];
			ic->getPDBs(Names, PDBs, clus->neigh, clus->mNumNeigh);

			// Refined Clustering
			double minDist, maxDist, mostFreqDist, xPercentileDist;
			int numDecoys = Names->size() > 2*RANDOM_DECOY_SIZE_FOR_THRESHOLD?
				RANDOM_DECOY_SIZE_FOR_THRESHOLD: Names->size()/2;
			ic->estimateDist(Names,
				NUM_TRIALS_FOR_THRESHOLD,
				numDecoys,
				0.5,
				&minDist,
				&maxDist,
				&mostFreqDist,
				&xPercentileDist);
			ic->reinitialize(Names, PDBs, xPercentileDist);
			ic->cluster();

			if ( ic->bestClusMargin < acceptMargin ) {
				cout << "MORE THAN ONE BEST DECOYS DETECTED!" << endl;
			}
		}

		ic->showClusters(3);
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
