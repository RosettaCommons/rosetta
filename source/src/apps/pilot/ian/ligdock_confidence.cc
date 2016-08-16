// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   demo/ian_test/ligdock_confidence.cc
///
/// @brief  Choose several ligand poses for further docking refinement.
/// @author Ian Davis (ian.w.davis@gmail.com)

// must be here to avoid VC++ ambiguous symbol w/ ObjexxFCL::byte
// for boinc builds - dek
#include <protocols/jobdist/JobDistributors.hh>


#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tools/make_vector1.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/option.hh>
#include <core/scoring/rms_util.hh>
#include <basic/Tracer.hh>

#include <protocols/jobdist/Jobs.hh>
#include <protocols/ligand_docking/LigandBaseProtocol.hh>


#include <fstream>
#include <set>
#include <sstream>


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>

#include <utility/excn/Exceptions.hh>


int
main( int argc, char * argv [] )
{
	try {

	using core::conformation::ResidueOP;
	using basic::options::option;
	using utility::vector1;
	using utility::file::FileName;
	using namespace basic::options::OptionKeys;
	using namespace protocols;
	using namespace protocols::jobdist;
	using namespace protocols::moves;
	basic::Tracer TR( "ligdock_confidence.main" );

	// Parses command line options and inits RNG.
	// Doesn't seem to hurt to do it again if already done once (?)
	devel::init(argc, argv);

	// A "native" pose for the diff reference point.
	// This used to be required but will be rarely used now.
	core::pose::PoseOP native_pose;
	if( option[ in::file::native ].user() ) {
		native_pose = new core::pose::Pose();
		core::import_pose::pose_from_file( *native_pose, option[ in::file::native ]().name() , core::import_pose::PDB_file);
	}

	vector1< FileName > atom_tree_diffs_file_names;
	if ( option[ in::file::silent ].user() )
		atom_tree_diffs_file_names = option[ in::file::silent ]().vector(); // make a copy
	else
		utility_exit_with_message("Must specify silent files via -in:file:silent");

	utility::io::ozstream out( "ligdock_confidence.tab" );
	out << "filename rms rmsnat lowclust_size lowclust_range lowclust_maxrms (num_overlap num_near_low+)+ dock_ok" << std::endl;
	for(core::Size file_idx = 1; file_idx <= atom_tree_diffs_file_names.size(); ++file_idx) {
		std::string const & atom_tree_diffs_file_name = atom_tree_diffs_file_names[file_idx];
		TR << "Reading atom_tree_diffs file " << atom_tree_diffs_file_name << " ... ";
		core::import_pose::atom_tree_diffs::AtomTreeDiff atdiff( atom_tree_diffs_file_name );
		core::import_pose::atom_tree_diffs::ScoresPairList const & scores_list = atdiff.scores();
		TR << scores_list.size() << " structures" << std::endl;

		// Keep only the top 5% by total score
		core::import_pose::atom_tree_diffs::ScoresPairList scores_list2;
		scores_list2.reserve( scores_list.size()/20 + 1 );
		protocols::ligand_docking::select_best_poses(atdiff, scores_list2);
		// Comes out already sorted by interface energy

		core::Real const min_rmsd = 1.0; //option[ docking::ligand::min_rms ];
		// Cluster heads -- sufficiently different from all other poses with lower interface_delta
		utility::vector1< core::pose::PoseOP > cluster_poses;
		utility::vector1< std::map< std::string, core::Real > > cluster_scores;
		utility::vector1< std::string > cluster_tags;
		// First N members of the first cluster (the one formed by the global minimum)
		core::Size const lowclust_N = 5;
		utility::vector1< core::pose::PoseOP > lowclust_poses;
		utility::vector1< std::map< std::string, core::Real > > lowclust_scores;
		// This is diagonal, to access rms_table[i][j], must have i > j
		utility::vector1< utility::vector1 < core::Real > > rms_table;
		// Loop over all selected poses, from lowest interface_delta on up
		for(core::Size i = 1; i <= scores_list2.size(); ++i) {
			std::string const & tag = scores_list2[i].first;
			core::pose::PoseOP a_pose = new core::pose::Pose();
			if( native_pose() == NULL ) atdiff.read_pose(tag, *a_pose);
			else atdiff.read_pose(tag, *a_pose, *native_pose);
			core::Size const last_rsd = a_pose->total_residue();
			// Compute rms to all existing cluster centers
			core::Real rms = 1e99;
			utility::vector1 < core::Real > rms_list;
			for(core::Size j = 1; j <= cluster_poses.size(); ++j) {
				core::Real this_rms = core::scoring::automorphic_rmsd(cluster_poses[j]->residue(last_rsd), a_pose->residue(last_rsd), false /*don't superimpose*/);
				rms_list.push_back(this_rms);
				rms = std::min(rms, this_rms);
				// Can't break early if we want to compute the full rms table (?)
				if(rms < min_rmsd) {
					TR << "Skip criteria: j = " << j << ", rms = " << rms << ", lowclust_poses.size() = " << lowclust_poses.size() << std::endl;
					if( j == 1 && lowclust_poses.size() < lowclust_N ) { // same cluster as the global minimum
						lowclust_poses.push_back( a_pose );
						lowclust_scores.push_back( scores_list2[i].second );
					}
					break;
				}
			}
			if( i == 1 ) {
				lowclust_poses.push_back( a_pose );
				lowclust_scores.push_back( scores_list2[i].second );
			}
			if(rms >= min_rmsd) {
				cluster_poses.push_back( a_pose );
				cluster_scores.push_back( scores_list2[i].second );
				cluster_tags.push_back( tag );
				rms_table.push_back( rms_list );
				TR << "Keeping  " << tag << " " << scores_list2[i].second["interface_delta"] << " " << scores_list2[i].second["ligand_auto_rms_no_super"] << std::endl;
			} else {
				TR << "Skipping " << tag << " " << scores_list2[i].second["interface_delta"] << " " << scores_list2[i].second["ligand_auto_rms_no_super"] << std::endl;
			}
		}

		// Calculate some statistics
		// We normalize all scores by dividing by the global minimum interface_delta
		core::Real const min_ifd = scores_list2[1].second["interface_delta"]; // same as lowclust_scores[1] and cluster_scores[1]
		core::Real const lowclust_max_ifd = lowclust_scores[ lowclust_scores.size() ]["interface_delta"];
		core::Real const lowclust_score_range = (1.0 - (lowclust_max_ifd / min_ifd)); // fractional range of scores in the lowest cluster
		out << atom_tree_diffs_file_name;
		out << " " << scores_list2[1].second["ligand_auto_rms_no_super"];
		out << " " << scores_list2[1].second["native_ligand_auto_rms_no_super"]; // defaults to zero if not present
		out << " " << lowclust_scores.size();
		out << " " << lowclust_score_range;
		core::Real max_rms = 0;
		for(core::Size i = 1; i <= lowclust_poses.size(); ++i) {
			for(core::Size j = i+1; j <= lowclust_poses.size(); ++j) {
				core::Size const last_rsd = lowclust_poses[j]->total_residue();
				core::Real const this_rms = core::scoring::automorphic_rmsd(lowclust_poses[j]->residue(last_rsd), lowclust_poses[i]->residue(last_rsd), false /*don't superimpose*/);
				max_rms = std::max(max_rms, this_rms);
			}
		}
		out << " " << max_rms; // max rms between any two poses in lowclust
		utility::vector1< core::Real > rms_cuts = utility::tools::make_vector1< core::Real >( 2.0, 4.0 );
		utility::vector1< core::Real > frac_cuts = utility::tools::make_vector1< core::Real >( 0.05, 0.10, 0.15, 0.20 );
		bool has_alt_minima = true; // one of my "by eye" criteria:  any clusters within 15% of the best and 2+A away?
		// How many other clusters score within __% of lowclust and have an rms of at least __ to it?
		for(core::Size i = 1; i <= rms_cuts.size(); ++i) {
			core::Real const rms_cut = rms_cuts[i];
			// Of those clusters, how many score better than the WORST member of lowclust?
			core::Size nbrs_below_lowclust = 0;
			for(core::Size k = 2; k <= cluster_scores.size(); ++k) {
				core::Real const cluster_rms = rms_table[k][1];
				if( cluster_rms > rms_cut && cluster_scores[k]["interface_delta"] < lowclust_max_ifd ) {
					nbrs_below_lowclust++;
					break;
				}
			}
			out << " " << nbrs_below_lowclust;
			// How many score within __% of the BEST member of lowclust (i.e. pose 1)
			for(core::Size j = 1; j <= frac_cuts.size(); ++j) {
				core::Real frac_cut = frac_cuts[j];
				core::Size other_clusters = 0;
				for(core::Size k = 2; k <= cluster_scores.size(); ++k) {
					core::Real const cluster_rms = rms_table[k][1];
					core::Real const frac_score = 1.0 - (cluster_scores[k]["interface_delta"] / min_ifd);
					if( cluster_rms > rms_cut && frac_score < frac_cut ) other_clusters++;
				}
				out << " " << other_clusters;
				if( rms_cut == 2.0 && frac_cut == 0.15 && other_clusters == 0 ) has_alt_minima = false;
			}
		}
		bool const well_converged = (lowclust_scores.size() >= 5 && lowclust_score_range < 0.10);
		out << " " << (well_converged || !has_alt_minima);
		out << std::endl;
	}
	out.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

