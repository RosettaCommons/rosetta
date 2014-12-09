// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/ian/cluster_ligand_poses.cc
///
/// @brief  Choose several ligand poses for further docking refinement.
/// @author Ian Davis (ian.w.davis@gmail.com)

// must be here to avoid VC++ ambiguous symbol w/ ObjexxFCL::byte
// for boinc builds - dek
#include <protocols/jobdist/JobDistributors.hh>


// AUTO-REMOVED #include <numeric/conversions.hh>
// AUTO-REMOVED #include <numeric/random/random_permutation.hh>
// AUTO-REMOVED #include <numeric/xyzVector.io.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1.io.hh>
#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>
#include <utility/pointer/owning_ptr.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <core/scoring/rms_util.hh>
#include <basic/Tracer.hh>

#include <protocols/cluster/APCluster.hh>
#include <protocols/jobdist/Jobs.hh>
// AUTO-REMOVED #include <protocols/jobdist/standard_mains.hh>
// AUTO-REMOVED #include <protocols/ligand_docking/LigandBaseProtocol.hh>


// AUTO-REMOVED #include <ctime>
#include <fstream>
#include <set>
#include <sstream>


// option key includes

#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>
#include <core/io/pdb/file_data.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <utility/excn/Exceptions.hh>




int
main( int argc, char * argv [] )
{
	try {

	using core::conformation::ResidueOP;
	using basic::options::option;
	using protocols::cluster::APCluster;
	using utility::vector1;
	using utility::file::FileName;
	using namespace basic::options::OptionKeys;
	using namespace protocols;
	using namespace protocols::jobdist;
	using namespace protocols::moves;
	basic::Tracer TR( "cluster_ligand_poses.main" );

	// Parses command line options and inits RNG.
	// Doesn't seem to hurt to do it again if already done once (?)
	devel::init(argc, argv);

	// A "native" pose for the diff reference point.
	// This used to be required but will be rarely used now.
	core::pose::PoseOP native_pose;
	if( option[ in::file::native ].user() ) {
		native_pose = new core::pose::Pose();
		core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]().name() );
	}

	vector1< FileName > atom_tree_diffs_file_names;
	if ( option[ in::file::silent ].user() )
		atom_tree_diffs_file_names = option[ in::file::silent ]().vector(); // make a copy
	else
		utility_exit_with_message("Must specify silent files via -in:file:silent");

	core::import_pose::atom_tree_diffs::ScoresPairList scores_list;
	vector1< ResidueOP > lig_residues;
	for(core::Size i = 1, i_end = atom_tree_diffs_file_names.size(); i <= i_end; ++i) {
		TR << "Reading atom_tree_diffs file " << atom_tree_diffs_file_names[i] << " ... ";
		core::import_pose::atom_tree_diffs::AtomTreeDiff atdiff( atom_tree_diffs_file_names[i] );
		core::import_pose::atom_tree_diffs::ScoresPairList const & new_scores_list = atdiff.scores();
		TR << new_scores_list.size() << " structures" << std::endl;
		//// Keep only the top 5% by total score
		//AtomTreeDiff::ScoresPairList scores_list2;
		//scores_list2.reserve( scores_list.size() );
		//protocols::ligand_docking::select_best_poses(atdiff, scores_list2);
		//// Comes out already sorted by interface energy
		for(core::Size j = 1, j_end = new_scores_list.size(); j <= j_end; ++j) {
			scores_list.push_back( new_scores_list[j] );
			core::pose::Pose a_pose;
			std::string const & tag = new_scores_list[j].first;
			if( native_pose() == NULL ) atdiff.read_pose(tag, a_pose);
			else atdiff.read_pose(tag, a_pose, *native_pose);
			// Keep a copy of the last residue in the pose and discard the rest to conserve memory
			lig_residues.push_back( a_pose.residue(a_pose.total_residue()).clone() );
			// Have to use std::cout because Tracers produce too much junk here:
			if( j % 10 == 0 ) std::cout << "." << std::flush;
		}
		std::cout << std::endl;
	}

	core::Size const max_poses = option[ docking::ligand::max_poses ];
	//core::Real const min_rmsd = option[ docking::ligand::min_rms ];
	//core::pose::PoseOP a_pose;
	//utility::vector1< core::pose::PoseOP > selected_poses;
	//utility::vector1< std::map< std::string, core::Real > > selected_scores;
	//utility::vector1< std::string > selected_tags;
	//// This is diagonal, requires funny indexing to access
	//utility::vector1< utility::vector1 < core::Real > > rms_table;

	// Calculate RMSDs between all ligand residues and store in a clustering object
	TR << "Calculating RMSDs for " << scores_list.size() << " total structures: " << std::endl;
	APCluster cluster(lig_residues.size(), max_poses);
	core::Real max_rms = 0, sum_rms = 0, min_ifd = -0.1;
	for(core::Size i = 1, i_end = lig_residues.size(); i <= i_end; ++i) {
		for(core::Size j = i+1, j_end = lig_residues.size(); j <= j_end; ++j) {
			core::Real const this_rms = core::scoring::automorphic_rmsd(*lig_residues[i], *lig_residues[j], false /*don't superimpose*/);
			cluster.set_sim(i, j, -this_rms);
			cluster.set_sim(j, i, -this_rms);
			sum_rms += this_rms;
			if( this_rms > max_rms ) max_rms = this_rms;
		}
		core::Real const this_ifd = scores_list[i].second["interface_delta"];
		if( this_ifd < min_ifd ) min_ifd = this_ifd;
			// Have to use std::cout because Tracers produce too much junk here:
		if(i%10 == 0) std::cout << "." << std::flush;
	}
	std::cout << std::endl;

	// Assign preferences (self-similarities) for each point
	// I'd rather know the median, but that can't be done in constant space
	// The purpose of this is to indicate that some poses (those with good energies)
	// are a priori better choices as cluster centers (exemplars).
	// Using the median (here, mean) similarity for all points would give a moderate
	// number of clusters, and using the minimum (max rmsd) would give few clusters.
	// Here, we interpolate linearly between the mean (for the pose with the absolute lowest interface delta)
	// and the minimum similarity (for poses with interface delta >= 0).
	// This hasn't been tested much and better choices are probably possible.
	// Other possibilities include ranking the poses by some metric (e.g. interface delta)
	// and using their ranks to compute their preferences.
	core::Real const mean_rms = sum_rms / ( lig_residues.size() * (lig_residues.size()-1) * 0.5 );
	TR << "Mean RMSD = " << mean_rms << ", max RMSD = " << max_rms << std::endl;
	for(core::Size i = 1, i_end = lig_residues.size(); i <= i_end; ++i) {
		core::Real const this_ifd = scores_list[i].second["interface_delta"];
		// a = 0 is best ifd; a = 1 is worst ifd
		core::Real a = (this_ifd - min_ifd) / (0.0 - min_ifd);
		if( a > 1 ) a = 1; // for things with a positive value of interface_delta
		core::Real const pref = a*max_rms + (1-a)*mean_rms;
		//TR << "ifd = " << this_ifd << ", a = " << a << ", pref = " << pref << std::endl;
		cluster.set_sim(i, i, -pref);
	}

	TR << "Clustering poses..." << std::endl;
	cluster.cluster(200, 10, 0.9); // these could probably be tuned some more...
	TR << cluster.get_num_exemplars() << " clusters selected" << std::endl;
	cluster.save_binary("ligand_clusters.bin");

	utility::io::ozstream out("ligand_clusters.tab");
	out << "#size exemplar low_interface_delta low_total_score ,member1,member2,...\n";
	vector1< core::Size > exemplars, members;
	cluster.get_all_exemplars( exemplars );
	for(core::Size i = 1, i_end = exemplars.size(); i <= i_end; ++i) {
		core::Size const ii = exemplars[i];
		//TR << "Exemplar: " << scores_list[ii].first << "; members";
		cluster.get_cluster_for(ii, members);
		core::Size low_ifd_idx = ii, low_tot_idx = ii;
		core::Real low_ifd_val = scores_list[ii].second["interface_delta"], low_tot_val =scores_list[ii].second["total_score"];
		for(core::Size j = 1, j_end = members.size(); j <= j_end; ++j) {
			core::Size const jj = members[j];
			core::Real const this_ifd = scores_list[jj].second["interface_delta"];
			core::Real const this_tot = scores_list[ii].second["total_score"];
			if( this_ifd < low_ifd_val ) {
				low_ifd_idx = jj;
				low_ifd_val = this_ifd;
			}
			if( this_tot < low_tot_val ) {
				low_tot_idx = jj;
				low_tot_val = this_tot;
			}
		}
		out << members.size() << " " << scores_list[ii].first;
		out << " " << low_ifd_val << " " << scores_list[low_ifd_idx].first;
		out << " " << low_tot_val << " " << scores_list[low_tot_idx].first;
		out << " ";
		for(core::Size j = 1, j_end = members.size(); j <= j_end; ++j) {
			core::Size const jj = members[j];
			out << "," << scores_list[jj].first;
			//TR << " " << scores_list[jj].first;
		}
		out << "\n";
		//TR << std::endl;
	}
	out.close();

	} catch ( utility::excn::EXCN_Base const & e ) { //YOU ADD
		std::cout << "caught exception " << e.msg() << std::endl; //YOU ADD
		return -1;
	}

	return 0;
}

