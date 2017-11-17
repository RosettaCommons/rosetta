// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   demo/ian_test/select_best_unique_ligand_poses.cc
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

#include <devel/init.hh>
#include <core/types.hh>
#include <basic/options/option.hh>
#include <core/scoring/rms_util.hh>
#include <basic/Tracer.hh>

#include <protocols/jobdist/Jobs.hh>
#include <protocols/ligand_docking/LigandBaseProtocol.hh>


#include <fstream>
#include <set>
#include <sstream>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>


int
main( int argc, char * argv [] )
{
	try {

		using basic::options::option;
		using namespace basic::options::OptionKeys;
		using namespace protocols;
		using namespace protocols::jobdist;
		using namespace protocols::moves;
		basic::Tracer TR( "select_best_unique_ligand_poses.main" );

		// Parses command line options and inits RNG.
		// Doesn't seem to hurt to do it again if already done once (?)
		devel::init(argc, argv);

		// A "native" pose for the diff reference point.
		// This used to be required but will be rarely used now.
		core::pose::PoseOP native_pose;
		if ( option[ in::file::native ].user() ) {
			native_pose = core::pose::PoseOP( new core::pose::Pose() );
			core::import_pose::pose_from_file( *native_pose, option[ in::file::native ]().name() , core::import_pose::PDB_file);
		}

		TR << "Reading silent file " << std::string( *(option[ in::file::silent ]().begin()) ) << " ... ";
		core::import_pose::atom_tree_diffs::AtomTreeDiff atdiff( *(option[ in::file::silent ]().begin()) );
		core::import_pose::atom_tree_diffs::ScoresPairList const & scores_list = atdiff.scores();
		TR << scores_list.size() << " structures" << std::endl;

		// Keep only the top 5% by total score
		core::import_pose::atom_tree_diffs::ScoresPairList scores_list2;
		scores_list2.reserve( scores_list.size() );
		protocols::ligand_docking::select_best_poses(atdiff, scores_list2);
		// Comes out already sorted by interface energy

		core::Size const max_poses = option[ docking::ligand::max_poses ];
		core::Real const min_rmsd = option[ docking::ligand::min_rms ];
		core::pose::PoseOP a_pose;
		utility::vector1< core::pose::PoseOP > selected_poses;
		utility::vector1< std::map< std::string, core::Real > > selected_scores;
		utility::vector1< std::string > selected_tags;
		// This is diagonal, requires funny indexing to access
		utility::vector1< utility::vector1 < core::Real > > rms_table;

		// Cycle through ranked poses, keeping them if they're at least min_rmsd
		// away from all other poses being kept so far.
		for ( core::Size i = 1; i <= scores_list2.size(); ++i ) {
			if ( selected_poses.size() >= max_poses ) break;
			std::string tag( scores_list2[i].first );
			a_pose = core::pose::PoseOP( new core::pose::Pose() );
			if ( native_pose == NULL ) atdiff.read_pose(tag, *a_pose);
			else atdiff.read_pose(tag, *a_pose, *native_pose);
			core::Size const last_rsd = a_pose->size();
			core::Real rms = 1e99;
			utility::vector1 < core::Real > rms_list;
			for ( core::Size j = 1; j <= selected_poses.size(); ++j ) {
				// Can't break early if we want to compute the full rms table
				//if(rms < min_rmsd) break;
				core::Real this_rms = core::scoring::automorphic_rmsd(selected_poses[j]->residue(last_rsd), a_pose->residue(last_rsd), false /*don't superimpose*/);
				rms_list.push_back(this_rms);
				rms = std::min(rms, this_rms);
			}
			if ( rms >= min_rmsd ) {
				selected_poses.push_back( a_pose );
				selected_scores.push_back( scores_list2[i].second );
				selected_tags.push_back( tag );
				rms_table.push_back( rms_list );
				TR << "Keeping  " << tag << std::endl;
			} else {
				TR << "Skipping " << tag << std::endl;
			}
		}

		// Write them out
		{
			utility::io::ozstream out( option[out::file::silent]() );
			std::map< std::string, core::Real > const empty_scores;
			core::pose::PoseCOP ref_pose = selected_poses[1];
			core::import_pose::atom_tree_diffs::dump_reference_pose(out, "REF", empty_scores, *ref_pose);
			for ( core::Size i = 1; i <= selected_poses.size(); ++i ) {
				core::import_pose::atom_tree_diffs::dump_atom_tree_diff(out, selected_tags[i], selected_scores[i], *ref_pose, *(selected_poses[i]), 6, 3, 1);
			}
			out.close();
		}

		// Write rmsd table
		{
			utility::io::ozstream out( "cluster_rms.tab" );
			// Column names
			for ( core::Size i = 1; i <= selected_poses.size(); ++i ) {
				out << ' ' << selected_tags[i];
			}
			out << '\n';
			// Rows
			for ( core::Size i = 1; i <= selected_poses.size(); ++i ) {
				out << selected_tags[i]; // row name
				for ( core::Size j = 1; j <= selected_poses.size(); ++j ) {
					if ( i < j ) out << ' ' << rms_table[j][i];
					else if ( i > j ) out << ' ' << rms_table[i][j];
					else out << " 0"; // i == j, self rms
				}
				out << '\n'; // end of row
			}
			out.close();
		}

		TR << "Selected " << selected_poses.size() << " poses overall" << std::endl;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

