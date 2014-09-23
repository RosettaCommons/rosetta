// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/public/ligand_dock/select_best_unique_ligand_poses.cc
///
/// @brief  Choose several ligand poses for further docking refinement.
/// @author Ian Davis (ian.w.davis@gmail.com)
/// @author Rocco Moretti (rmoretti@u.washington.edu)

// must be here to avoid VC++ ambiguous symbol w/ ObjexxFCL::byte
// for boinc builds - dek
//#include <protocols/jobdist/JobDistributors.hh>

#include <devel/init.hh>

#include <protocols/ligand_docking/LigandBaseProtocol.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/AtomTreeDiffJobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>
#include <core/scoring/rms_util.hh>

#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

//Utility Includes

#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <utility/io/ozstream.hh>

// option key includes

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

// C++ Includes

#include <map>

#include <utility/excn/Exceptions.hh>

static thread_local basic::Tracer TR( "apps.public.ligand_dock.select_best_unique_ligand_poses" );

void
go_main() {
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using namespace protocols::jd2;

	/// INPUT

	JobInputterOP inputter( JobDistributorFactory::create_job_inputter() );
	protocols::jd2::Jobs inputs;
	inputter->fill_jobs( inputs );

	core::import_pose::atom_tree_diffs::ScoresPairList scores_list;
	std::map< std::string, JobOP> tag_job_map;
	for( Jobs::const_iterator job_iter(inputs.begin()); job_iter != inputs.end(); ++job_iter ) {
		std::string const & tag( (*job_iter)->input_tag() );
		tag_job_map[ tag ] = *job_iter;
		core::import_pose::atom_tree_diffs::Scores scoremap( (*job_iter)->output_string_real_pairs_begin(), (*job_iter)->output_string_real_pairs_end() );
		// Need to add data from those tucked into the pose itself (for binary silent file support).
		// We're actually only interested in "ligand_is_touching", "interface_delta",  and "total_score"
		if( scoremap.find("ligand_is_touching") == scoremap.end() || scoremap.find("interface_delta") == scoremap.end() || scoremap.find("total_score") == scoremap.end() ) {
			core::pose::PoseCOP in_pose = (*job_iter)->get_pose();
			runtime_assert( in_pose != 0 );
			core::Real value;
			if( scoremap.find("ligand_is_touching") == scoremap.end() && getPoseExtraScore( *in_pose, "ligand_is_touching", value ) ) {
				scoremap["ligand_is_touching"] = value;
			}
			if( scoremap.find("interface_delta") == scoremap.end() && getPoseExtraScore( *in_pose, "interface_delta", value ) ) {
				scoremap["interface_delta"] = value;
			}
			if( scoremap.find("total_score") == scoremap.end() ) {
				if( getPoseExtraScore( *in_pose, "total_score", value ) ) {
					scoremap["total_score"] = value;
				} else {
					core::scoring::EnergyMap const & emap( in_pose->energies().total_energies() );
					scoremap["total_score"] = emap[ core::scoring::total_score ];
				}
			}

		}

		core::import_pose::atom_tree_diffs::ScoresPair pair( tag, scoremap );
		scores_list.push_back(pair);
	}

	TR << scores_list.size() << " structures" << std::endl;

	/// PROCESSING
	// Keep only the top 5% by total score
	core::import_pose::atom_tree_diffs::ScoresPairList scores_list2;
	scores_list2.reserve( scores_list.size() );
	protocols::ligand_docking::select_best_poses(scores_list, scores_list2, option[ docking::ligand::subset_to_keep ] );
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
	for(core::Size i = 1; i <= scores_list2.size(); ++i) {
		if( selected_poses.size() >= max_poses ) break;
		std::string tag( scores_list2[i].first );
		a_pose = new core::pose::Pose();
		inputter->pose_from_job( *a_pose, tag_job_map[ tag ] );
		core::Size const last_rsd = a_pose->total_residue();
		core::Real rms = 1e99;
		utility::vector1 < core::Real > rms_list;
		for(core::Size j = 1; j <= selected_poses.size(); ++j) {
			// Can't break early if we want to compute the full rms table
			//if(rms < min_rmsd) break;
			core::Real this_rms = core::scoring::automorphic_rmsd(selected_poses[j]->residue(last_rsd), a_pose->residue(last_rsd), false /*don't superimpose*/);
			rms_list.push_back(this_rms);
			rms = std::min(rms, this_rms);
		}
		if(rms >= min_rmsd) {
			selected_poses.push_back( a_pose );
			selected_scores.push_back( scores_list2[i].second );
			selected_tags.push_back( tag );
			rms_table.push_back( rms_list );
			TR << "Keeping  " << tag << std::endl;
		} else {
			TR << "Skipping " << tag << std::endl;
		}
	}

	// OUTPUT
	// Write them out
	{
		AtomTreeDiffJobOutputterOP default_outputter( new AtomTreeDiffJobOutputter );
		default_outputter->set_precision(6, 3, 1);
		JobOutputterOP outputter( JobDistributorFactory::create_job_outputter( default_outputter ) );

		for(core::Size i = 1; i <= selected_poses.size(); ++i) {
			outputter->final_pose( tag_job_map[ selected_tags[i] ], *(selected_poses[i]) );
		}
	}

	// Write rmsd table
	{
		utility::io::ozstream out( "cluster_rms.tab" );
		// Column names
		for(core::Size i = 1; i <= selected_poses.size(); ++i) {
			out << ' ' << selected_tags[i];
		}
		out << '\n';
		// Rows
		for(core::Size i = 1; i <= selected_poses.size(); ++i) {
			out << selected_tags[i]; // row name
			for(core::Size j = 1; j <= selected_poses.size(); ++j) {
				if( i < j ) out << ' ' << rms_table[j][i];
				else if( i > j ) out << ' ' << rms_table[i][j];
				else out << " 0"; // i == j, self rms
			}
			out << '\n'; // end of row
		}
		out.close();
	}

	TR << "Selected " << selected_poses.size() << " poses overall" << std::endl;
}

int
main( int argc, char * argv [] ) {
	try {
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	// Parses command line options and inits RNG.
	// Doesn't seem to hurt to do it again if already done once (?)
	devel::init(argc, argv);

	if( option[ in::file::silent ].user() && ! option[ in::file::atom_tree_diff ].user() && ! option[ in::file::silent_struct_type ].user() ) {
		//Backwards compatability options munging
		std::string first_silent( *(option[ in::file::silent ]().begin()) );
		if( core::import_pose::atom_tree_diffs::file_is_atom_tree_diff( first_silent ) ) {
			TR.Warning << "WARNING: File " << first_silent << " looks to be AtomTypeDiff - reinterpreting -in:file:silent as -in:file:atom_tree_diff." << std::endl;
			TR.Warning << "         Explicitly set -in:file:silent_struct_type to override." << std::endl;
			option[ in::file::atom_tree_diff ]( option[ in::file::silent ] );
			option[ in::file::silent ].deactivate();
		}
	}
	if( option[ out::file::silent ].user() && ! option[ out::file::atom_tree_diff ].user() && ! option[ out::file::silent_struct_type ].user() ) {
		//Backwards compatability options munging
		TR.Warning << "WARNING: For backward compatibility, by default select_best_unique_ligand_poses will output Atom Tree Diff format files to -out:file:silent" << std::endl;
		TR.Warning << "         To output regular format silent file format, explicitly set -out:file:silent_struct_type" << std::endl;
		option[ out::file::atom_tree_diff ].value( option[ out::file::silent ] );
		option[ out::file::silent ].deactivate();
	}
	//Save users from themselves
	if( option[ in::file::keep_input_scores ] == false ) {
		TR.Warning << "WARNING: The program uses input scores, but -in:file:keep_input_scores is false. Resetting." << std::endl;
		option[in::file::keep_input_scores]( true );
	}
	if( option[ out::nstruct ] != 1 ) {
		TR.Warning << "WARNING: -out::nstruct other than one doesn't make sense. Resetting." << std::endl;
		option[out::nstruct]( 1 );
	}
	//Since nstruct is 1 ...
	if( ! option[ out::no_nstruct_label ] ) {
		TR << "Turning off nstruct labeling on output. Explicitly set -out:no_nstruct_label false to override." << std::endl;
		option[out::no_nstruct_label]( true );
	}

	go_main();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
