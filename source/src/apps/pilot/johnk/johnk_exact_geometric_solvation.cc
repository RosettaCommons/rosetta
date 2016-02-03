// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author jk

#include <iostream>
#include <iomanip>

#include <core/pose/Pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>

#include <basic/options/keys/score.OptionKeys.gen.hh>

//#include <core/scoring/geometric_solvation/ExactOccludedHbondSolEnergy.hh>
#include <core/scoring/geometric_solvation/non_scorefxn_exact_model.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_KEY( Boolean, report_LK )
OPT_KEY( Boolean, report_fitted )
OPT_KEY( Boolean, report_orig_geosol )
OPT_KEY( Boolean, report_exact )

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.johnk_exact_geometric_solvation.main" );

void dump_data( std::string const & fname, utility::vector1 <core::Real> & data ) {
  utility::io::ozstream outstream;
	outstream.open(fname, std::ios::out);
	for ( Size i = 1; i <= data.size(); i++ ) {
		//		outstream << curr_pose.pdb_info()->number( i ) << ' ' << residue_energies[i] << "\n";
		outstream << data[i] << "\n";
	}
	outstream.close();
	outstream.clear();
}

/// General testing code
int
main( int argc, char * argv [] )
{
    try {
	NEW_OPT( report_LK, "whether to report solvation energies from LK", false );
	NEW_OPT( report_fitted, "whether to report solvation energies from fitted occ", false );
	NEW_OPT( report_orig_geosol, "whether to old (Hbond energy directly) solvation energies", false );
	NEW_OPT( report_exact, "report the exact calculation", false );

	devel::init(argc, argv);

	TR << "Starting to compute various solvation energies" << std::endl;

	// scoring functions
	scoring::ScoreFunctionOP LK_scorefxn( get_score_function() );

	scoring::ScoreFunctionOP orig_geosol_scorefxn( get_score_function() );
	orig_geosol_scorefxn->set_weight( core::scoring::fa_sol, 0. );
	orig_geosol_scorefxn->set_weight( core::scoring::geom_sol, 0.65 );

	scoring::ScoreFunctionOP occ_pw_scorefxn( get_score_function() );
	occ_pw_scorefxn->set_weight( core::scoring::fa_sol, 0. );
	occ_pw_scorefxn->set_weight( core::scoring::occ_sol_fitted, 0.65 );

	scoring::ScoreFunctionOP occ_exact_scorefxn( get_score_function() );
	occ_exact_scorefxn->set_weight( core::scoring::fa_sol, 0. );
	occ_exact_scorefxn->set_weight( core::scoring::occ_sol_exact, 0.65 );

	std::string format_for_exact = "exact";
	if ( option[ score::exact_occ_pairwise ] ) format_for_exact += "_pw";
	if ( option[ score::exact_occ_split_between_res ] ) format_for_exact += "_splitres";

	for (core::Size f=1; f <= basic::options::start_files().size(); f++) {

		std::string const curr_decoy_fname = basic::options::start_files().at(f);
		TR << "Processing decoy " << curr_decoy_fname << std::endl;

		pose::Pose curr_pose;
		core::import_pose::pose_from_file( curr_pose, curr_decoy_fname , core::import_pose::PDB_file);

		// score with LK, turning off geosol because we change the weights later in this loop
		(*LK_scorefxn)(curr_pose);

		if ( option[report_LK] ) {
			// get LK solvation scores as a fxn of residue number
			utility::vector1 <core::Real> LK_energies;
			for ( Size i = 1; i <= curr_pose.total_residue(); i++ ) {
				LK_energies.push_back( curr_pose.energies().residue_total_energies(i)[ fa_sol ] );
			}
			dump_data( curr_decoy_fname+".LK.out", LK_energies );
		}

		if ( option[report_fitted] ) {
			// get pairwise geometric solvation scores as a fxn of residue number
			(*occ_pw_scorefxn)(curr_pose);
			utility::vector1 <core::Real> occ_pw_energies;
			for ( Size i = 1; i <= curr_pose.total_residue(); i++ ) {
				occ_pw_energies.push_back( curr_pose.energies().residue_total_energies(i)[ occ_sol_fitted ] );
			}
			dump_data( curr_decoy_fname+".fitted.out", occ_pw_energies );
		}

		if ( option[report_orig_geosol] ) {
			// get solvation scores as a fxn of residue number
			(*orig_geosol_scorefxn)(curr_pose);
			utility::vector1 <core::Real> orig_geosol_energies;
			for ( Size i = 1; i <= curr_pose.total_residue(); i++ ) {
				orig_geosol_energies.push_back( curr_pose.energies().residue_total_energies(i)[ geom_sol ] );
			}
			dump_data( curr_decoy_fname+".orig_geosol.out", orig_geosol_energies );
		}

		if ( option[report_exact] ) {
			// get exact ODO scores as a fxn of residue number
			(*occ_exact_scorefxn)(curr_pose);
			utility::vector1 <core::Real> occ_exact_energies;
			for ( Size i = 1; i <= curr_pose.total_residue(); i++ ) {
				occ_exact_energies.push_back( curr_pose.energies().residue_total_energies(i)[ occ_sol_exact ] );
			}
			dump_data( curr_decoy_fname+"."+format_for_exact+".out", occ_exact_energies );

			// compare to old exact code
			// jk note: this works, no need to keep doing this!
			//			utility::vector1 <core::Real> old_exact;
			//			core::scoring::geometric_solvation::compute_exact_geosol( curr_pose, false, option[ score::exact_occ_pairwise ],
			//				option[ score::exact_occ_split_between_res ], old_exact );
			//			dump_data( curr_decoy_fname+"."+format_for_exact+".OLD.out", old_exact );

		}

		TR << "moving on to next decoy" << std::endl;

	}

	TR << "Done computing solvation energies" << std::endl;

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}


