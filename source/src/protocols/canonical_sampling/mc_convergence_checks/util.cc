// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author


#include <protocols/canonical_sampling/mc_convergence_checks/Pool_ConvergenceCheck.hh>
#include <protocols/canonical_sampling/mc_convergence_checks/Heat_ConvergenceCheck.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/mc.OptionKeys.gen.hh>

#include <utility/sys_util.hh>
#include <utility/io/izstream.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace canonical_sampling {
namespace mc_convergence_checks {

using namespace core;
using namespace basic::options;
void setup_convergence_checks_from_cmdline( moves::MonteCarlo& mc ) {

	if ( option[ OptionKeys::mc::known_structures ].user() ) {
		/* this block is to take care of a
		particular problem with the archive... the file can be gone for short-whiles (during IO of Archive-process ) */
		Size trial( 100 );
		utility::io::izstream testin;
		do {
			utility::sys_sleep( 1 );
			testin.open( option[ OptionKeys::mc::known_structures ]() );
		} while ( !testin.good() && trial-- > 0 );
		// if ( trial == 0 ) throw CREATE_EXCEPTION(BadInput,  "can't open "+option[ OptionKeys::mc::known_structures ]() );
		/* okay file is now good, or will never be... */
		Pool_RMSD_OP pool_ptr( new Pool_RMSD(
			option[ OptionKeys::mc::known_structures ]()
			) );

		if ( option[ OptionKeys::mc::excluded_residues_from_rmsd ].user() ) {
			pool_ptr->set_excluded_residues( option[ OptionKeys::mc::excluded_residues_from_rmsd ]() );
		}

		mc.push_back( moves::MonteCarloExceptionConvergeOP( new Pool_ConvergenceCheck( pool_ptr, option[ OptionKeys::mc::max_rmsd_against_known_structures ]()  ) ));
		protocols::jd2::JobDistributor::get_instance()->job_outputter()->add_evaluation( evaluation::PoseEvaluatorOP( new Pool_Evaluator( pool_ptr ) ) );

	}

	/* other convergence checker */
	if ( option[ basic::options::OptionKeys::mc::heat_convergence_check ].user() ) {
		mc.push_back( moves::MonteCarloExceptionConvergeOP( new Heat_ConvergenceCheck() ) );
		mc.set_heat_after_cycles( option[ basic::options::OptionKeys::mc::heat_convergence_check ] );
	}

	/* add more checkers here ... */
}

}
}
}
