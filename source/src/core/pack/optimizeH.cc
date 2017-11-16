// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/optimizeH.cc
/// @brief  standard hydrogen optimization subroutine
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/pack/optimizeH.hh>

// Package headers
#include <core/conformation/Residue.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Project headers
#include <basic/options/option.hh>
#include <basic/Tracer.hh> // RENAME THIS FILE WITH A CAPITAL T
#include <basic/basic.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <core/id/AtomID_Map.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {

/// Tracer instance for this file
/// Named after the original location of this code
static basic::Tracer TR( "core.io.pdb.file_data" );

void
optimize_H_and_notify(
	pose::Pose & pose,
	id::AtomID_Mask const & missing
)
{
	using namespace scoring;
	using namespace conformation;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ScoreFunctionOP sfxn;
	if ( !option[ score::optH_weights ].user() &&  !option[ score::optH_patch ].user() ) {
		sfxn = get_score_function();
	} else if ( option[ score::optH_weights ].user() && option[ score::optH_patch ].user() ) {
		sfxn = ScoreFunctionFactory::create_score_function( option[ score::optH_weights ](), option[ score::optH_patch ]() );
	} else if ( option[ score::optH_weights ].user() ) {
		sfxn = ScoreFunctionFactory::create_score_function( option[ score::optH_weights ]() );
	} else {
		// patch.user() and ! weights.user()  -- this is an odd request.
		sfxn = ScoreFunctionFactory::create_score_function( PRE_TALARIS_2013_STANDARD_WTS, option[ score::optH_patch ]() );
	}


	(*sfxn)(pose); // structure must be scored before pack_rotamers can be called.
	pose::Pose const start_pose( pose );
	pack::optimizeH( pose, *sfxn );
	// warn about positions that changed
	for ( Size i=1; i<= pose.size(); ++i  ) {
		Residue const & old_rsd( start_pose.residue(i) );
		Residue const & new_rsd(       pose.residue(i) );
		debug_assert( old_rsd.nchi() == new_rsd.nchi() && old_rsd.type().n_proton_chi() == new_rsd.type().n_proton_chi() );
		for ( Size chino=1; chino<= old_rsd.nchi(); ++chino ) {
			Real const chidev( std::abs( basic::subtract_degree_angles( old_rsd.chi( chino ), new_rsd.chi( chino ) ) ) );
			bool const chi_atom_was_missing( missing[ id::AtomID( old_rsd.chi_atoms( chino )[4], i ) ] );
			if ( chidev > 0.1 && !chi_atom_was_missing ) {
				using namespace ObjexxFCL::format;
				if ( old_rsd.type().is_proton_chi( chino ) ) {
					TR.Warning << "OPT-H: proton chi angle change: chidev= " << F(9,3,chidev) << " chino= " << I(2,chino) <<
						" position: " << I(4,i) << ' ' << old_rsd.name() << ' ' << new_rsd.name() << std::endl;
				} else {
					TR.Warning << "OPT-H: heavyatom chi angle change: chidev= " << F(9,3,chidev) << " chino= " <<
						I(2,chino) << " position: " << I(4,i) << ' ' << old_rsd.name() << ' ' << new_rsd.name() << std::endl;
				}
			}
		}
	}
}

void
optimizeH(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn
)
{
	using namespace task;
	using namespace basic::options;

	PackerTaskOP task = TaskFactory::create_packer_task( pose );

	task->initialize_from_command_line();
	task->or_optimize_h_mode( true );
	task->or_include_current( true );
	task->or_flip_HNQ( option[ OptionKeys::packing::flip_HNQ ] );
	task->or_multi_cool_annealer( option[ OptionKeys::packing::optH_MCA ]() );

	pack_rotamers( pose, sfxn, task );

}

}
}

