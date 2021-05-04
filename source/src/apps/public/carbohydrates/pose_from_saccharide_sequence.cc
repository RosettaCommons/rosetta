// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/public/carbohydrates/pose_from_saccharide_sequence.cc
/// @brief Convert a saccharide sequence into a minimized Pose structure
/// @author Morgan Nance (@mlnance)

// Unit Headers
#include <devel/init.hh>

// Project Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/analysis/GlycanInfoMover.hh>
#include <protocols/moves/Mover.hh>

#include <core/pose/carbohydrates/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/minimization_packing/MinMover.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/carbohydrates.OptionKeys.gen.hh>

#include <utility/tag/Tag.hh>

// Utility Headers
#include <cstdlib>

static basic::Tracer TR( "apps.public.carbohydrates.pose_from_saccharide_sequence" );

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init( argc, argv );

		// Create the pose object from a saccharide sequence
		TR << "Creating Pose from saccharide sequence" << std::endl;
		core::pose::PoseOP pose =
			core::pose::pose_from_saccharide_sequence
			( option[ carbohydrates::saccharide_sequence ].value() );
		// And assign all residues in the Pose to chain X
		pose->pdb_info()->set_chains('X');

		// Minimize the glycan conformation once
		core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function();
		core::select::residue_selector::TrueResidueSelectorCOP all_res_stor =
			utility::pointer::make_shared
			< core::select::residue_selector::TrueResidueSelector >();
		core::kinematics::MoveMapOP minimizer_mm =
			core::pose::carbohydrates::create_glycan_movemap_from_residue_selector
			( *pose, all_res_stor,
			true, // include_iupac_chi
			false, // include_glycan_ring_torsions
			true, // include_bb_torsions
			false ); // cartesian
		protocols::minimization_packing::MinMoverOP minimizer =
			utility::pointer::make_shared< protocols::minimization_packing::MinMover >
			( minimizer_mm, sf, "lbfgs_armijo_nonmonotone", 0.01, true );
		TR << "Minimizing" << std::endl;
		minimizer->apply( *pose );
		sf->show( TR, *pose );

		// Show the user some glycan info and output the structure
		protocols::analysis::GlycanInfoMover().apply( *pose );
		TR << "Writing Pose to 'out_saccharide.pdb'" << std::endl;
		pose->dump_pdb("out_saccharide.pdb");

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
