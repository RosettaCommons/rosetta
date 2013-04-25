// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file MinMover.cc
/// @brief
/// @author Ingemar Andre

// Unit headers
#include <protocols/symmetric_docking/SymFoldandDockSlideTrialMover.hh>

// Package headers
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymDockingInitialPerturbation.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>


// options
#include <basic/options/option.hh>
#include <basic/options/keys/fold_and_dock.OptionKeys.gen.hh>

// ObjexxFCL Headers

// C++ Headers

// Utility Headers
#include <basic/Tracer.hh>

#include <core/conformation/symmetry/SymDof.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace symmetric_docking {

static basic::Tracer TR("protocols.moves.symmetry.SymFoldandDockSlideTrialMover");

SymFoldandDockSlideTrialMover::SymFoldandDockSlideTrialMover()
	: Mover("SymFoldandDockSlideTrialMover") {}

SymFoldandDockSlideTrialMover::~SymFoldandDockSlideTrialMover(){}

void
SymFoldandDockSlideTrialMover::apply( core::pose::Pose & pose )
{
	protocols::simple_moves::symmetry::SetupForSymmetryMover setup;
	setup.apply( pose );

	using namespace core::conformation::symmetry;
	using namespace protocols::symmetric_docking;
	using namespace basic::options;

	assert( core::pose::symmetry::is_symmetric( pose ));
	SymmetricConformation & symm_conf (
        dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

	SymSlideInfo const & slide_info( symm_conf.Symmetry_Info()->get_slide_info() );

	TR.Debug << "Slide into contact mover..." << std::endl;
	if ( option[ OptionKeys::fold_and_dock::rotate_anchor_to_x ].user() ) {
		TR.Debug << "Rotate anchor to x axis.." << std::endl;
		core::pose::symmetry::rotate_anchor_to_x_axis( pose );
	}

	if ( slide_info.get_slide_type() == SEQUENTIAL ) {
		simple_moves::symmetry::SequentialSymmetrySlider symm_slider = simple_moves::symmetry::SequentialSymmetrySlider( pose,
																																			 slide_info.get_SlideCriteriaType(),
																																			 slide_info.get_SlideCriteriaVal() );
		symm_slider.apply( pose );
	}
	if ( slide_info.get_slide_type() == ORDERED_SEQUENTIAL ) {
		simple_moves::symmetry::OrderedSequentialSymmetrySlider symm_slider = simple_moves::symmetry::OrderedSequentialSymmetrySlider( pose,
																																									 slide_info.get_SlideCriteriaType(),
																																									 slide_info.get_SlideCriteriaVal(),
																																									 slide_info.get_slide_order() );
		symm_slider.apply( pose );
	}
	if ( slide_info.get_slide_type() == RANDOM ) {
		simple_moves::symmetry::RandomSymmetrySlider symm_slider = simple_moves::symmetry::RandomSymmetrySlider( pose,
																														 slide_info.get_SlideCriteriaType(),
																														 slide_info.get_SlideCriteriaVal() );
		symm_slider.apply( pose );
	}

}

std::string
SymFoldandDockSlideTrialMover::get_name() const {
	return "SymFoldandDockSlideTrialMover";
}

} // symmetric_docking
} // protocols
