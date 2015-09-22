// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief protocols for folding into density
/// @details
/// @author Frank DiMaio

#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/electron_density/SetupForDensityScoringMoverCreator.hh>
#include <protocols/electron_density/util.hh>


#include <core/scoring/electron_density/util.hh>

#include <core/types.hh>
#include <utility/tag/Tag.hh>
#include <core/chemical/AA.hh>
#include <numeric/xyzVector.hh>

#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Residue.hh>

// Symmetry


#include <core/pose/util.hh>


#include <basic/options/option.hh>

// option key includes
#include <basic/options/keys/edensity.OptionKeys.gen.hh>


#include <basic/Tracer.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;


namespace protocols {
namespace electron_density {

static THREAD_LOCAL basic::Tracer TR( "protocols.electron_density.util" );

using namespace protocols;
using namespace core;

std::string
SetupForDensityScoringMoverCreator::keyname() const
{
	return SetupForDensityScoringMoverCreator::mover_name();
}

protocols::moves::MoverOP
SetupForDensityScoringMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetupForDensityScoringMover );
}

std::string
SetupForDensityScoringMoverCreator::mover_name()
{
	return "SetupForDensityScoring";
}

///////////////////////////////////////
///////////////////////////////////////

SetupForDensityScoringMover::SetupForDensityScoringMover() : Mover() {
	using namespace basic::options;
	dock_into_dens_strategy_ = option[ OptionKeys::edensity::realign ]();
}

void SetupForDensityScoringMover::apply( core::pose::Pose & pose ) {
	core::pose::addVirtualResAsRoot( pose );
	core::scoring::electron_density::getDensityMap().maskResidues( mask_reses_ );
	last_score = dockPoseIntoMap( pose, dock_into_dens_strategy_ );
	core::scoring::electron_density::getDensityMap().clearMask(  );
}

protocols::moves::MoverOP
SetupForDensityScoringMover::clone() const {
	return( protocols::moves::MoverOP( new SetupForDensityScoringMover( *this ) ) );
}

std::string
SetupForDensityScoringMover::get_name() const {
	return SetupForDensityScoringMoverCreator::mover_name();
}

void SetupForDensityScoringMover::mask( protocols::loops::Loops const & loops ) {
	mask_reses_.clear();
	for ( protocols::loops::Loops::LoopList::const_iterator it=loops.loops().begin(), it_end=loops.loops().end(); it != it_end; ++it ) {
		for ( core::Size i=it->start(), i_end=it->stop(); i<i_end; ++i ) {
			mask_reses_.push_back(i);
		}
	}
}

void SetupForDensityScoringMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const & ) {

	TR << "Parsing SetupForDensityScoringMover----" << std::endl;
	dock_into_dens_strategy_ = tag->getOption<std::string>( "realign", "no" );

	// TO DO: make mask parsable
}

///////////////////////////////////////
///////////////////////////////////////

void set_pose_and_scorefxn_for_edens_scoring( core::pose::Pose & pose, core::scoring::ScoreFunction &scorefxn ) {
	core::pose::addVirtualResAsRoot( pose );
	core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( scorefxn );
}


}
}

