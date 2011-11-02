// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file MinMover.cc
/// @brief
/// @author ashworth

// Unit headers
#include <protocols/moves/symmetry/SetupForSymmetryMoverCreator.hh>
#include <protocols/moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/symmetric_docking/SymDockingInitialPerturbation.hh>

// AUTO-REMOVED #include <protocols/moves/DataMap.hh>
// AUTO-REMOVED #include <protocols/rosetta_scripts/util.hh>
#include <utility/tag/Tag.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

#include <core/scoring/symmetry/SymmetricScoreFunction.hh>


// ObjexxFCL Headers

// C++ Headers

// Utility Headers
#include <basic/Tracer.hh>

#include <protocols/jobdist/Jobs.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace moves {
namespace symmetry {

static basic::Tracer TR("protocols.moves.symmetry.SetupForSymmetryMover");

// creators
std::string
SetupForSymmetryMoverCreator::keyname() const {
	return SetupForSymmetryMoverCreator::mover_name();
}

protocols::moves::MoverOP
SetupForSymmetryMoverCreator::create_mover() const {
	return new SetupForSymmetryMover;
}

std::string
SetupForSymmetryMoverCreator::mover_name() {
	return "SetupForSymmetry";
}

std::string
ExtractAsymmetricUnitMoverCreator::keyname() const {
	return ExtractAsymmetricUnitMoverCreator::mover_name();
}

protocols::moves::MoverOP
ExtractAsymmetricUnitMoverCreator::create_mover() const {
	return new ExtractAsymmetricUnitMover;
}

std::string
ExtractAsymmetricUnitMoverCreator::mover_name() {
	return "ExtractAsymmetricUnit";
}

////////////////////

SetupForSymmetryMover::SetupForSymmetryMover()
	: Mover("SetupForSymmetryMover"), slide_(false), symmdef_file_("") { }

SetupForSymmetryMover::SetupForSymmetryMover( std::string const & symmdef_file)
	: Mover("SetupForSymmetryMover"), slide_(false), symmdef_file_(symmdef_file) { }


SetupForSymmetryMover::~SetupForSymmetryMover(){}

void
SetupForSymmetryMover::apply( core::pose::Pose & pose )
{
	// If we are alredy symmetric do nothing
	if ( core::pose::symmetry::is_symmetric( pose ) ) return;

	core::pose::symmetry::make_symmetric_pose( pose, symmdef_file_ );
	assert( core::pose::symmetry::is_symmetric( pose ) );

	//fpd  explicitly update disulfide lr energy container
	if ( pose.is_fullatom() ) {
		core::scoring::symmetry::SymmetricScoreFunction disulf_score;
		disulf_score.set_weight( core::scoring::dslf_ss_dst, 1.0 );
		disulf_score.setup_for_scoring( pose );
	}

	// (Optionally) set rigid-body dofs from file
	//    SymDockingInitialPerturbation's behavior is controlled by flags and does nothing by default
	protocols::moves::MoverOP symdock =
		new protocols::symmetric_docking::SymDockingInitialPerturbation(slide_);
	symdock->apply( pose );
}

void SetupForSymmetryMover::parse_my_tag( 
			utility::tag::TagPtr const tag,
			moves::DataMap & /*data*/,
			filters::Filters_map const & /*filters*/,
			moves::Movers_map const & /*movers*/,
			core::pose::Pose const & /*pose*/ ) { 
	symmdef_file_ = tag->getOption<std::string>("definition", "");
}

std::string
SetupForSymmetryMover::get_name() const {
	return SetupForSymmetryMoverCreator::mover_name();
}

////////////////////

ExtractAsymmetricUnitMover::ExtractAsymmetricUnitMover()
	: Mover("ExtractAsymmetricUnitMover") { }

ExtractAsymmetricUnitMover::~ExtractAsymmetricUnitMover(){}

void
ExtractAsymmetricUnitMover::apply( core::pose::Pose & pose )
{
	// If we are not symmetric do nothing
	if ( !core::pose::symmetry::is_symmetric( pose ) ) return;

	core::pose::Pose pose_asu;
	core::pose::symmetry::extract_asymmetric_unit(pose, pose_asu);
	pose = pose_asu;
}

void ExtractAsymmetricUnitMover::parse_my_tag( 
			utility::tag::TagPtr const tag,
			moves::DataMap & /*data*/,
			filters::Filters_map const & /*filters*/,
			moves::Movers_map const & /*movers*/,
			core::pose::Pose const & /*pose*/ ) { }

std::string
ExtractAsymmetricUnitMover::get_name() const {
	return ExtractAsymmetricUnitMoverCreator::mover_name();
}


} // symmetry
} // moves
} // protocols
