// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MinMover.cc
/// @brief
/// @author ashworth

// Unit headers
#include <protocols/simple_moves/symmetry/SetupForSymmetryMoverCreator.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymDockingInitialPerturbation.hh>
#include <core/conformation/symmetry/SymmData.hh>

#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>

#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>


#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

// ObjexxFCL Headers

// C++ Headers
#include <string>

// Utility Headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves{
namespace symmetry{

static thread_local basic::Tracer TR( "protocols.simple_moves.symmetry.SetupForSymmetryMover" );

// creators
std::string
SetupForSymmetryMoverCreator::keyname() const {
	return SetupForSymmetryMoverCreator::mover_name();
}

protocols::moves::MoverOP
SetupForSymmetryMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetupForSymmetryMover );
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
	return protocols::moves::MoverOP( new ExtractAsymmetricUnitMover );
}

std::string
ExtractAsymmetricUnitMoverCreator::mover_name() {
	return "ExtractAsymmetricUnit";
}

////////////////////

std::string
ExtractAsymmetricPoseMoverCreator::keyname() const {
	return ExtractAsymmetricPoseMoverCreator::mover_name();
}

protocols::moves::MoverOP
ExtractAsymmetricPoseMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ExtractAsymmetricPoseMover );
}

std::string
ExtractAsymmetricPoseMoverCreator::mover_name() {
	return "ExtractAsymmetricPose";
}

////////////////////


SetupForSymmetryMover::SetupForSymmetryMover() :
	protocols::moves::Mover("SetupForSymmetryMover"),
	slide_(false),
	symmdef_()
{
}

SetupForSymmetryMover::SetupForSymmetryMover( core::conformation::symmetry::SymmDataOP symmdata ) :
	protocols::moves::Mover("SetupForSymmetryMover"),
	slide_(false),
	symmdef_( symmdata )
{}

SetupForSymmetryMover::SetupForSymmetryMover( std::string const & symmdef_file) :
	protocols::moves::Mover("SetupForSymmetryMover"),
	slide_(false),
	symmdef_()
{
	symmdef_ = core::conformation::symmetry::SymmDataOP( new core::conformation::symmetry::SymmData() );
	symmdef_->read_symmetry_data_from_file(symmdef_file);
}

SetupForSymmetryMover::~SetupForSymmetryMover(){}

void
SetupForSymmetryMover::apply( core::pose::Pose & pose )
{
	using namespace basic::options;

	// If we are alredy symmetric do nothing
	if ( core::pose::symmetry::is_symmetric( pose ) ) return;

	if(!symmdef_){
		if(option[ OptionKeys::symmetry::symmetry_definition].user()){
			symmdef_ = core::conformation::symmetry::SymmDataOP( new core::conformation::symmetry::SymmData() );
			symmdef_->read_symmetry_data_from_file(
				option[OptionKeys::symmetry::symmetry_definition]);

		} else {
			throw utility::excn::EXCN_BadInput(
				"The -symmetry:symmetry_definition command line option "
				"was not specified.");
		}
	}


	core::pose::symmetry::make_symmetric_pose( pose, *symmdef_ );
	assert( core::pose::symmetry::is_symmetric( pose ) );

	//fpd  explicitly update disulfide lr energy container
	if ( pose.is_fullatom() ) {
		core::scoring::symmetry::SymmetricScoreFunction disulf_score;
		disulf_score.set_weight( core::scoring::dslf_ss_dst, 1.0 );
		disulf_score.setup_for_scoring( pose );
	}

	// (Optionally) set rigid-body dofs from file
	//    SymDockingInitialPerturbation's behavior is controlled by flags and does nothing by default
	protocols::moves::MoverOP symdock( new protocols::simple_moves::symmetry::SymDockingInitialPerturbation(slide_) );
	symdock->apply( pose );
}

void SetupForSymmetryMover::parse_my_tag(
			utility::tag::TagCOP tag,
			basic::datacache::DataMap & /*data*/,
			filters::Filters_map const & /*filters*/,
			moves::Movers_map const & /*movers*/,
			core::pose::Pose const & /*pose*/ ) {

	using namespace basic::options;
	using namespace basic::resource_manager;

	if(tag->hasOption("definition") && tag->hasOption("resource_description")){
		throw utility::excn::EXCN_BadInput(
			"SetupForSymmetry takes either a 'definition' OR "
			"a 'resource_description' tag but not both.");
	}

	if(tag->hasOption("definition")){
		symmdef_ = core::conformation::symmetry::SymmDataOP( new core::conformation::symmetry::SymmData() );
		symmdef_->read_symmetry_data_from_file(
			tag->getOption<std::string>("definition"));
		option[OptionKeys::symmetry::symmetry_definition].value( "dummy" );

	} else if(tag->hasOption("resource_description")){
		symmdef_ = get_resource< core::conformation::symmetry::SymmData >(
			tag->getOption<std::string>("resource_description"));
		option[OptionKeys::symmetry::symmetry_definition].value( "dummy" );

	} else if(option[ OptionKeys::symmetry::symmetry_definition].user()){
		symmdef_ = core::conformation::symmetry::SymmDataOP( new core::conformation::symmetry::SymmData() );
		symmdef_->read_symmetry_data_from_file(
			option[OptionKeys::symmetry::symmetry_definition]);

	} else {
		throw utility::excn::EXCN_BadInput(
			"To use SetupForSymmetryMover with rosetta scripts please supply either a 'definition' tag, a 'resource_decription' tag or specify -symmetry:symmetry_definition the command line.");
	}
}

std::string
SetupForSymmetryMover::get_name() const {
	return SetupForSymmetryMoverCreator::mover_name();
}

////////////////////

ExtractAsymmetricUnitMover::ExtractAsymmetricUnitMover()
	: protocols::moves::Mover("ExtractAsymmetricUnitMover") { }

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
	utility::tag::TagCOP const /*tag*/,
	basic::datacache::DataMap & /*data*/,
	filters::Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/ ) { }

std::string
ExtractAsymmetricUnitMover::get_name() const {
	return ExtractAsymmetricUnitMoverCreator::mover_name();
}

/////////////////

ExtractAsymmetricPoseMover::ExtractAsymmetricPoseMover()
	: protocols::moves::Mover("ExtractAsymmetricPoseMover") { }

ExtractAsymmetricPoseMover::~ExtractAsymmetricPoseMover(){}


void
ExtractAsymmetricPoseMover::apply( core::pose::Pose & pose )
{
	//using namespace basic::options;
	// If we are not symmetric do nothing
	if ( !core::pose::symmetry::is_symmetric( pose ) ) return;

	core::pose::symmetry::make_asymmetric_pose( pose );
	//clear the symmetry deffinition option so it doesn't interfere with repack and minimization of
	//the assymetric pose.
	//option[OptionKeys::symmetry::symmetry_definition].clear();
}

void ExtractAsymmetricPoseMover::parse_my_tag(
	utility::tag::TagCOP const /*tag*/,
	basic::datacache::DataMap & /*data*/,
	filters::Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/ ) { }

std::string
ExtractAsymmetricPoseMover::get_name() const {
	return ExtractAsymmetricPoseMoverCreator::mover_name();
}


} //symmetry
} // simple_moves
} // protocols
