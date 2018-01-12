// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SaneMinMover.cc
/// @brief
/// @author James Thompson

#include <protocols/minimization_packing/SaneMinMover.hh>
#include <protocols/minimization_packing/SaneMinMoverCreator.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/prof.hh>

#include <utility>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace minimization_packing {

static basic::Tracer TR( "protocols.minimization_packing.SaneMinMover" );

// XRW TEMP std::string
// XRW TEMP SaneMinMoverCreator::keyname() const {
// XRW TEMP  return SaneMinMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SaneMinMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SaneMinMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SaneMinMover::mover_name() {
// XRW TEMP  return "SaneMinMover";
// XRW TEMP }

// default constructor
SaneMinMover::SaneMinMover() : protocols::moves::Mover("SaneMinMover") {
	set_defaults_();
}

SaneMinMover::SaneMinMover( std::string const & name ) :
	protocols::moves::Mover(name) {
	set_defaults_();
}

SaneMinMover::~SaneMinMover() = default;

// constructor with arguments
SaneMinMover::SaneMinMover(
	core::kinematics::MoveMapOP movemap_in,
	core::scoring::ScoreFunctionOP scorefxn_in,
	core::optimization::MinimizerOptionsOP min_options_in,
	bool cartesian_in
) : protocols::moves::Mover("SaneMinMover"),
	cartesian_(cartesian_in),
	movemap_(std::move(movemap_in)),
	scorefxn_(std::move(scorefxn_in)),
	min_options_(std::move(min_options_in))
{}

bool SaneMinMover::cartesian() const {
	return cartesian_;
}

core::kinematics::MoveMapOP SaneMinMover::move_map() const {
	return movemap_;
}

core::scoring::ScoreFunctionOP SaneMinMover::score_function() const {
	return scorefxn_;
}

core::optimization::MinimizerOptionsOP SaneMinMover::min_options() const {
	return min_options_;
}

void
SaneMinMover::apply( core::pose::Pose & pose ) {
	using namespace core::optimization;

	PROF_START( basic::MINMOVER_APPLY );
	(*scorefxn_)(pose);
	if ( cartesian() ) {
		CartesianMinimizer minimizer;
		minimizer.run( pose, *movemap_, *scorefxn_, *min_options_ );
	} else {
		AtomTreeMinimizer minimizer;
		minimizer.run( pose, *movemap_, *scorefxn_, *min_options_ );
	}
	PROF_STOP( basic::MINMOVER_APPLY );

	scorefxn_->show(TR, pose);
	TR.flush();
}

// XRW TEMP std::string
// XRW TEMP SaneMinMover::get_name() const {
// XRW TEMP  return SaneMinMover::mover_name();
// XRW TEMP }

protocols::moves::MoverOP SaneMinMover::clone() const { return protocols::moves::MoverOP( new protocols::minimization_packing::SaneMinMover( *this ) ); }

void SaneMinMover::set_defaults_() {
	cartesian_   = false;
	movemap_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap );
	scorefxn_    = core::scoring::get_score_function();
	min_options_ = core::optimization::MinimizerOptionsOP( new core::optimization::MinimizerOptions( "lbfgs_armijo_nonmonotone", 1e-2, true, false, false ) );
}

std::string SaneMinMover::get_name() const {
	return mover_name();
}

std::string SaneMinMover::mover_name() {
	return "SaneMinMover";
}

void SaneMinMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Minimizes a Pose to a local energy minimum by "
		"performing energy minimization of a ScoreFunction over the allowable degrees "
		"of freedom, defined by a MoveMap. Unlike the classic MinMover, the only "
		"method for setting Minimization options is via the MinimizerOptions class.",
		attlist);
}

std::string SaneMinMoverCreator::keyname() const {
	return SaneMinMover::mover_name();
}

protocols::moves::MoverOP
SaneMinMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SaneMinMover );
}

void SaneMinMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SaneMinMover::provide_xml_schema( xsd );
}


} // moves
} // protocols
