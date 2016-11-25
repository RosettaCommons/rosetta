// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/canonical_sampling/TrialCounterObserver.cc
/// @brief protocols::canonical_sampling::TrialCounterObserver methods implemented
/// @author


// Unit Headers
#include <protocols/canonical_sampling/TrialCounterObserver.hh>
#include <protocols/canonical_sampling/TrialCounterObserverCreator.hh>

#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
// Package Headers

// Project Headers
#include <core/pose/Pose.hh>


// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <cmath>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer tr( "protocols.canonical_sampling.TrialCounter" );

namespace protocols {
namespace canonical_sampling {


// XRW TEMP std::string
// XRW TEMP TrialCounterObserverCreator::keyname() const {
// XRW TEMP  return TrialCounterObserver::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP TrialCounterObserverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new TrialCounterObserver );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP TrialCounterObserver::mover_name() {
// XRW TEMP  return "TrialCounterObserver";
// XRW TEMP }


TrialCounterObserver::TrialCounterObserver(
) : ThermodynamicObserver()
{
	Mover::type( "TrialCounterObserver" );
}

TrialCounterObserver::~TrialCounterObserver() = default;

// XRW TEMP std::string TrialCounterObserver::get_name() const {
// XRW TEMP  return "TrialCounterObserver";
// XRW TEMP }

protocols::moves::MoverOP
TrialCounterObserver::clone() const {
	return protocols::moves::MoverOP( new TrialCounterObserver( *this ) );
}


void
TrialCounterObserver::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) { //no options ...
	file_ = tag->getOption< std::string >("file","trial.stats");
	io_stride_ = tag->getOption< core::Size >("stride", 10000 );
}


void
TrialCounterObserver::initialize_simulation(
	core::pose::Pose & /*pose*/,
	MetropolisHastingsMover const& mhm, /*metropolis_hastings_mover*/
	core::Size //default=0; non-zero if trajectory is restarted
)
{
	counters_.set_temperature_observer( mhm.tempering().get() );
	counters_.reset();
}

void
TrialCounterObserver::observe_after_metropolis(
	MetropolisHastingsMover const & mhm
) {
	std::string const& move_type( mhm.last_move().type() );
	counters_.count_trial( move_type );
	if ( mhm.last_accepted() ) {
		counters_.count_accepted( move_type );
	}
	if ( mhm.current_trial() && mhm.current_trial() % io_stride_ == 0 ) {
#ifdef WIN32
		core::Size output_ct( floor( (double)mhm.current_trial() / (double)io_stride_ ) );
#else
		core::Size output_ct( floor( mhm.current_trial() / io_stride_ ) );
#endif
		counters_.write_to_file( file_, mhm.output_name() + utility::to_string( output_ct ) );
	}
}

void
TrialCounterObserver::finalize_simulation(
	core::pose::Pose & /*pose*/,
	MetropolisHastingsMover const &mhm /*metropolis_hastings_mover*/
)
{
	counters_.show( tr.Info );
	counters_.write_to_file( file_, mhm.output_name() );
}

std::string TrialCounterObserver::get_name() const {
	return mover_name();
}

std::string TrialCounterObserver::mover_name() {
	return "TrialCounterObserver";
}

void TrialCounterObserver::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "file", xs_string, "Output file for this observer", "trial.stats" )
		+ XMLSchemaAttribute::attribute_w_default( "stride", xsct_non_negative_integer, "How many steps between outputs", "10000" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "observes acceptance rates of various moves at different temperatures", attlist );
}

std::string TrialCounterObserverCreator::keyname() const {
	return TrialCounterObserver::mover_name();
}

protocols::moves::MoverOP
TrialCounterObserverCreator::create_mover() const {
	return protocols::moves::MoverOP( new TrialCounterObserver );
}

void TrialCounterObserverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TrialCounterObserver::provide_xml_schema( xsd );
}


} //moves
} //protocols

