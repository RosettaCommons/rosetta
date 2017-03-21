// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/moves/IteratedConvergenceMover.cc
/// @brief  Implementation of IteratedConvergenceMover, class to repeatedly apply a submover until filter convergence is reached
/// @author Rocco Moretti (rmoretti@u.washington.edu)

// Unit headers
#include <protocols/moves/IteratedConvergenceMover.hh>
#include <protocols/moves/IteratedConvergenceMoverCreator.hh>

#include <protocols/filters/Filter.hh>

#include <basic/Tracer.hh>

// Utility Headers
#include <utility>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

// option key includes

#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace moves {

static THREAD_LOCAL basic::Tracer TR( "protocols.moves.IteratedConvergenceMover" );

IteratedConvergenceMover::IteratedConvergenceMover() :
	Mover("IteratedConvergenceMover"),
	delta_(0.1),
	cycles_(1),
	maxcycles_(1000)
{}

IteratedConvergenceMover::IteratedConvergenceMover( MoverOP submover, filters::FilterCOP filter, core::Real delta, core::Size cycles, core::Size maxcycles ) :
	Mover("IteratedConvergenceMover"),
	submover_(std::move(submover)),
	filter_(std::move(filter)),
	cycles_(cycles),
	maxcycles_(maxcycles)
{
	this->delta(delta);
}

IteratedConvergenceMover::~IteratedConvergenceMover()= default;

IteratedConvergenceMover::IteratedConvergenceMover( IteratedConvergenceMover const & other ) :
	//utility::pointer::ReferenceCount(),
	Mover(other),
	submover_(other.submover_),
	filter_(other.filter_),
	delta_(other.delta_),
	cycles_(other.cycles_),
	maxcycles_(other.maxcycles_)
{}

void
IteratedConvergenceMover::apply( Pose & pose )
{
	runtime_assert(submover_ != nullptr);
	runtime_assert(filter_ != nullptr);

	core::Real refval(filter_->report_sm(pose));
	core::Size ncyc(0);
	for ( core::Size ii(1); ii <= maxcycles_; ++ii ) {
		submover_->apply(pose);
		core::Real trialval(filter_->report_sm(pose));
		if ( (trialval < refval - delta_) || ( trialval > refval + delta_ ) ) {
			//Outside window, reset and continue
			TR << "After "<<ii<<" applications, reseting reference to "<<trialval<<" which is "<<(trialval-refval)<<" away from previous reference."<<std::endl;
			refval=trialval;
			ncyc=0;
			continue;
		}
		//within tolerance
		++ncyc;
		if ( ncyc >= cycles_ ) {
			TR << "Exiting IteratedConvergence with convergence after "<<ii<<" applications of mover."<<std::endl;
			return;
		}
		TR<<"After "<<ii<<" applications, at cycle "<<ncyc<<" of "<<cycles_<<", filter delta in range at "<<(trialval-refval)<<std::endl;
	}
	TR << "Exiting IteratedConvergence WITHOUT convergence after "<<maxcycles_<<" applications of mover."<<std::endl;
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
IteratedConvergenceMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const & filters,
	Movers_map const & movers,
	Pose const &
)
{
	delta( tag->getOption< core::Real >( "delta", 0.1 ) );
	cycles( tag->getOption< core::Size >( "cycles", 1 ) );
	maxcycles( tag->getOption< core::Size >( "maxcycles", 1000 ) );

	std::string mover_name("");
	std::string filter_name("");
	if ( tag->hasOption("mover") ) mover_name = tag->getOption< std::string >( "mover" );
	if ( tag->hasOption("mover_name") ) mover_name = tag->getOption< std::string >( "mover_name" );
	if ( tag->hasOption("filter") ) filter_name = tag->getOption< std::string >( "filter" );
	if ( tag->hasOption("filter_name") ) filter_name = tag->getOption< std::string >( "filter_name" );

	auto  find_mover ( movers.find( mover_name ));
	auto find_filter( filters.find( filter_name ));

	if ( find_mover == movers.end() ) {
		TR.Error << "mover not found in map: \n" << tag << std::endl;
		runtime_assert( find_mover != movers.end() );
	}
	if ( find_filter == filters.end() ) {
		TR.Error << "filter not found in map: \n" << tag << std::endl;
		runtime_assert( find_filter != filters.end() );
	}
	submover(find_mover->second);
	filter(find_filter->second);

	TR << "Setting IteratedConvergence for "<<cycles()<<" cycles (max "<<maxcycles()<<") with mover '"<<mover_name<<"' and filter '"<<filter_name<<"' and delta tolerance "<<delta()<<std::endl;
}

/// @brief required in the context of the parser/scripting scheme
MoverOP
IteratedConvergenceMover::fresh_instance() const
{
	return MoverOP( new IteratedConvergenceMover );
}

/// @brief required in the context of the parser/scripting scheme
MoverOP
IteratedConvergenceMover::clone() const
{
	return MoverOP( new IteratedConvergenceMover( *this ) );
}

// setters
void IteratedConvergenceMover::submover( MoverOP mover )
{
	runtime_assert( mover != nullptr );
	submover_ = mover;
}

void IteratedConvergenceMover::filter( filters::FilterCOP filter )
{
	runtime_assert( filter != nullptr );
	filter_ = filter;
}

void IteratedConvergenceMover::delta( core::Real delta ) {
	if ( delta < 0 ) {
		utility_exit_with_message("Delta value given to IteratedConvergenceMover must be non-negative");
	}
	delta_ = delta;
}

void IteratedConvergenceMover::cycles( core::Size cycles ) { cycles_ = cycles; }
void IteratedConvergenceMover::maxcycles( core::Size maxcycles ) { maxcycles_ = maxcycles; }

// accessors
MoverOP IteratedConvergenceMover::submover() const { return submover_; }
filters::FilterCOP IteratedConvergenceMover::filter() const { return filter_; }
core::Real IteratedConvergenceMover::delta() const { return delta_; }
core::Size IteratedConvergenceMover::cycles() const { return cycles_; }
core::Size IteratedConvergenceMover::maxcycles() const { return maxcycles_; }

std::string IteratedConvergenceMover::get_name() const {
	return mover_name();
}

std::string IteratedConvergenceMover::mover_name() {
	return "IteratedConvergence";
}

void IteratedConvergenceMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute(
		"mover", xs_string,
		"the mover to repeatedly apply");
	attlist + XMLSchemaAttribute(
		"mover_name", xs_string,
		"the mover to repeatedly apply");
	attlist + XMLSchemaAttribute(
		"filter", xs_string,
		"the filter to use when assaying for convergence "
		"(should return a reasonable value from report_sm())");
	attlist + XMLSchemaAttribute(
		"filter_name", xs_string,
		"the filter to use when assaying for convergence "
		"(should return a reasonable value from report_sm())");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"delta", xsct_real,
		"how close do the filter values have to be to count as converged", "0.1");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"cycles", xsct_non_negative_integer,
		"for how many mover applications does the filter value have to "
		"fall within delta of the reference value before counting as converged. "
		"If the filter is outside of the range, the reference value is "
		"reset to the new filter value", "1");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"maxcycles", xsct_non_negative_integer,
		"exit regardless if filter doesn't converge within this many "
		"applications of the mover - intended mainly as a safety "
		"check to prevent infinite recursion", "1000");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Repeatedly applies a sub-mover until the given filter returns "
		"a value within the given delta for the given number of cycles",
		attlist );
}

std::string IteratedConvergenceMoverCreator::keyname() const {
	return IteratedConvergenceMover::mover_name();
}

protocols::moves::MoverOP
IteratedConvergenceMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new IteratedConvergenceMover );
}

void IteratedConvergenceMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	IteratedConvergenceMover::provide_xml_schema( xsd );
}


} // moves
} // protocols

