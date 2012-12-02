// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/SigmoidFilter.cc
/// @brief
/// @author Gabi Pszolla & Sarel Fleishman


//Unit Headers
#include <protocols/simple_filters/SigmoidFilter.hh>
#include <protocols/simple_filters/SigmoidFilterCreator.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <protocols/rosetta_scripts/util.hh>
namespace protocols{
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.Sigmoid" );

protocols::filters::FilterOP
SigmoidFilterCreator::create_filter() const { return new Sigmoid; }

std::string
SigmoidFilterCreator::keyname() const { return "Sigmoid"; }

//default ctor
Sigmoid::Sigmoid() :
protocols::filters::Filter( "Sigmoid" ),
filter_( NULL ),
steepness_( 1.0 ),
offset_( 0.0 ),
baseline_( 0.0 ),
negate_( false ),
threshold_( 0 )
{
}

Sigmoid::~Sigmoid() {}

void
Sigmoid::reset_baseline( core::pose::Pose const & pose ){
	baseline_ = filter()->report_sm( pose );
	TR<<"Resetting baseline to: "<<baseline_;
}

void
Sigmoid::parse_my_tag( utility::tag::TagPtr const tag, moves::DataMap &, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & )
{
	steepness( tag->getOption< core::Real >( "steepness", 1.0 ) );
	offset( tag->getOption< core::Real >( "offset", 0 ));
	negate( tag->getOption< bool >( "negate", false ) );
	threshold( tag->getOption< core::Real >( "threshold", 0 ) );
	filter( protocols::rosetta_scripts::parse_filter( tag->getOption< std::string >( "filter" ), filters ) );
	TR<<"Sigmoid with options: steepness "<<steepness()<<" offset "<<offset()<<" negate "<<negate()<<" threshold "<<threshold()<<" filter: "<<tag->getOption< std::string >( "filter" ) << std::endl;
}

bool
Sigmoid::apply( core::pose::Pose const & pose ) const {
	core::Real const val ( compute( pose ) );
	return( val >= threshold() );
}

void
Sigmoid::report( std::ostream &o, core::pose::Pose const & pose ) const {
	core::Real const val = compute( pose );
	o << "Sigmoid returns "<<val<<std::endl;
}

core::Real
Sigmoid::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
Sigmoid::compute(
	core::pose::Pose const & pose
) const {
  core::Real const val( ( negate() ? -filter()->report_sm( pose ) : filter()->report_sm( pose ) ) - baseline_ );
  core::Real const transform( 1.0 / ( ( 1.0 + std::exp( ( val - offset_ ) * steepness_ ) ) ) );
  TR<<"filter val/transform: "<<val<<" "<<transform<<std::endl;
  TR<<"returning: "<<transform<<std::endl;
  return( transform );
}

protocols::filters::FilterOP
Sigmoid::filter() const{ return filter_; }

void
Sigmoid::filter( protocols::filters::FilterOP f ){ filter_ = f; }
}
}
