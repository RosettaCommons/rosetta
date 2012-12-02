// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/OperatorFilter.cc
/// @brief
/// @author Gabi Pszolla & Sarel Fleishman


//Unit Headers
#include <protocols/simple_filters/OperatorFilter.hh>
#include <protocols/simple_filters/OperatorFilterCreator.hh>
#include <protocols/simple_filters/SigmoidFilter.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <utility/string_util.hh>
#include <protocols/filters/BasicFilters.hh>
namespace protocols{
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.Operator" );

protocols::filters::FilterOP
OperatorFilterCreator::create_filter() const { return new Operator; }

std::string
OperatorFilterCreator::keyname() const { return "Operator"; }

//default ctor
Operator::Operator() :
protocols::filters::Filter( "Operator" ),
operation_( PRODUCT ),
threshold_( 0.0 )
{
	filters_.clear();
}

Operator::~Operator() {}

void
Operator::reset_baseline( core::pose::Pose const & pose ){
	foreach( protocols::filters::FilterOP comp_statement_filt, filters() ){///all RosettaScripts user-defined filters are compoundstatements
    protocols::filters::CompoundFilterOP comp_filt_op( dynamic_cast< protocols::filters::CompoundFilter * >( comp_statement_filt() ) );
    runtime_assert( comp_filt_op );
    for( protocols::filters::CompoundFilter::CompoundStatement::iterator cs_it = comp_filt_op->begin(); cs_it != comp_filt_op->end(); ++cs_it ){
       protocols::filters::FilterOP f( cs_it->first );
			if( f->get_type() == "Sigmoid" ){
				SigmoidOP sigmoid_filter( dynamic_cast< Sigmoid * >( f() ) );
				runtime_assert( sigmoid_filter );
				sigmoid_filter->reset_baseline( pose );
				TR<<"Resetting Sigmoid filter's baseline"<<std::endl;
			}
		}//for cs_it
	}//foreach
}

void
Operator::parse_my_tag( utility::tag::TagPtr const tag, moves::DataMap &, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & )
{
	std::string const op( tag->getOption< std::string >( "operation" ) );
	if( op=="SUM" )
		operation( SUM );
	if( op=="PRODUCT" )
		operation( PRODUCT );
	if( op=="NORMALIZED_SUM" )
		operation( NORMALIZED_SUM );
	if( op=="MAX" )
		operation( MAX );
	if( op=="MIN" )
		operation( MIN );
	if( op != "SUM" && op != "PRODUCT" && op != "NORMALIZED_SUM" && op != "MAX" && op != "MIN" )
		utility_exit_with_message( "Operation " + op + " not recognized" );
	threshold( tag->getOption< core::Real >( "threshold", 0 ) );

  utility::vector1< std::string > const filter_names( utility::string_split( tag->getOption< std::string >( "filters" ), ',' ) );
	foreach( std::string const fname, filter_names ){
		add_filter( protocols::rosetta_scripts::parse_filter( fname, filters ) );
		TR<<"Adding filter "<<fname<<std::endl;
	}
	TR<<" using operator "<<op<<std::endl;
}

bool
Operator::apply( core::pose::Pose const & pose ) const {
	core::Real const val ( compute( pose ) );
	TR<<"Filter returns "<<val<<std::endl;
	return( val >= threshold() );
}

void
Operator::report( std::ostream &o, core::pose::Pose const & pose ) const {
	core::Real const val = compute( pose );
	o << "Operator returns "<<val<<std::endl;
}

core::Real
Operator::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
Operator::compute(
	core::pose::Pose const & pose
) const {
	core::Real val( 0.0 );
	if( operation() == PRODUCT )
		val = 1.0;
	if( operation() == MIN )
		val = 9999999.999;
	if( operation() == MAX )
		val = -99999999.999;
	foreach( protocols::filters::FilterOP f, filters() ){
		core::Real const filter_val( f()->report_sm( pose ) );
		if( operation() == SUM || operation() == NORMALIZED_SUM )
			val += filter_val;
		if( operation() == PRODUCT )
			val *= filter_val;
		if( operation() == MIN ){
			if( filter_val <= val )
				val = filter_val;
		}
		if( operation() == MAX ){
			if( filter_val >= val )
				val = filter_val;
		}
	}
	if( operation() == NORMALIZED_SUM )
		val /= (core::Real) filters().size();
  return( val );
}

utility::vector1< protocols::filters::FilterOP >
Operator::filters() const{ return filters_; }

void
Operator::add_filter( protocols::filters::FilterOP f ){ filters_.push_back( f ); }
}
}
