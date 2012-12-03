// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/OperatorFilter.cc
/// @brief
/// @author Gabi Pszolla & Sarel Fleishman


//Unit Headers
#include <protocols/filters/OperatorFilter.hh>
#include <protocols/filters/OperatorFilterCreator.hh>
#include <protocols/filters/SigmoidFilter.hh>
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
namespace filters {

static basic::Tracer TR( "protocols.filters.Operator" );

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
Operator::reset_baseline( core::pose::Pose const & pose, bool const attempt_read_from_checkpoint ){
	foreach( FilterOP comp_statement_filt, filters() ){///all RosettaScripts user-defined filters are compoundstatements
    CompoundFilterOP comp_filt_op( dynamic_cast< CompoundFilter * >( comp_statement_filt() ) );
    runtime_assert( comp_filt_op );
    for( CompoundFilter::CompoundStatement::iterator cs_it = comp_filt_op->begin(); cs_it != comp_filt_op->end(); ++cs_it ){
       protocols::filters::FilterOP f( cs_it->first );
			if( f->get_type() == "Sigmoid" ){
				SigmoidOP sigmoid_filter( dynamic_cast< Sigmoid * >( f() ) );
				runtime_assert( sigmoid_filter );
				sigmoid_filter->reset_baseline( pose, attempt_read_from_checkpoint );
				TR<<"Resetting Sigmoid filter's baseline"<<std::endl;
			}
			else if( f->get_type() == "Operator" ){ //recursive call
				OperatorOP operator_filter( dynamic_cast< Operator * >( f() ) );
				runtime_assert( operator_filter );
				operator_filter->reset_baseline( pose, attempt_read_from_checkpoint );
				TR<<"Resetting Operator filter's baseline"<<std::endl;
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
	if( op=="SUBTRACT" )
		operation( SUBTRACT );
	if( op=="ABS" )
		operation( ABS );
	if( op != "SUM" && op != "PRODUCT" && op != "NORMALIZED_SUM" && op != "MAX" && op != "MIN" && op != "SUBTRACT" && op != "ABS" )
		utility_exit_with_message( "Operation " + op + " not recognized" );
	threshold( tag->getOption< core::Real >( "threshold", 0 ) );

  utility::vector1< std::string > const filter_names( utility::string_split( tag->getOption< std::string >( "filters" ), ',' ) );
	foreach( std::string const fname, filter_names ){
		add_filter( protocols::rosetta_scripts::parse_filter( fname, filters ) );
		TR<<"Adding filter "<<fname<<std::endl;
	}
	if( operation() == SUBTRACT && filters_.size() != 2 )
		utility_exit_with_message( "Operation SUBTRACT requested, but the number of filters provided is different than 2. I only know how to subtract one filter from another" );
	if( operation() == ABS && filters_.size() != 1 )
		utility_exit_with_message( "Operation ABS requested, but the number of filters provided is different than 1. I only know how to take the absolute value of one filter" );
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
	core::Real const val( compute( pose ) );
	TR<<"Operator returning: "<<val<<std::endl;
	return( val );
}

core::Real
Operator::compute(
	core::pose::Pose const & pose
) const {
	core::Real val( 0.0 );
	if( operation() == SUBTRACT ){
		runtime_assert( filters().size() == 2 );
		core::Real const val1( filters()[ 1 ]->report_sm( pose ) );
		core::Real const val2( filters()[ 2 ]->report_sm( pose ) );
		core::Real const difference( val1 - val2 );
		TR<<"Filters' values "<<val1<<" "<<val2<<" and the difference is "<<difference<<std::endl;
		return difference;
	}
	if( operation() == ABS ){
		runtime_assert( filters().size() == 1 );
		core::Real const val( filters()[ 1 ]->report_sm( pose ) );
		core::Real const abs_val( std::abs( val ) );
		TR<<"Filter returns "<<val<<" and its absolute value is "<<abs_val<<std::endl;
		return abs_val;
	}
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
