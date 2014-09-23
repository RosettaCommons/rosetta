// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ResetBaselineMover.cc
/// @brief

// Unit headers
#include <protocols/simple_moves/ResetBaselineMover.hh>
#include <protocols/simple_moves/ResetBaselineMoverCreator.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/simple_filters/OperatorFilter.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
static thread_local basic::Tracer TR( "protocols.simple_moves.ResetBaselineMover" );
#include <utility/tag/Tag.hh>

#include <core/pose/Pose.hh>

namespace protocols {
namespace simple_moves {

std::string
ResetBaselineMoverCreator::keyname() const
{
	return ResetBaselineMoverCreator::mover_name();
}

protocols::moves::MoverOP
ResetBaselineMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ResetBaselineMover );
}

std::string
ResetBaselineMoverCreator::mover_name()
{
	return "ResetBaseline";
}

ResetBaselineMover::ResetBaselineMover()
	: moves::Mover("ResetBaseline"),
	filter_( /* NULL */ )
{}

void
ResetBaselineMover::apply( Pose & pose )
{
	using namespace protocols::filters;
	using namespace protocols::simple_filters;

  std::string const filter_type( filter()->get_type() );
  runtime_assert( filter_type == "Operator" || filter_type == "CompoundStatement" );
	if( filter_type == "Operator" ){
		TR<<"resetting Operator filter's baseline"<<std::endl;
		OperatorOP operator_filter( utility::pointer::dynamic_pointer_cast< protocols::simple_filters::Operator > ( filter() )) ;
		runtime_assert( operator_filter != 0 );
		operator_filter->reset_baseline( pose, false/*reset baselines from checkpoint*/ );
	}
	else if( filter_type == "CompoundStatement" ){
		CompoundFilterOP comp_filt_op( utility::pointer::dynamic_pointer_cast< protocols::filters::CompoundFilter > ( filter() ) );
		runtime_assert( comp_filt_op != 0 );
    for( CompoundFilter::CompoundStatement::iterator cs_it = comp_filt_op->begin(); cs_it != comp_filt_op->end(); ++cs_it ){
      FilterOP filt( cs_it->first );
      if( filt->get_type() == "Operator" ){
        TR<<"Resetting Operator filter's baseline"<<std::endl;
        OperatorOP operator_filter( utility::pointer::dynamic_pointer_cast< protocols::simple_filters::Operator > ( filt ) );
				runtime_assert( operator_filter != 0 );
        operator_filter->reset_baseline( pose, false/*read baselines from checkpoint*/ );
      }// fi Operator
    }// for cs_it
  }//elseif CompoundStatement
}

std::string
ResetBaselineMover::get_name() const {
	return ResetBaselineMoverCreator::mover_name();
}

moves::MoverOP
ResetBaselineMover::clone() const
{
	return moves::MoverOP( new ResetBaselineMover( *this ) );
}

moves::MoverOP
ResetBaselineMover::fresh_instance() const
{
	return moves::MoverOP( new ResetBaselineMover );
}

void
ResetBaselineMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &filters,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	std::string const filter_name( tag->getOption< std::string >("filter" ) );
	filter( protocols::rosetta_scripts::parse_filter( filter_name, filters ) );
	std::string const filter_type( filter()->get_type());
	TR<<"filter: "<<filter_name<<" of type: "<<filter_type<<std::endl;
}

protocols::filters::FilterOP ResetBaselineMover::filter() const{ return filter_; }

void ResetBaselineMover::filter( protocols::filters::FilterOP f ){
	filter_ = f;
	std::string const filter_type( filter_->get_type() );
	runtime_assert( filter_type == "Operator" || filter_type == "CompoundStatement" );
}

} // simple_moves
} // protocols
