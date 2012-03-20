// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/BuriedUnsatHbondFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#include <protocols/simple_filters/BuriedUnsatHbondFilter.hh>
#include <protocols/simple_filters/BuriedUnsatHbondFilterCreator.hh>

#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/MetricValue.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer buried_unsat_hbond_filter_tracer( "protocols.simple_filters.BuriedUnsatHbondFilter" );

protocols::filters::FilterOP
BuriedUnsatHbondFilterCreator::create_filter() const { return new BuriedUnsatHbondFilter; }

std::string
BuriedUnsatHbondFilterCreator::keyname() const { return "BuriedUnsatHbonds"; }

BuriedUnsatHbondFilter::BuriedUnsatHbondFilter( core::Size const upper_threshold, core::Size const jump_num ) :
	Filter( "BuriedUnsatHbonds" ),
	upper_threshold_( upper_threshold ),
	jump_num_( jump_num )
{ }

BuriedUnsatHbondFilter::~BuriedUnsatHbondFilter(){}

void
BuriedUnsatHbondFilter::parse_my_tag( utility::tag::TagPtr const tag, moves::DataMap & datamap, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	jump_num_ = tag->getOption<core::Size>( "jump_number", 1 );
	upper_threshold_ = tag->getOption<core::Size>( "cutoff", 20 );

	std::string const scorefxn_key( tag->getOption<std::string>("scorefxn") );
	if ( datamap.has( "scorefxns", scorefxn_key ) ) {
		sfxn_ = datamap.get< core::scoring::ScoreFunction* >( "scorefxns", scorefxn_key );
	} else {
		sfxn_ = core::scoring::getScoreFunction();
	}

	buried_unsat_hbond_filter_tracer<<"Buried Unsatisfied Hbond filter over jump number " << jump_num_ << " with cutoff " << upper_threshold_ << std::endl;
}

bool
BuriedUnsatHbondFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const unsat_hbonds( compute( pose ) );

	buried_unsat_hbond_filter_tracer<<"# unsatisfied hbonds: "<<unsat_hbonds<<". ";
	if( unsat_hbonds <= upper_threshold_ ){
		buried_unsat_hbond_filter_tracer<<"passing." <<std::endl;
		return true;
	}
	else {
		buried_unsat_hbond_filter_tracer<<"failing."<<std::endl;
		return false;
	}
}

void
BuriedUnsatHbondFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const unsat_hbonds( compute( pose ));
	out<<"# unsatisfied hbonds: "<< unsat_hbonds<<'\n';
}

core::Real
BuriedUnsatHbondFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const unsat_hbonds( compute( pose ));
	return( unsat_hbonds );
}

core::Real
BuriedUnsatHbondFilter::compute( core::pose::Pose const & pose ) const {

	runtime_assert( jump_num_ <= pose.num_jump() );

	// score the pose to populate the 10A neighborgraph
	core::pose::Pose bound( pose );

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
	( *scorefxn )( bound );
	bound.update_residue_neighbors();

	core::pose::Pose unbound( bound );
	if( jump_num_ ) {
		core::Real const unbound_dist = 1000.0;
		protocols::rigid::RigidBodyTransMover trans_mover( unbound, jump_num_ );
		trans_mover.trans_axis( trans_mover.trans_axis() );
		trans_mover.step_size(unbound_dist);
		trans_mover.apply( unbound );
		//unbound.update_residue_neighbors();
		(*sfxn_)(unbound ); // score the new pose, or we get assertion erros.
	}

	basic::MetricValue< core::Size > mv_bound, mv_unbound;
	basic::MetricValue< std::string > mv_bound_str, mv_unbound_str;

	using namespace protocols::toolbox::pose_metric_calculators;
	// Despite the name, it's counting H-bonders, not any old polars.
	BuriedUnsatisfiedPolarsCalculator calc_bound("default", "default"), calc_unbound("default", "default");
	calc_bound.get("all_bur_unsat_polars", mv_bound, bound);
	std::string bound_string = calc_bound.get( "residue_bur_unsat_polars", bound );
	buried_unsat_hbond_filter_tracer << "BOUND: " << bound_string << std::endl;

	core::Real unsat_hbonds( 0.0 );
	if( jump_num_ ) {
		calc_unbound.get("all_bur_unsat_polars", mv_unbound, unbound);
		unsat_hbonds = mv_bound.value() - mv_unbound.value();
		buried_unsat_hbond_filter_tracer << "unbound_unsat=" << mv_unbound.value() << "    " << "bound_unsat=" << mv_bound.value() << std::endl;
		std::string unbound_string = calc_unbound.get( "residue_bur_unsat_polars", unbound );
		buried_unsat_hbond_filter_tracer << "UNBOUND: " << unbound_string << std::endl;
	}
	else unsat_hbonds = mv_bound.value();


	return( unsat_hbonds );
}
}
}
