// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Zhe Zhang

#include <devel/replica_docking/InteractionScoreFilter.hh>
#include <devel/replica_docking/InteractionScoreFilterCreator.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <basic/MetricValue.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/vector1.hh>

// Project Headers
#include <utility/excn/Exceptions.hh>
#include <core/types.hh>

namespace devel {
namespace replica_docking {

static basic::Tracer TR( "devel.replica_docking.InteractionScoreFilter" );

protocols::filters::FilterOP
InteractionScoreFilterCreator::create_filter() const { return new InteractionScoreFilter; }

std::string
InteractionScoreFilterCreator::keyname() const { return "I_sc"; }


InteractionScoreFilter::InteractionScoreFilter() :
	Filter( "I_sc" ),
	lower_threshold_( -30.0 ),
	upper_threshold_(0.0),
	jump_( 1 ),
	scorefxn_( NULL )
{
	scorefxn_ = core::scoring::getScoreFunction();
}

InteractionScoreFilter::InteractionScoreFilter( std::string const scorefxn_name, core::Size const rb_jump,core::Real const lower_threshold, core::Real const upper_threshold ) :
	Filter( "I_sc" ),
	lower_threshold_( lower_threshold ),
	upper_threshold_(upper_threshold),
	jump_( rb_jump )
{
	//	scorefxn_->initialize_from_file( scorefxn_name );
 	scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( scorefxn_name );
}

InteractionScoreFilter::InteractionScoreFilter( core::scoring::ScoreFunctionCOP scorefxn, core::Size const rb_jump, core::Real const lower_threshold, core::Real const upper_threshold ) :
	Filter( "I_sc" ),
	lower_threshold_( lower_threshold ),
	upper_threshold_(upper_threshold),
	jump_( rb_jump )
{
	scorefxn_ = scorefxn->clone();
}

InteractionScoreFilter::~InteractionScoreFilter(){}

protocols::filters::FilterOP
InteractionScoreFilter::clone() const{
	return new InteractionScoreFilter( *this );
}

protocols::filters::FilterOP
InteractionScoreFilter::fresh_instance() const{
	return new InteractionScoreFilter;
}

void
InteractionScoreFilter::parse_my_tag( utility::tag::TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose )
{
 	std::string const scorefxn_name(
		protocols::rosetta_scripts::get_score_function_name(tag) );
	scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( scorefxn_name );
// 	scorefxn_ = new core::scoring::ScoreFunction( *(data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name )) );

//	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );

	lower_threshold_ = tag->getOption<core::Real>( "threshold", -30 );
	upper_threshold_ = tag->getOption<core::Real>( "upper_threshold", 0);
	jump( tag->getOption< core::Size >( "jump", 1 ));

	if( !pose.is_fullatom() )
		throw utility::excn::EXCN_RosettaScriptsOption( "ERROR: it doesn't make sense to calculate the interaction score on low-res pose since in I_sc=bound-A-B, A&B are constant.\n It is totally a waste of time" );

	TR<<"InterfaceScoreFilter with lower threshold of "<<lower_threshold_<<" and jump "<<jump()<<'\n';
	TR.flush();
}

bool
InteractionScoreFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const I_sc( compute( pose ) );

	TR<<"I_sc is "<<I_sc<<". ";
	if( I_sc >= lower_threshold_ && I_sc <= upper_threshold_ ){
		TR<<"passing." <<std::endl;
		return true;
	}
	else {
		TR<<"failing."<<std::endl;
		return false;
	}
}

void
InteractionScoreFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const I_sc( compute( pose ));
	out<<"I_sc= "<< I_sc<<'\n';
}

core::Real
InteractionScoreFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const I_sc( compute( pose ));
	return( I_sc );
}

void
InteractionScoreFilter::jump( core::Size const jump )
{
	jump_ = jump;
}

core::Size
InteractionScoreFilter::jump() const
{
	return jump_;
}

core::Real
InteractionScoreFilter::compute( core::pose::Pose const & pose ) const {

	core::pose::Pose split_pose( pose );
	scorefxn_->set_weight( core::scoring::atom_pair_constraint, 0.0 );
	core::Real bound_energy = (*scorefxn_)( split_pose );

	core::kinematics::Jump bound_pose_jump = split_pose.jump( jump_ );
	core::kinematics::Jump unbound_pose_jump( bound_pose_jump );
	core::Real trans = 10000;
	core::kinematics::Stub upstream_stub = split_pose.conformation().upstream_jump_stub( jump_ );

	core::Vector dummy_axis(1,1,1);
	unbound_pose_jump.translation_along_axis( upstream_stub, dummy_axis, trans );
	split_pose.set_jump( jump_, unbound_pose_jump );
	core::Real unbound_energy = (*scorefxn_)( split_pose );

	core::Real interaction_energy = ( bound_energy - unbound_energy );
	TR.Debug << "unbound_energy " << unbound_energy << " bound_energy " << bound_energy << " I_sc " << interaction_energy << std::endl;

	return( interaction_energy );
}

}
}
