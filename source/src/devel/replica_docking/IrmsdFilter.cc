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

#include <devel/replica_docking/IrmsdFilter.hh>
#include <devel/replica_docking/IrmsdFilterCreator.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <basic/MetricValue.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/DataMap.hh>
#include <utility/vector1.hh>

#include <core/import_pose/import_pose.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Project Headers
#include <core/types.hh>

static basic::Tracer TR( "devel.replica_docking.IrmsdFilter" );

namespace devel {
namespace replica_docking {

protocols::filters::FilterOP
IrmsdFilterCreator::create_filter() const { return new IrmsdFilter; }

std::string
IrmsdFilterCreator::keyname() const { return "Irms"; }


IrmsdFilter::IrmsdFilter() :
	Filter( "Irms" ),
	lower_threshold_( 0.0 ),
	upper_threshold_(9999),
	jump_( 1 )
{
	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		core::import_pose::pose_from_pdb( *native_pose_, basic::options::option[ basic::options::OptionKeys::inf::file::native ]);
	} else {
		utility_exit_with_message("need to specify native pdb to calculate Irms");
	}
	scorefxn_ = core::scoring::getScoreFunction();

}

IrmsdFilter::IrmsdFilter( core::scoring::ScoreFunctionOP sfxn, core::Size const rb_jump,core::Real const lower_threshold, core::Real const upper_threshold ) :
	Filter( "Irms" ),
	lower_threshold_( lower_threshold ),
	upper_threshold_(upper_threshold),
	jump_( rb_jump )
{
	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		core::import_pose::pose_from_pdb( *native_pose_, basic::options::option[ basic::options::OptionKeys::inf::file::native ]);
	} else {
		utility_exit_with_message("need to specify native pdb to calculate Irms");
	}

	if( !sfxn ) {
		scorefxn_ = core::scoring::getScoreFunction();
	} else {
		scorefxn_ = sfxn->clone();
	}
	tr.Info <<"IrmsdEvaluator: "<<"score" << std::endl;
	scorefxn_->show(tr.Info);

}

IrmsdFilter::~IrmsdFilter(){}

protocols::filters::FilterOP
IrmsdFilter::clone() const{
	return new IrmsdFilter( *this );
}

protocols::filters::FilterOP
IrmsdFilter::fresh_instance() const{
	return new IrmsdFilter;
}

void
IrmsdFilter::parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose )
{

//  	std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn", "score12" ) );
// 	scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( scorefxn_name );
// // 	scorefxn_ = new core::scoring::ScoreFunction( *(data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name )) );

// //	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );

	lower_threshold_ = tag->getOption<core::Real>( "threshold", 0.0 );
	upper_threshold_ = tag->getOption<core::Real>( "upper_threshold", 9999);
	jump( tag->getOption< core::Size >( "jump", 1 ));

}

bool
IrmsdFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const Irms( compute( pose ) );

	TR<<"Irms is "<<Irms<<". ";
	if( Irms >= lower_threshold_ && Irms <= upper_threshold_ ){
		TR<<"passing." <<std::endl;
		return true;
	}
	else {
		TR<<"failing."<<std::endl;
		return false;
	}
}

void
IrmsdFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const Irms( compute( pose ));
	out<<"Irms= "<< Irms<<'\n';
}

core::Real
IrmsdFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const Irms( compute( pose ));
	return( Irms );
}

void
IrmsdFilter::jump( core::Size const jump )
{
	jump_ = jump;
}

core::Size
IrmsdFilter::jump() const
{
	return jump_;
}

core::Real
IrmsdFilter::compute( core::pose::Pose const & pose ) const {

	irms = protocols::docking::Calc_Irms( pose, *native_pose_, scorefxn_, jump_ );

	return( interaction_energy );
}

}
}
