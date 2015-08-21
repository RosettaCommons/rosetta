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
#include <core/kinematics/Jump.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/docking/metrics.hh>


#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/MetricValue.hh>
// Project Headers


static thread_local basic::Tracer TR( "devel.replica_docking.IrmsdFilter" );

namespace devel {
namespace replica_docking {

protocols::filters::FilterOP
IrmsdFilterCreator::create_filter() const { return protocols::filters::FilterOP( new IrmsdFilter ); }

std::string
IrmsdFilterCreator::keyname() const { return "Irms"; }

void IrmsdFilter::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( in::file::native );
}


IrmsdFilter::IrmsdFilter() :
	Filter( "Irms" ),
	lower_threshold_( 0.0 ),
	upper_threshold_(9999)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( option[ in::file::native ].user() ) {
		core::pose::PoseOP native_pose( new core::pose::Pose() );
		core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]());
		native_pose_ = native_pose;
	} else {
		utility_exit_with_message("need to specify native pdb to calculate Irms");
	}
	scorefxn_ = core::scoring::get_score_function();
	// scorefxn_->show(TR.Info);
	movable_jumps_ = utility::tools::make_vector1<core::Size>(1);
	TR << "End constructer"<<std::endl;

}

IrmsdFilter::IrmsdFilter( core::scoring::ScoreFunctionOP sfxn, core::Size const rb_jump,core::Real const lower_threshold, core::Real const upper_threshold ) :
	Filter( "Irms" ),
	lower_threshold_( lower_threshold ),
	upper_threshold_(upper_threshold)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( option[ in::file::native ].user() ) {
		core::pose::PoseOP native_pose( new core::pose::Pose() );
		core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]);
		native_pose_ = native_pose;
	} else {
		utility_exit_with_message("need to specify native pdb to calculate Irms");
	}

	if ( !sfxn ) {
		scorefxn_ = core::scoring::get_score_function();
	} else {
		scorefxn_ = sfxn->clone();
	}
	TR.Info <<"IrmsdEvaluator: "<<"score" << std::endl;
	// scorefxn_->show(TR.Info);
	movable_jumps_.push_back( rb_jump );
	TR << "End constructer"<<std::endl;

}

IrmsdFilter::~IrmsdFilter(){}

protocols::filters::FilterOP
IrmsdFilter::clone() const{
	return protocols::filters::FilterOP( new IrmsdFilter( *this ) );
}

protocols::filters::FilterOP
IrmsdFilter::fresh_instance() const{
	return protocols::filters::FilterOP( new IrmsdFilter );
}

void
IrmsdFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	std::string const scorefxn_name(
		protocols::rosetta_scripts::get_score_function_name(tag) );
	scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( scorefxn_name );
	// //  scorefxn_ = new core::scoring::ScoreFunction( *(data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name )) );

	// // scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );

	lower_threshold_ = tag->getOption<core::Real>( "threshold", 0.0 );
	upper_threshold_ = tag->getOption<core::Real>( "upper_threshold", 9999);
	jump( tag->getOption< core::Size >( "jump", 1 ));

}

bool
IrmsdFilter::apply( core::pose::Pose const & pose ) const {

	core::Real const Irms( compute( pose ) );

	TR<<"Irms is "<<Irms<<". ";
	if ( Irms >= lower_threshold_ && Irms <= upper_threshold_ ) {
		TR<<"passing." <<std::endl;
		return true;
	} else {
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
IrmsdFilter::jump( core::Size const jump_id )
{
	movable_jumps_.push_back( jump_id );
}


core::Real
IrmsdFilter::compute( core::pose::Pose const & pose ) const {
	TR<<"compute Irms"<< std::endl;

	core::Real irms = protocols::docking::calc_Irmsd( pose, *native_pose_, scorefxn_, movable_jumps_ );
	return( irms );
}

}
}
