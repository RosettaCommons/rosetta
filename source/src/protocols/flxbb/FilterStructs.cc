// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/flxbb/FilterStructs.cc
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Header
#include <protocols/flxbb/FilterStructs.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility>
#include <utility/exit.hh>

#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.flxbb.FilterStructs" );


using namespace core;
using namespace protocols::flxbb;

namespace protocols {
namespace flxbb {

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief default constructor
FilterStructs::FilterStructs():
	name_( "BaseFilterStructs" ),
	filter_on_( true ),
	ntrial_( 5 ),
	current_trial_( 0 ),
	best_pose_( /* NULL */ )
{}

/// @brief value constructor
FilterStructs::FilterStructs( String  name ):
	name_(std::move( name )),
	filter_on_( true ),
	ntrial_( 5 ),
	current_trial_( 0 ),
	best_pose_( /* NULL */ )
{}

/// @brief value constructor
FilterStructs::FilterStructs( String  name, Size const ntrial ):
	name_(std::move( name )),
	filter_on_( true ),
	ntrial_( ntrial ),
	current_trial_( 0 ),
	best_pose_( /* NULL */ )
{}


/// @brief value constructor
FilterStructs::FilterStructs( String  name, Pose const & pose, Size const ntrial ):
	name_(std::move( name )),
	filter_on_( true ),
	ntrial_( ntrial ),
	current_trial_( 0 ),
	best_pose_( PoseOP( new Pose(pose) ) )
{}

/// @brief copy constructor
FilterStructs::FilterStructs( FilterStructs const & rval ):
	ReferenceCount(),
	name_( rval.name_ ),
	filter_on_( rval.filter_on_ ),
	ntrial_( rval.ntrial_ ),
	best_pose_( rval.best_pose_ )
{
	current_trial_ = 0;
}

/// @brief destructor
FilterStructs::~FilterStructs()= default;

/// @brief clone this object
FilterStructsOP
FilterStructs::clone() const
{
	utility_exit_with_message( "clone has been called on a Mover which has not overriden the base class implementation.  Probably you tried to pass a Mover to the job distributor or parser which does not have clone implemented.  Implement the function and try again.\n");
	return FilterStructsOP( nullptr );
}

/// @brief create a new instance of this type of object
FilterStructsOP
FilterStructs::fresh_instance() const
{
	utility_exit_with_message("fresh_instance has been called on a Mover which has not overriden the base class implementation.  Probably you tried to pass a Mover to the job distributor which does not have fresh_instance implemented.  Implement the function and try again.\n");
	return FilterStructsOP( nullptr );
}

/// @brief return best pose
pose::PoseOP
FilterStructs::get_bestpose() const
{
	return best_pose_;
}

/// @brief
void
FilterStructs::name( String const & name )
{
	name_ = name;
}

/// @brief set ntrial
void
FilterStructs::set_ntrial( Size const ntrial )
{
	ntrial_ = ntrial;
}

/// @brief
void
FilterStructs::initialize( Pose const & pose )
{
	set_filter_on();
	current_trial_ = 0;
	best_pose_ = PoseOP( new Pose( pose ) );
}

/// @brief set ntrial
void
FilterStructs::count_ntrial()
{
	++current_trial_;
}

/// @brief
bool
FilterStructs::filter_is_over()
{
	if ( current_trial_ >= ntrial_ ) {
		return true;
	} else {
		return false;
	}
}

/// @brief set best pose
void
FilterStructs::set_bestpose( Pose const & pose )
{
	best_pose_ = PoseOP( new Pose( pose ) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief value constructor
FilterStructs_Packstat::FilterStructs_Packstat( Size const ntrial ):
	FilterStructs( "FilterStructs_Packstat", ntrial ),
	best_packscore_( 0.0 )
{}

/// @brief value constructor
FilterStructs_Packstat::FilterStructs_Packstat( Pose const & pose, Size const ntrial ):
	FilterStructs( "FilterStructs_Packstat", pose, ntrial ),
	best_packscore_( 0.0 )
{}

/// @brief copy constructor
FilterStructs_Packstat::FilterStructs_Packstat( FilterStructs_Packstat const & )= default;

/// @brief destructor
FilterStructs_Packstat::~FilterStructs_Packstat()= default;

/// @brief clone
FilterStructsOP FilterStructs_Packstat::clone() const
{
	return FilterStructsOP( new FilterStructs_Packstat( *this ) );
}

/// @brief fresh instance
FilterStructsOP FilterStructs_Packstat::fresh_instance() const
{
	return FilterStructsOP( new FilterStructs_Packstat() );
}

/// @brief
void FilterStructs_Packstat::reset( Pose const & pose )
{
	best_packscore_ = 0.0;
	initialize( pose );
}

/// @brief filter apply
void FilterStructs_Packstat::apply( Pose const & pose )
{
	using namespace core::scoring::packstat;
	count_ntrial();
	Real packscore;
	packscore = compute_packing_score( pose );

	if (  packscore > best_packscore_ ) {
		set_bestpose( pose );
		best_packscore_ = packscore;
		TR << " Packscore : " << best_packscore_ << std::endl;
	}
	if ( filter_is_over() ) {
		set_filter_off();
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief constructor
FilterStructs_TotalCharge::FilterStructs_TotalCharge( Size const ntrial ):
	FilterStructs( "FilterStructs_TotalCharge", ntrial ),
	disallowed_value_( 1 )
{}

/// @brief constructor
FilterStructs_TotalCharge::FilterStructs_TotalCharge( Pose const & pose, Size const ntrial ):
	FilterStructs( "FilterStructs_TotalCharge", pose, ntrial ),
	disallowed_value_( 1 )
{}

/// @brief copy constructor
FilterStructs_TotalCharge::FilterStructs_TotalCharge( FilterStructs_TotalCharge const & )= default;

/// @brief destructor
FilterStructs_TotalCharge::~FilterStructs_TotalCharge() = default;

/// @brief clone
FilterStructsOP FilterStructs_TotalCharge::clone() const
{
	return FilterStructsOP( new FilterStructs_TotalCharge( *this ) );
}

/// @brief clone
FilterStructsOP FilterStructs_TotalCharge::fresh_instance() const
{
	return FilterStructsOP( new FilterStructs_TotalCharge() );
}

/// @brief
void FilterStructs_TotalCharge::reset( Pose const & pose )
{
	initialize( pose );
}

/// @brief filter apply
void FilterStructs_TotalCharge::apply( Pose const & pose )
{

	count_ntrial();
	Real total_charge( 0.0 );
	for ( Size i=1; i<=pose.size(); ++i ) {

		chemical::AA aa = pose.aa( i  );

		if ( aa == chemical::aa_from_name( "GLU" ) ) {
			total_charge -= 1.0;
		} else if ( aa == chemical::aa_from_name( "ASP" ) ) {
			total_charge -= 1.0;
		} else if ( aa == chemical::aa_from_name( "ARG" ) ) {
			total_charge += 1.0;
		} else if ( aa == chemical::aa_from_name( "LYS" ) ) {
			total_charge += 1.0;
		}
	}

	if ( total_charge != disallowed_value_ || filter_is_over() ) {
		set_filter_off();
		set_bestpose( pose );
		TR << " Total charge : " << total_charge << std::endl;
	}

}


} // namespace flxbb
} // namespace protocols
