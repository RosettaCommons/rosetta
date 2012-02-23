// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/DdgFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#include <protocols/simple_filters/DdgFilter.hh>
#include <protocols/simple_filters/DdgFilterCreator.hh>

#include <protocols/filters/Filter.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/scoring/Interface.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/Energies.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <protocols/simple_moves/ddG.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/rigid/RigidBodyMover.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.DdgFilter" );

protocols::filters::FilterOP
DdgFilterCreator::create_filter() const { return new DdgFilter; }

std::string
DdgFilterCreator::keyname() const { return "Ddg"; }


DdgFilter::DdgFilter() :
	filters::Filter( "Ddg" ),
	ddg_threshold_( -15.0 ),
	scorefxn_( NULL ),
	rb_jump_( 1 ),
	repeats_( 1 ),
	symmetry_(false),
	repack_( true ),
	relax_mover_( NULL )
{}

DdgFilter::~DdgFilter() {}

void
DdgFilter::repack( bool const repack )
{
	repack_ = repack;
}

bool
DdgFilter::repack() const
{
	return repack_;
}

void
DdgFilter::parse_my_tag( utility::tag::TagPtr const tag, moves::DataMap & data, filters::Filters_map const & , moves::Movers_map const & movers, core::pose::Pose const & )
{
	using namespace core::scoring;

	std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn", "score12" ) );
	scorefxn_ = new ScoreFunction( *(data.get< ScoreFunction * >( "scorefxns", scorefxn_name )) );
	ddg_threshold_ = tag->getOption<core::Real>( "threshold", -15 );
	rb_jump_ = tag->getOption< core::Size >( "jump", 1 );
	repeats( tag->getOption< core::Size >( "repeats", 1 ) );
	repack( tag->getOption< bool >( "repack", 1 ) );
	symmetry_ = tag->getOption<bool>( "symmetry", 0 );
	if( tag->hasOption( "relax_mover" ) )
		relax_mover( protocols::rosetta_scripts::parse_mover( tag->getOption< std::string >( "relax_mover" ), movers ) );


	if( repeats() > 1 && !repack() )
		utility_exit_with_message( "ERROR: it doesn't make sense to have repeats if repack is false, since the values converge very well." );

	if ( symmetry_ )
		TR<<"ddg filter with threshold "<< ddg_threshold_<<" repeats="<<repeats()<<" and scorefxn "<<scorefxn_name<<" with symmetry " <<std::endl;
	else
		TR<<"ddg filter with threshold "<< ddg_threshold_<<" repeats="<<repeats()<<" and scorefxn "<<scorefxn_name<<" over jump "<<rb_jump_<<" and repack "<<repack()<<std::endl;
}

bool
DdgFilter::apply( core::pose::Pose const & pose ) const
{
	core::Real const pose_ddg( compute( pose ) );
	TR<<"ddg is "<<pose_ddg<<" ";
	if( pose_ddg <= ddg_threshold_ ) {
		TR<<"passing"<<std::endl;
		return true;
	}
	TR<<"failing"<<std::endl;
	TR.flush();
	return false;
}

void
DdgFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const pose_ddg( compute( pose ) );
	out<<"ddg "<<pose_ddg<<'\n';
}

core::Real
DdgFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const pose_ddg( compute( pose ) );
	return( pose_ddg );
}

core::Size
DdgFilter::repeats() const
{
	return( repeats_ );
}

void
DdgFilter::repeats( core::Size const repeats )
{
	repeats_ = repeats;
}

core::Real
DdgFilter::compute( core::pose::Pose const & pose ) const {
	if( repack() ){
		protocols::simple_moves::ddG ddg( scorefxn_, rb_jump_, symmetry_ );
		ddg.relax_mover( relax_mover() );
		core::Real average( 0.0 );
		for( core::Size i = 1; i<=repeats_; ++i ){
			ddg.calculate( pose );
			average += ddg.sum_ddG();
			ddg.report_ddG( TR );
		}
		return average / (core::Real)repeats_;
	}
	else{
		if( repeats() > 1 && !repack() )
			utility_exit_with_message( "ERROR: it doesn't make sense to have repeats if repack is false, since the values converge very well." );
		using namespace protocols::moves;

		simple_filters::ScoreTypeFilter const stf( scorefxn_, core::scoring::total_score, 10000/*threshold*/ );
		core::pose::Pose split_pose( pose );
		rigid::RigidBodyTransMoverOP translate( new rigid::RigidBodyTransMover( split_pose, rb_jump_ ) );
		translate->step_size( 1000.0 );
		translate->apply( split_pose );
		core::Real const bound_energy( stf.compute( pose ));
		core::Real const unbound_energy( stf.compute( split_pose ));
		core::Real const dG( bound_energy - unbound_energy );
		return( dG );
	}
}



DdgFilter::DdgFilter( core::Real const ddg_threshold, core::scoring::ScoreFunctionCOP scorefxn, core::Size const rb_jump/*=1*/, core::Size const repeats/*=1*/, bool const symmetry /*=false*/ ) :
	Filter("Ddg" ),
	repack_( true ),
	relax_mover_( NULL )
{
	ddg_threshold_ = ddg_threshold;
	scorefxn_ = scorefxn->clone();
	rb_jump_ = rb_jump;
	repeats_ = repeats;
	symmetry_ = symmetry;
}

void
DdgFilter::relax_mover( protocols::moves::MoverOP m ){
	relax_mover_ = m;
}

protocols::moves::MoverOP
DdgFilter::relax_mover() const{ return relax_mover_; }


}
}
