// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/FNatFilter.cc
/// @brief filtering on fraction of native contacts
/// @author Jacob Corn (jecorn@u.washington.edu)
#include <protocols/protein_interface_design/filters/FNatFilter.hh>
#include <protocols/protein_interface_design/filters/FNatFilterCreator.hh>
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <core/pose/PDBInfo.hh>
#include <protocols/docking/metrics.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>

#include <algorithm>
#include <list>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

FNatFilter::FNatFilter() :
  protocols::filters::Filter( "FNat" ),
  threshold_( 5.0 ),
  reference_pose_( NULL )
{}

FNatFilter::FNatFilter(protocols::docking::DockJumps const movable_jumps,
			 core::Real const threshold,
			 core::pose::PoseOP reference_pose)
  : protocols::filters::Filter( "FNat" ),
    threshold_(threshold),
    reference_pose_(reference_pose),
    movable_jumps_(movable_jumps)
{}

FNatFilter::~FNatFilter() {}

protocols::filters::FilterOP
FNatFilter::clone() const {
	return new FNatFilter( *this );
}

static thread_local basic::Tracer TR( "protocols.protein_interface_design.filters.FNatFilter" );
core::Real
FNatFilter::compute( core::pose::Pose const & pose ) const
{
  return protocols::docking::calc_Fnat(pose, *reference_pose_, scorefxn_, movable_jumps_);
}

bool
FNatFilter::apply( core::pose::Pose const & pose ) const {

	core::Real const f_nat( compute( pose ));
	TR << "f(nat): " << f_nat;
	if( f_nat <= threshold_ )
	{
		TR<<" passing."<<std::endl;
		return( true );
	}
	else TR<<" failing." << std::endl;
	return( false );
}

void
FNatFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const f_nat( compute( pose ));
	out<<"f(nat): " << f_nat <<'\n';
}

core::Real
FNatFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const f_nat( compute( pose ));
	return( (core::Real) f_nat );
}

void
FNatFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & reference_pose )
{
	/// @details
	///if the save pose mover has been instantiated, this filter can calculate the rms
	///against the ref pose
	if( tag->hasOption("reference_name") ){
		reference_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag,data_map );
	}
	else{
		reference_pose_ = new core::pose::Pose( reference_pose );
		if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() )
			core::import_pose::pose_from_pdb( *reference_pose_, basic::options::option[ basic::options::OptionKeys::in::file::native ] );
	}

	threshold_ = tag->getOption<core::Real>( "threshold", 5 );

	//TODO: support multiple jumps
	core::Size jump_num = tag->getOption<core::Size>( "jump", 1);

	//TODO: convert jump_num to movable_jumps_ (vector0?)
	movable_jumps_.push_back(jump_num);

	scorefxn_ = rosetta_scripts::parse_score_function( tag, data_map )->clone();

	TR<<"Built FNatFilter with threshold " << threshold_ << ", scorefxn " <<
		rosetta_scripts::get_score_function_name(tag) <<", jump "<< jump_num << std::endl;

}

protocols::filters::FilterOP
FNatFilterCreator::create_filter() const { return new FNatFilter; }

std::string
FNatFilterCreator::keyname() const { return "FNat"; }


} // filters
} // protein_interface_design
} // devel


