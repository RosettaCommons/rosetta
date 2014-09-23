// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/RmsdFilter.cc
/// @brief rmsd filtering
/// @author Jacob Corn (jecorn@u.washington.edu)
#include <protocols/protein_interface_design/filters/LRmsdFilter.hh>
#include <protocols/protein_interface_design/filters/LRmsdFilterCreator.hh>
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
// AUTO-REMOVED #include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <core/pose/PDBInfo.hh>
#include <protocols/docking/metrics.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>

#include <algorithm>
#include <list>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

LRmsdFilter::LRmsdFilter() :
  protocols::filters::Filter( "LRmsd" ),
  threshold_( 5.0 ),
  reference_pose_( /* NULL */ )
{}

LRmsdFilter::LRmsdFilter(protocols::docking::DockJumps const movable_jumps,
			 core::Real const threshold,
			 core::pose::PoseOP reference_pose)
  : protocols::filters::Filter( "LRmsd" ),
    threshold_(threshold),
    reference_pose_(reference_pose),
    movable_jumps_(movable_jumps)
{}

LRmsdFilter::~LRmsdFilter() {}

protocols::filters::FilterOP
LRmsdFilter::clone() const {
	return protocols::filters::FilterOP( new LRmsdFilter( *this ) );
}

static thread_local basic::Tracer TR( "protocols.protein_interface_design.filters.LRmsdFilter" );
core::Real
LRmsdFilter::compute( core::pose::Pose const & pose ) const
{
  return protocols::docking::calc_Lrmsd(pose, *reference_pose_, movable_jumps_);
}

bool
LRmsdFilter::apply( core::pose::Pose const & pose ) const {

	core::Real const rmsd( compute( pose ));
	TR << "L_rmsd: " << rmsd ;
	if( rmsd <= threshold_ )
	{
		TR<<" passing."<<std::endl;
		return( true );
	}
	else TR<<" failing." << std::endl;
	return( false );
}

void
LRmsdFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ));
	out<<"RMSD: " << rmsd <<'\n';
}

core::Real
LRmsdFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ));
	return( (core::Real) rmsd );
}

void
LRmsdFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & reference_pose )
{
	/// @details
	///if the save pose mover has been instantiated, this filter can calculate the rms
	///against the ref pose
	if( tag->hasOption("reference_name") ){
		reference_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag,data_map );
	}
	else{
		reference_pose_ = core::pose::PoseOP( new core::pose::Pose( reference_pose ) );
		if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() )
			core::import_pose::pose_from_pdb( *reference_pose_, basic::options::option[ basic::options::OptionKeys::in::file::native ] );
	}

	threshold_ = tag->getOption<core::Real>( "threshold", 5 );

	//TODO: support multiple jumps
	core::Size jump_num = tag->getOption<core::Size>( "jump", 1);

	//TODO: convert jump_num to movable_jumps_ (vector0?)
	movable_jumps_.push_back(jump_num);
}

protocols::filters::FilterOP
LRmsdFilterCreator::create_filter() const { return protocols::filters::FilterOP( new LRmsdFilter ); }

std::string
LRmsdFilterCreator::keyname() const { return "LRmsd"; }


} // filters
} // protein_interface_design
} // devel


