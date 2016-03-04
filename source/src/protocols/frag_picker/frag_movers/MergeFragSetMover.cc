// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/frag_picker/frag_movers/MergeFragSetMover
///
/// @brief      Merges two fragment sets
///
/// @author     TJ Brunette (tjbrunette@gmail.com)
/// @note

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>


// Core Headers
#include <protocols/frag_picker/frag_movers/MergeFragSetMover.hh>
#include <protocols/frag_picker/frag_movers/MergeFragSetMoverCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.frag_picker.frag_movers.MergeFragSetMover" );

namespace protocols {
namespace frag_picker {
namespace frag_movers {

std::string MergeFragSetMoverCreator::keyname() const
{
	return MergeFragSetMoverCreator::mover_name();
}

std::string MergeFragSetMoverCreator::mover_name(){
	return "MergeFragSetMover";
}

protocols::moves::MoverOP
MergeFragSetMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MergeFragSetMover );
}

MergeFragSetMover::MergeFragSetMover():moves::Mover("MergeFragSetMover"){}

void MergeFragSetMover::apply(core::pose::Pose & pose) {
	TR << pose.total_residue();
}

std::string MergeFragSetMover::get_name() const {
	return "MergeFragSetMover";
}

// void
// MergeFragSetMover::parse_my_tag(
//  utility::tag::TagCOP tag,
//  basic::datacache::DataMap & datamap,
//  protocols::filters::Filters_map const &,
//  protocols::moves::Movers_map const &,
//  core::pose::Pose const & ){
// }

void MergeFragSetMover::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & ){
}

}//frag_movers
}//frag_picker
}//protocols
