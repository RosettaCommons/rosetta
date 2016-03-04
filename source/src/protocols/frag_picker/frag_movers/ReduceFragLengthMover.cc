// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/frag_picker/frag_movers/ReduceFragLengthMover
///
/// @brief      Converts a larger frag set such as 9mer to a smaller one such as 3mers
/// @details
///
/// @author     TJ Brunette (tjbrunette@gmail.com)
/// @note
///
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>


// Core Headers
#include <protocols/frag_picker/frag_movers/ReduceFragLengthMover.hh>
#include <protocols/frag_picker/frag_movers/ReduceFragLengthMoverCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.frag_picker.frag_movers.ReduceFragLengthMover" );

namespace protocols {
namespace frag_picker {
namespace frag_movers {

std::string ReduceFragLengthMoverCreator::keyname() const
{
	return ReduceFragLengthMoverCreator::mover_name();
}

std::string ReduceFragLengthMoverCreator::mover_name(){
	return "ReduceFragLengthMover";
}

protocols::moves::MoverOP
ReduceFragLengthMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ReduceFragLengthMover );
}

ReduceFragLengthMover::ReduceFragLengthMover():moves::Mover("ReduceFragLengthMover"){}

void ReduceFragLengthMover::apply(core::pose::Pose & pose) {
	TR << pose.total_residue() << std::endl;
}

std::string ReduceFragLengthMover::get_name() const {
	return "ReduceFragLengthMover";
}

// void
// ReduceFragLengthMover::parse_my_tag(
//  utility::tag::TagCOP tag,
//  basic::datacache::DataMap & datamap,
//  protocols::filters::Filters_map const &,
//  protocols::moves::Movers_map const &,
//  core::pose::Pose const & ){
// }
void
ReduceFragLengthMover::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & ){
}

}//frag_movers
}//frag_picker
}//protocols
