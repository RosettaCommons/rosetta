// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./protocols/moves/DsspMover.cc
/// @brief performs dssp and set secondary structure to pose
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit headers
#include <protocols/moves/DsspMover.hh>
#include <protocols/moves/DsspMoverCreator.hh>

// Project Headers
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/moves/Mover.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Parser headers


#include <string>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.DsspMover" );

namespace protocols {
namespace moves {

std::string
DsspMoverCreator::keyname() const
{
	return DsspMoverCreator::mover_name();
}

protocols::moves::MoverOP
DsspMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new DsspMover );
}

std::string
DsspMoverCreator::mover_name()
{
	return "Dssp";
}

DsspMover::DsspMover():
	Mover( DsspMoverCreator::mover_name() ),
	reduced_IG_as_L_(0)
{}

DsspMover::~DsspMover()
= default;

/// @brief clone this object
DsspMover::MoverOP DsspMover::clone() const {

	return DsspMover::MoverOP( new DsspMover( *this ) );
}

/// @brief create this type of object
DsspMover::MoverOP DsspMover::fresh_instance() const {
	return DsspMover::MoverOP( new DsspMover() );
}

/// @details virtual main
void DsspMover::apply( Pose & pose )
{
	if ( !reduced_IG_as_L_ ) {
		core::scoring::dssp::Dssp dssp( pose );
		dssp.insert_ss_into_pose( pose );
		TR << dssp.get_dssp_secstruct() << std::endl;
	} else {
		TR << "reduce IG as L" << std::endl;
		core::scoring::dssp::Dssp dssp( pose );
		dssp.insert_ss_into_pose_no_IG_helix( pose );
		TR << dssp.get_dssp_secstruct() << std::endl;
	}
}

std::string
DsspMover::get_name() const {
	return DsspMoverCreator::mover_name();
}

/// @brief parse xml
void
DsspMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	reduced_IG_as_L_ = tag->getOption<bool>( "reduced_IG_as_L" , 0 );
	TR << "DsspMover loaded." << std::endl;
}

}  // namespace moves
}  // namespace protocols
