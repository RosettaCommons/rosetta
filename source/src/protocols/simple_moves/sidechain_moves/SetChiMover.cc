// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/sidechain_moves/SetChiMover.cc
/// @brief  A mover to change one chi angle
/// @author Noah Ollikanen

// Unit headers
#include <protocols/simple_moves/sidechain_moves/SetChiMover.hh>
#include <protocols/simple_moves/sidechain_moves/SetChiMoverCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
//parsing
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace simple_moves {
namespace sidechain_moves {

using namespace core;
using namespace core::chemical;
using namespace std;

using core::pose::Pose;
using core::conformation::Residue;

static thread_local basic::Tracer TR( "protocols.simple_moves.sidechain_moves.SetChiMover" );

std::string
SetChiMoverCreator::keyname() const
{
	return SetChiMoverCreator::mover_name();
}

protocols::moves::MoverOP
SetChiMoverCreator::create_mover() const {
	return new SetChiMover;
}

std::string
SetChiMoverCreator::mover_name()
{
	return "SetChiMover";
}

SetChiMover::~SetChiMover() {}

///@brief default ctor
SetChiMover::SetChiMover() :
	parent(),
	angle_( 0 ),
	resnum_( 0 ),
	chinum_( 0 )
{}

void SetChiMover::apply( Pose & pose ) {
	runtime_assert( resnum() > 0 );
	runtime_assert( resnum() <= pose.total_residue() );

	if ( chinum() <= pose.residue(resnum()).nchi() ) {
		pose.set_chi(chinum(), resnum(), angle());
		TR<<"Set chi"<<chinum()<<" of residue "<<resnum()<<" to "<<angle()<<std::endl;
	}
	
	pose.update_residue_neighbors();
}

std::string
SetChiMover::get_name() const {
	return SetChiMoverCreator::mover_name();
}

void SetChiMover::parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & pose)
{
	angle( tag->getOption< core::Real >( "angle" ) );
	resnum( core::pose::parse_resnum( tag->getOption< std::string >( "resnum" ), pose ) );
	chinum( tag->getOption< core::Size >( "chinum" ) );
	
}


} // sidechain_moves
} // simple_moves
} // protocols
