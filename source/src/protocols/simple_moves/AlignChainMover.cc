// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AlignChainMover.cc
/// @brief

// Unit headers
#include <protocols/simple_moves/AlignChainMover.hh>
#include <protocols/simple_moves/AlignChainMoverCreator.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.AlignChainMover" );

#include <utility/tag/Tag.hh>
#include <boost/foreach.hpp>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>
#include <numeric/xyzVector.hh>
#include <protocols/toolbox/superimpose.hh>

namespace protocols {
namespace simple_moves {

std::string
AlignChainMoverCreator::keyname() const
{
	return AlignChainMoverCreator::mover_name();
}

protocols::moves::MoverOP
AlignChainMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AlignChainMover );
}

std::string
AlignChainMoverCreator::mover_name()
{
	return "AlignChain";
}

AlignChainMover::AlignChainMover()
: moves::Mover("AlignChain"),
	pose_( /* NULL */ ),
	source_chain_( 0 ),
	target_chain_( 0 )
{
}

core::pose::PoseOP
AlignChainMover::pose() const{ return pose_; }

void
AlignChainMover::pose( core::pose::PoseOP pose ){ pose_ = pose; }

utility::vector1< numeric::xyzVector< core::Real > >
Ca_coord( core::pose::Pose const & pose, utility::vector1< core::Size > const positions ){
	utility::vector1< numeric::xyzVector< core::Real > > coords;

	coords.clear();
	BOOST_FOREACH ( core::Size const pos, positions ) {
		coords.push_back( pose.residue( pos ).xyz( "CA" ) );
	}
	return coords;
}

void
AlignChainMover::apply( Pose & in_pose )
{
	utility::vector1< core::Size > in_pose_positions, target_positions;
	in_pose_positions.clear(); target_positions.clear();
	runtime_assert( pose() != 0 );
	core::Size const in_pose_chain_begin( source_chain() == 0 ? 1 : in_pose.conformation().chain_begin( source_chain() ) );
	core::Size const in_pose_chain_end  ( source_chain() == 0 ? in_pose.total_residue() : in_pose.conformation().chain_end( source_chain() ) );
	core::Size const target_pose_chain_begin( target_chain() == 0 ? 1 : pose()->conformation().chain_begin( target_chain() ) );
	core::Size const target_pose_chain_end  ( target_chain() == 0 ? pose()->total_residue() : pose()->conformation().chain_end( target_chain() ) );
	TR<<"In_pose from residue: "<<in_pose_chain_begin<<" to_residue: "<<in_pose_chain_end<<"\ntarget_pose from residue: "<<target_pose_chain_begin<<" to_residue: "<<target_pose_chain_end<<std::endl;
	for ( core::Size i = in_pose_chain_begin; i<=in_pose_chain_end; ++i ) {
		in_pose_positions.push_back( i );
	}
	for ( core::Size i = target_pose_chain_begin; i<=target_pose_chain_end; ++i ) {
		target_positions.push_back( i );
	}

	utility::vector1< numeric::xyzVector< core::Real > > init_coords( Ca_coord( in_pose, in_pose_positions ) ), ref_coords( Ca_coord( *pose(), target_positions ) );

	numeric::xyzMatrix< core::Real > rotation;
	numeric::xyzVector< core::Real > to_init_center, to_fit_center;

	using namespace protocols::toolbox;

	superposition_transform( init_coords, ref_coords, rotation, to_init_center, to_fit_center );

	apply_superposition_transform( in_pose, rotation, to_init_center, to_fit_center );
}

std::string
AlignChainMover::get_name() const {
	return AlignChainMoverCreator::mover_name();
}

moves::MoverOP
AlignChainMover::clone() const
{
	return moves::MoverOP( new AlignChainMover( *this ) );
}

moves::MoverOP
AlignChainMover::fresh_instance() const
{
	return moves::MoverOP( new AlignChainMover );
}

void
AlignChainMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	source_chain( tag->getOption< core::Size >( "source_chain", 0 ) );
	target_chain( tag->getOption< core::Size >( "target_chain", 0 ) );
	std::string const fname( tag->getOption< std::string >( "target_name" ) );

	pose( core::import_pose::pose_from_file( fname ) );

	TR<<"source_chain: "<<source_chain()<<" target_chain: "<<target_chain()<<" pdb name: "<<fname<<std::endl;
}

} // simple_moves
} // protocols
