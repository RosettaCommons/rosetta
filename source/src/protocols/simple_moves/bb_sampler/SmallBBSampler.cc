// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/simple_moves/bb_sampler/SmallBBSampler.cc
/// @brief A bb sampler that samples within a range of a starting angle.  Similar to small mover.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/simple_moves/bb_sampler/SmallBBSampler.hh>
#include <protocols/simple_moves/bb_sampler/BBDihedralSampler.hh>

#include <core/pose/util.hh>

#include <numeric/random/random.hh>

#include <basic/Tracer.hh>
#include <basic/basic.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.bb_sampler.SmallBBSampler" );


namespace protocols {
namespace simple_moves {
namespace bb_sampler {

SmallBBSampler::SmallBBSampler():
	BBDihedralSampler()
{
	set_angle_max( 360 );
}

SmallBBSampler::SmallBBSampler( core::id::MainchainTorsionType torsion_type):
	BBDihedralSampler()
{
	torsion_type_ = torsion_type;
	set_angle_max( 360 );
}

SmallBBSampler::SmallBBSampler( core::id::MainchainTorsionType torsion_type, core::Real angle_max):
	BBDihedralSampler()

{
	torsion_type_ = torsion_type;
	set_angle_max( angle_max );
}

void SmallBBSampler::set_angle_max( core::Real const angle )
{
	set_angle_max( 'H', angle ); // helix
	set_angle_max( 'L', angle ); // other
	set_angle_max( 'E', angle ); // strand
}

void SmallBBSampler::set_angle_max( char const type, core::Real const angle ) {
	angle_max_[ type ] = angle;
}

void SmallBBSampler::set_angle_max( std::map< char, core::Real > angle_max_in ) {
	angle_max_.swap( angle_max_in );
}

core::Real
SmallBBSampler::get_angle_max(char const type) const
{
	return angle_max_.find(type)->second;
}


SmallBBSampler::~SmallBBSampler(){}

SmallBBSampler::SmallBBSampler( SmallBBSampler const & src ):
	BBDihedralSampler(src),
	angle_max_(src.angle_max_)
{

}

core::Real
SmallBBSampler::get_torsion(core::pose::Pose const & pose, core::Size resnum) const {

	using namespace core::pose;

	char const ss( pose.secstruct( resnum ) );
	//core::Angle const mx( angle_max_.find( ss )->second );
	core::Angle current_angle = get_bb_torsion( core::Size( torsion_type_ ), pose, resnum);
	//TR << "current: " << current_angle;

	core::Angle max_dev = get_angle_max( ss )/2;

	core::Angle new_angle = basic::periodic_range( current_angle - max_dev + numeric::random::rg().uniform() * max_dev * 2, 360.0 );
	//TR << "new angle: " << new_angle;

	return new_angle;

}

void
SmallBBSampler::set_torsion_to_pose(core::pose::Pose &pose, core::Size resnum) const {

	core::Angle new_angle = get_torsion( pose, resnum );
	core::pose::set_bb_torsion( core::Size( torsion_type_ ), pose, resnum, new_angle);

}

SmallBBSamplerOP
SmallBBSampler::clone() const {
	return SmallBBSamplerOP( new SmallBBSampler( *this ) );
}


} //protocols
} //simple_moves
} //bb_sampler






