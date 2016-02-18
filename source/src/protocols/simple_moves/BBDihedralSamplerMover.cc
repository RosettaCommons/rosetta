// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/carbohydrates/BBDihedralSamplerMover.cc
/// @brief Mover interface to BBDihedralSampler.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/simple_moves/BBDihedralSamplerMover.hh>
#include <protocols/simple_moves/BBDihedralSamplerMoverCreator.hh>
#include <protocols/simple_moves/bb_sampler/BBDihedralSampler.hh>

#include <core/pose/Pose.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh>

#include <core/id/types.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.BBDihedralSamplerMover" );


namespace protocols {
namespace simple_moves {
using namespace protocols::simple_moves::bb_sampler;

BBDihedralSamplerMover::BBDihedralSamplerMover():
	protocols::moves::Mover( "BBDihedralSamplerMover" ),
	sampler_(/* NULL */)
{

}

BBDihedralSamplerMover::BBDihedralSamplerMover( BBDihedralSamplerOP sampler):
	protocols::moves::Mover( "BBDihedralSamplerMover" )

{
	sampler_ = sampler;
}

BBDihedralSamplerMover::BBDihedralSamplerMover( BBDihedralSamplerOP sampler, core::kinematics::MoveMapCOP movemap):
	protocols::moves::Mover( "BBDihedralSamplerMover" )
{
	sampler_ = sampler;
	set_movemap( movemap );
}

BBDihedralSamplerMover::~BBDihedralSamplerMover(){}

BBDihedralSamplerMover::BBDihedralSamplerMover( BBDihedralSamplerMover const & src ):
	protocols::moves::Mover( src ),
	sampler_(src.sampler_),
	bb_residues_(src.bb_residues_)
{

}

void
BBDihedralSamplerMover::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

}

protocols::moves::MoverOP
BBDihedralSamplerMover::clone() const{
	return protocols::moves::MoverOP( new BBDihedralSamplerMover( *this ) );
}

/*
BBDihedralSamplerMover & BBDihedralSamplerMoveroperator=( BBDihedralSamplerMover const & src){
return BBDihedralSamplerMover( src );
}
*/


moves::MoverOP
BBDihedralSamplerMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new BBDihedralSamplerMover );
}

std::string
BBDihedralSamplerMover::get_name() const {
	return "BBDihedralSamplerMover";
}

void
BBDihedralSamplerMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

std::ostream &operator<< (std::ostream &os, BBDihedralSamplerMover const &mover)
{
	mover.show(os);
	return os;
}

void
BBDihedralSamplerMover::set_movemap( core::kinematics::MoveMapCOP movemap){
	using namespace core::kinematics;
	bb_residues_ = get_residues_from_movemap_with_id( core::id::BB, *movemap); //This is done here so we dont have to do this at each apply and waste time.

}

void
BBDihedralSamplerMover::set_single_resnum( core::Size resnum ){
	bb_residues_.clear();
	bb_residues_.push_back( resnum );
}

void
BBDihedralSamplerMover::apply( core::pose::Pose & pose ){
	if ( ! sampler_ ) {
		utility_exit_with_message(" No Sampler set for BBDihedralSamplerMover!");
	}

	if ( bb_residues_.size() == 0 ) {
		TR << "No Movemap Set.  Using all residues." << std::endl;
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			bb_residues_.push_back( i );
		}
	}

	core::Size index = numeric::random::rg().random_range( 1, bb_residues_.size() );
	core::Size resnum = bb_residues_[ index ];

	//TR << "Optimizing "<< resnum << " with " << sampler_->name() << " torsion " << core::Size(sampler_->get_torsion_type()) << std::endl;

	try {
		sampler_->set_torsion_to_pose( pose, resnum );
		set_last_move_status(protocols::moves::MS_SUCCESS);

	} catch ( utility::excn::EXCN_Base& excn ) {
		TR.Debug << "Could not set torsion for resnum "<< resnum << std::endl;
		set_last_move_status(protocols::moves::MS_FAIL);

	}

}


/////////////// Creator ///////////////

protocols::moves::MoverOP
BBDihedralSamplerMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new BBDihedralSamplerMover );
}

std::string
BBDihedralSamplerMoverCreator::keyname() const {
	return BBDihedralSamplerMoverCreator::mover_name();
}

std::string
BBDihedralSamplerMoverCreator::mover_name(){
	return "BBDihedralSamplerMover";
}

} //protocols
} //carbohydrates


