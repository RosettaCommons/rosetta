// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
using namespace core::kinematics;

BBDihedralSamplerMover::BBDihedralSamplerMover():
	protocols::moves::Mover( "BBDihedralSamplerMover" ),
	movemap_(/* NULL */)
{

}

BBDihedralSamplerMover::BBDihedralSamplerMover( BBDihedralSamplerOP sampler):
	protocols::moves::Mover( "BBDihedralSamplerMover" ),
	movemap_(/* NULL */)

{
	set_sampler( sampler );
}

BBDihedralSamplerMover::BBDihedralSamplerMover( BBDihedralSamplerOP sampler, core::kinematics::MoveMapCOP movemap):
	protocols::moves::Mover( "BBDihedralSamplerMover" )
{
	set_sampler( sampler );
	set_movemap( movemap );
}

BBDihedralSamplerMover::~BBDihedralSamplerMover()= default;

BBDihedralSamplerMover::BBDihedralSamplerMover( BBDihedralSamplerMover const & src ):
	protocols::moves::Mover( src ),
	sampler_torsion_types_(src.sampler_torsion_types_),
	bb_residues_(src.bb_residues_),
	sampler_movemap_union_(src.sampler_movemap_union_)
{
	samplers_.clear();
	for ( auto const & kv : src.samplers_ ) {
		utility::vector1< BBDihedralSamplerOP > sampler_list;
		for ( BBDihedralSamplerOP sampler : kv.second ) {
			BBDihedralSamplerOP new_sampler = sampler->clone();
			sampler_list.push_back( new_sampler );
		}
		samplers_[kv.first] = sampler_list;
	}
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
	movemap_ = movemap;
	sampler_movemap_union_.clear();
	bb_residues_.clear();

}

void
BBDihedralSamplerMover::set_single_resnum( core::Size resnum ){
	bb_residues_.clear();
	sampler_movemap_union_.clear();
	bb_residues_.push_back( resnum );
}

void
BBDihedralSamplerMover::set_sampler( bb_sampler::BBDihedralSamplerOP sampler ){
	samplers_.clear();
	sampler_torsion_types_.clear();
	sampler_movemap_union_.clear();
	samplers_[ sampler->get_torsion_type() ]; //Initialize vector
	samplers_[ sampler->get_torsion_type() ].push_back( sampler );


	sampler_torsion_types_.push_back( sampler->get_torsion_type() );

}

void
BBDihedralSamplerMover::add_sampler( bb_sampler::BBDihedralSamplerOP sampler ){
	sampler_movemap_union_.clear();
	if ( samplers_.count( sampler->get_torsion_type() ) != 0 ) {
		samplers_[ sampler->get_torsion_type() ].push_back( sampler );
		//sampler_torsion_types_.push_back( sampler->get_torsion_type() );
	} else {
		samplers_[ sampler->get_torsion_type() ];
		samplers_[ sampler->get_torsion_type() ].push_back( sampler );

		sampler_torsion_types_.push_back( sampler->get_torsion_type() );
	}

}

void
BBDihedralSamplerMover::setup_all_bb_residues( core::pose::Pose const & pose) {

	if ( TR.Debug.visible() ) {
		TR.Debug << "Setting up all bb residues" << std::endl;
	}

	MoveMapOP mm = MoveMapOP( new MoveMap );

	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		mm->set_bb(i, true);
		bb_residues_.push_back( i );

	}
	movemap_ = mm;

}

void
BBDihedralSamplerMover::setup_sampler_movemap_union( core::pose::Pose const & pose ) {

	utility::vector1< core::Size > pruned_resnums;
	if ( TR.Debug.visible() ) {
		TR.Debug << "Initializing sampler and movemap union." << std::endl;
	}

	if ( bb_residues_.size() == 0 ) {
		setup_all_bb_residues( pose );
	}

	for ( core::Size i = 1; i <= bb_residues_.size(); ++i ) {
		core::Size resnum = bb_residues_[ i ];
		sampler_movemap_union_[ resnum ];
		for ( core::Size x = 1; x <= sampler_torsion_types_.size(); ++x ) {
			//TR << "Resnum " << resnum << " Torsion " << sampler_torsion_types_[x] << " ON? " << movemap_->get_bb( resnum, sampler_torsion_types_[x]) << std::endl;
			if ( movemap_->get_bb( resnum, sampler_torsion_types_[x]) ) {
				sampler_movemap_union_[ resnum ].push_back( sampler_torsion_types_[x] );
			}
		}
	}

	//Update BB Residues to only those residues which both have torsions set on in the movemap and we have samplers to actually sample those torsions.
	typedef std::map<core::Size, utility::vector1< core::Size > >::const_iterator iter_type;
	for ( iter_type it = sampler_movemap_union_.begin(); it != sampler_movemap_union_.end(); ++it ) {
		if ( it->second.size() > 0 ) {
			pruned_resnums.push_back( it->first );
		}
	}
	bb_residues_ = pruned_resnums;
}

void
BBDihedralSamplerMover::apply( core::pose::Pose & pose ){
	if ( movemap_ && bb_residues_.size() == 0 ) {
		bb_residues_ = get_residues_from_movemap_bb_any_torsion( *movemap_, pose.size() ); //This is done here so we dont have to do this at each apply and waste time.
	}

	if ( samplers_.size() == 0 ) {
		utility_exit_with_message(" No Sampler set for BBDihedralSamplerMover!");
	}

	if ( bb_residues_.size() == 0 ) {
		setup_all_bb_residues( pose );
	}

	if ( sampler_movemap_union_.empty() ) {
		setup_sampler_movemap_union( pose );
	}


	core::Size index = numeric::random::rg().random_range( 1, bb_residues_.size() );
	core::Size resnum = bb_residues_[ index ];



	//Get the one sampler or choose a sampler:
	bb_sampler::BBDihedralSamplerCOP sampler;
	if ( sampler_torsion_types_.size() == 1 && samplers_[ sampler_torsion_types_[ 1 ] ].size() == 1 ) {
		//Check to make sure this is in the movemap.
		if ( ! movemap_->get_bb( resnum, sampler_torsion_types_[ 1 ] ) ) {
			TR << "Resnum set as true, but torsion set to false in Movemap.  Cannot sample. "<< resnum <<","<< sampler_torsion_types_[ 1 ] << std::endl;
			set_last_move_status(protocols::moves::MS_FAIL);
			return;
		} else {
			sampler = samplers_[ sampler_torsion_types_[ 1 ] ][ 1 ]; //First and only torsion type, first and only sampler that samples that torsion type.
		}

	} else if ( sampler_movemap_union_[ resnum ].size() == 0 ) {
		utility_exit_with_message(" BBDihedralSamplerMover - a chosen residue num has no dihedral union between movemaps and set bb samplers.  We should never be here!");
	} else {

		//Choose a TorsionType based on what kind of samplers we have and the MoveMap of the residue.
		core::Size torsion_index = numeric::random::rg().random_range( 1, sampler_movemap_union_[ resnum ].size() );
		core::Size torsion = sampler_movemap_union_[ resnum ][ torsion_index ];


		//Choose a sampler from the samplers available for that torsion type.
		core::Size sampler_index = numeric::random::rg().random_range( 1, samplers_[ torsion ].size() );
		sampler = samplers_[ torsion ][ sampler_index ];
	}


	// Apply the sampler

	TR<< "Optimizing "<< resnum << " with " << sampler->get_name() << " at torsion " << core::Size(sampler->get_torsion_type()) << std::endl;


	try {
		sampler->set_torsion_to_pose( pose, resnum );
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


