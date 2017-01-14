// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/sampler/MC_Comb.cc
/// @brief Ensemble of Markov chain samplers for modeler combinatorially.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/recces/sampler/MC_Comb.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/id/TorsionID.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.recces.sampler.MC_Comb" );

using namespace core;

namespace protocols {
namespace recces {
namespace sampler {
///////////////////////////////////////////////////////////////////////////
MC_Comb::MC_Comb():
	MC_Sampler()
{
	set_name( "MC_Comb" );
}
///////////////////////////////////////////////////////////////////////////
void MC_Comb::init() {
	//	runtime_assert( !rotamer_list_.empty() );
	for ( auto const & rotamer : rotamer_list_ ) {
		rotamer->init();
	}
	set_init( true );
	if ( update_pose_ != 0 ) set_update_pose( update_pose_ ); // in case any rotamers got created by init()
	reset();
}
///////////////////////////////////////////////////////////////////////////
void MC_Comb::reset() {
	runtime_assert( is_init() );
	for ( auto const & rotamer : rotamer_list_ ) {
		rotamer->reset();
	}
}
///////////////////////////////////////////////////////////////////////////
void MC_Comb::operator++() {
	for ( auto const & rotamer : rotamer_list_ ) ++( *rotamer );
}
///////////////////////////////////////////////////////////////////////////
void MC_Comb::update() {
	runtime_assert( is_init() );
	for ( auto const & rotamer : rotamer_list_ ) {
		rotamer->update();
	}
}
///////////////////////////////////////////////////////////////////////////
void MC_Comb::apply( Pose & pose ) {
	runtime_assert( is_init() );
	found_move_ = false;
	for ( auto const & rotamer : rotamer_list_ ) {
		rotamer->apply( pose );
		if ( rotamer->found_move() ) found_move_ = true;
	}
}
///////////////////////////////////////////////////////////////////////////
void MC_Comb::set_uniform_modeler( bool const setting ) {
	MC_Sampler::set_uniform_modeler( setting );
	for ( auto const & rotamer : rotamer_list_ ) {
		rotamer->set_uniform_modeler( setting );
	}
}
///////////////////////////////////////////////////////////////////////////
void MC_Comb::set_update_pose( core::pose::PoseCOP setting ) {
	update_pose_ = setting;
	for ( auto const & rotamer : rotamer_list_ ) {
		rotamer->set_update_pose( setting );
	}
}
///////////////////////////////////////////////////////////////////////////
void MC_Comb::show( std::ostream & out, Size const indent ) const {
	SamplerPlusPlus::show( out, indent );
	for ( Size k = 1; k <= rotamer_list_.size(); k++ ) rotamer_list_[k]->show( out, indent + 1 );
}
///////////////////////////////////////////////////////////////////////////
MC_SamplerOP
MC_Comb::find( core::id::TorsionID const & torsion_id ) {
	for ( auto rotamer : rotamer_list_ ) {
		if ( rotamer->find( torsion_id ) != 0 ) return rotamer->find( torsion_id );
	}
	return 0;
}

} //sampler
} //recces
} //protocols
