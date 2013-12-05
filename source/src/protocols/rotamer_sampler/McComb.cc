// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/McComb.cc
/// @brief Ensemble of Markov chain samplers for sampling combinatorially.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/rotamer_sampler/McComb.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.rotamer_sampler.McComb" );

using namespace core;

namespace protocols {
namespace rotamer_sampler {
///////////////////////////////////////////////////////////////////////////
McComb::McComb():
	McRotamer()
{}
///////////////////////////////////////////////////////////////////////////
void McComb::init() {
	runtime_assert( !rotamer_list_.empty() );
	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) {
		rotamer_list_[i]->init();
		if ( !rotamer_list_[i]->not_end() )
				TR << "Got a null rotamer sampler!" << std::endl;
	}
	set_init( true );
	reset();
}
///////////////////////////////////////////////////////////////////////////
void McComb::reset() {
	runtime_assert( is_init() );
	for ( Size i = 1; i <= rotamer_list_.size(); ++i )
			rotamer_list_[i]->reset();
}
///////////////////////////////////////////////////////////////////////////
void McComb::operator++() {
	runtime_assert( not_end() );
	for ( Size i = 1; i <= rotamer_list_.size(); ++i ) ++( *rotamer_list_[i] );
}
///////////////////////////////////////////////////////////////////////////
void McComb::update() {
	runtime_assert( is_init() );
	for ( Size i = 1; i <= rotamer_list_.size(); ++i )
			rotamer_list_[i]->update();
}
///////////////////////////////////////////////////////////////////////////
void McComb::apply( Pose & pose ) {
	runtime_assert( is_init() );
	for ( Size i = 1; i <= rotamer_list_.size(); ++i )
			rotamer_list_[i]->apply( pose );
}
///////////////////////////////////////////////////////////////////////////
void McComb::set_uniform_sampling( bool const setting ) {
	McRotamer::set_uniform_sampling( setting );
	for ( Size i = 1; i <= rotamer_list_.size(); ++i )
			rotamer_list_[i]->set_uniform_sampling( setting );
}
///////////////////////////////////////////////////////////////////////////
}
}
