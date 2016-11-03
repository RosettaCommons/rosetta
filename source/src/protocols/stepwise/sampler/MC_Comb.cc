// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/MC_Comb.cc
/// @brief Ensemble of Markov chain samplers for modeler combinatorially.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/stepwise/sampler/MC_Comb.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.sampler.MC_Comb" );

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampler {
///////////////////////////////////////////////////////////////////////////
MC_Comb::MC_Comb():
	MC_StepWiseSampler()
{}
///////////////////////////////////////////////////////////////////////////
void MC_Comb::init() {
	runtime_assert( !rotamer_list_.empty() );
	for ( auto const & rotamer : rotamer_list_ ) {
		rotamer->init();
		if ( !rotamer->not_end() ) {
			TR << "Got a null rotamer sampler!" << std::endl;
		}
	}
	set_init( true );
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
	runtime_assert( not_end() );
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
	for ( auto const & rotamer : rotamer_list_ ) {
		rotamer->apply( pose );
	}
}
///////////////////////////////////////////////////////////////////////////
void MC_Comb::set_uniform_modeler( bool const setting ) {
	MC_StepWiseSampler::set_uniform_modeler( setting );
	for ( auto const & rotamer : rotamer_list_ ) {
		rotamer->set_uniform_modeler( setting );
	}
}
///////////////////////////////////////////////////////////////////////////
} //sampler
} //stepwise
} //protocols
